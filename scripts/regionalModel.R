####################################################
##### Bayesian Regional Flood Frequency Analysis
####################################################

# Load yearly flood time series and covariates
load('data/norwayFlood.RData')

# Load functions 
source('scripts/functions.R')

# Define location for map plotting
S <- cbind(covariates$longitude,covariates$latitude)

#Plot histograms for some covariates 
par(mar=c(5,4,1,1),mfrow=c(2,3),cex.axis=1.3,cex.lab=1.5)
hist(rowSums(!is.na(floodData)),xlab='Record length (years)',ylab='Stations',main = '')
hist(log(covariates$area.total),xlab='Catchment area (km²)',ylab='',main = '',xaxt='n')
axis(1,at=log(c(40,10^(2:6))),labels=c(40,10^(2:6)))
hist(covariates$pct.eff.lake,xlab='Lake percentage',ylab='',main = '')
hist(covariates$avg.rain.april,xlab='Average rain April (mm)',ylab='Stations',main = '')
hist(covariates$avg.rain.august,xlab='Mean annual precipitation (mm)',ylab='',main = '')
hist(covariates$avg.frac.rain,xlab='Rain contribution (%)',ylab='',main = '')

#Plot map for some covariates 
library(maps)
all.colors <- colorRampPalette(c("yellow","forestgreen", "blue"))(7) #heat.colors(7)[7:1]
par(mar=c(5,4,1,1),mfrow=c(1,1),cex.axis=1.3,cex.lab=1.5)
var1 <- covariates$avg.frac.rain
mugr <- rep(NA,length(var1))
gr <- seq(0, 1, length.out=6)
#gr <- c(0,0.1,0.2,0.3,0.5,0.7,1)
for (g in 1:(length(gr)-1)){
  mugr[var1>=gr[g] & var1<gr[g+1]] <- g
}
gr <- 100*gr
ltxt <-  paste0("[",format(gr[-length(gr)]),",",format(gr[-1]),']') # text for legend

par(mar=c(0,0,3,0))  
maps::map("world",regions="Norway",fill=TRUE,col="gray80",ylim=c(58,71),xlim=c(5,30.5),border=FALSE)
title("Rain contribution (%)")
box()
points(S[,1], S[,2], pch=16, col=all.colors[mugr])
legend("bottomright",legend=ltxt,fill=,col=all.colors[1:7],pch=16,cex=1.4)

####
# set the response variable
####
#datalist <- as.list(as.data.frame(t(floodData/covariates$area_total)))
datalist <- as.list(as.data.frame(t(floodData)))
# remove all NA's from datalist
for (i in 1:length(datalist)) {
  datalist[[i]] <- (datalist[[i]])[!is.na(datalist[[i]])]
}

#Construct matrix of covariates
covMat <- as.matrix(covariates[,2:13])

# Standardize covariates
for (i in 2:dim(covMat)[2]) {
  covMat[,i] <- standardiseVar(covMat[,i])
}

##########################################
#### Run model
##########################################
nsim <- 100

# Specify prior
prior.user <- NULL
prior.user$eta$beta.0 <- c(log(1000),rep(0,dim(covMat)[2]-1))

library(SpatGEVBMA)

starttime <- Sys.time()
res <- spatial.gev.bma(datalist,covMat,S,nsim,prior.user=prior.user,nonspatial = TRUE, log.kappa=TRUE,print.every=5000)
endtime <- Sys.time()
res$cov.names <- colnames(covariates)

#Process MCMC samples and remove burn-in
burn <- 0.2*nsim
ch <- gev.process.results(res,burn=burn)

#Print inclusion probabilities for covariates
mat <- cbind(mu=ch$tbl.mu[,"Prob"],kappa=ch$tbl.eta[,"Prob"],xi=ch$tbl.xi[,"Prob"])
rownames(mat)<- c("Constant",colnames(covMat)[2:dim(covMat)[2]])
round(mat,2)*100

##########################################
###### #Plot return values plots
##########################################

#Selection of station
st.no <- rownames(floodData)[1:2]
T <- c(seq(1.05,1.9,by=0.05),2,3,5,10,20,50,100,200,500,1000)
return.level <- list()
len <- dim(res$THETA[-(1:burn),,])[1] #len=nsim-burn
zp <- matrix(0,nrow=len,ncol=length(T))
p <- 1/T

# Calculated GEV distribution based on MCMC samples
for (i in st.no) {
  mu <- res$THETA[-(1:burn),,1]%*%covMat[i,] + rnorm(len,0,sd = 1/sqrt(res$ALPHA[-(1:burn),1]))
  ka <- res$THETA[-(1:burn),,2]%*%covMat[i,] + rnorm(len,0,sd = 1/sqrt(res$ALPHA[-(1:burn),2]))
  xi <- res$THETA[-(1:burn),,3]%*%covMat[i,] + rnorm(len,0,sd = 1/sqrt(res$ALPHA[-(1:burn),3]))
  if(res$log.kappa) ka <- exp(ka)
  for(j in 1:length(T)) {
    zp[,j] <- returnLevel(mu,ka,xi,T[j],useGumbel=FALSE)
  }
  
  tmplist <- list()
  tmplist$MEAN <- apply(zp,2,median)
  tmplist$lower <- apply(zp,2,quantile,0.1)
  tmplist$upper <- apply(zp,2,quantile,0.9)  
  return.level[[i]] <- tmplist    
}

#Plot return values plots
plotReturnLevel(return.level = return.level,floodData = floodData,
                floodCov = covariates,
                T=T,writeToFile=FALSE,ymin=-0.01,ymax=NULL)


