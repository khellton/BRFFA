standardiseVar <- function(x)
{
  return((x-mean(x,na.rm=TRUE))/sqrt(var(x,na.rm=TRUE)))
} 

returnLevel <- function(MU,KA,XI,T,useGumbel=FALSE) {
  p <- 1-1/T
  if (!useGumbel) { # GEV
    tmp <-  (-log(p))^{-XI} - 1
    y.p <- MU + tmp / (XI*KA)
  }
  else {
    y.p <- MU - log( -log(p) )/ (KA) # GUMBEL   evd::qgumbel(1-1/1000,loc=MU.0,scale=1/KA.0)
  }
  return(y.p)
}

plotMap.MAP <- function(MAP,long,lat,k=1,
                        all.colors=NULL)
{
  if (is.null(all.colors)) {
    all.colors <- heat.colors(7)[7:1]
  }
  
  # mu
  mugr <- rep(NA,length(MAP[,"mu"]))
  gr <- round(seq(min(MAP[,"mu"]),max(MAP[,"mu"]),length.out=7),2)
  gr <- c(0,0.1,0.2,0.3,0.5,0.7,1)
  for (g in 1:(length(gr)-1)){
    mugr[MAP[,"mu"]>=gr[g] & MAP[,"mu"]<gr[g+1]] <- g
  }
  ltxt <-  paste(gr[-length(gr)],"-",gr[-1],sep="") # text for legend
  par(mar=c(1,1,1,1)+0.01)  
  map("world",regions="Norway", fill=TRUE,col="gray80",ylim=c(58,71),xlim=c(5,30.5),border=FALSE)
  box()
  points(long, lat, pch=16, col=all.colors[mugr])
  legend("bottomright",legend=ltxt,fill=,col=all.colors[1:7],pch=16)
  legend("topleft","mu")
  #qplot(long, lat)
  par(new=TRUE)
  xc <- k
  yc <- 1
  par(mfg=c(xc,yc))
  par(mar=c(11,15,8,3))
  hist(MAP[,"mu"],breaks=50,xlab="",ylab="",main="",axes=FALSE)
  axis(1,cex.axis=0.5,line=-1.8,tick=FALSE)   
  
  #kappa
  kappagr <- rep(NA,length(MAP[,"kappa"]))
  gr <- c(1,10,12,14,16,50,150)
  for (g in 1:(length(gr)-1)){
    kappagr[MAP[,"kappa"]>=gr[g] & MAP[,"kappa"]<gr[g+1]] <- g
  }
  ltxt <-  paste(gr[-length(gr)],"-",gr[-1],sep="") # text for legend
  par(mar=c(1,1,1,1)+0.01)
  map("world",regions="Norway", fill=TRUE,col="gray80",ylim=c(58,71),xlim=c(5,30.5),border=FALSE)
  box()
  points(long, lat, pch=16, col=all.colors[kappagr])
  legend("bottomright",legend=ltxt,fill=,col=all.colors[1:7],pch=16)
  legend("topleft","kappa")
  #qplot(long, lat)
  par(new=TRUE)
  xc <- k
  yc <- 2
  par(mfg=c(xc,yc))
  par(mar=c(11,15,8,3))
  hist(MAP[,"kappa"],breaks=50,xlab="",ylab="",main="",axes=FALSE)
  axis(1,cex.axis=0.5,line=-1.8,tick=FALSE)
  
  #xi
  xigr <- rep(NA,length(MAP[,"xi"]))
  gr <- c(-0.6,-0.2,-0.1,0,0.1,0.2,0.6)
  for (g in 1:(length(gr)-1)){
    xigr[MAP[,"xi"]>=gr[g] & MAP[,"xi"]<gr[g+1]] <- g
  }
  ltxt <-  paste(gr[-length(gr)],"-",gr[-1],sep="") # text for legend
  par(mar=c(1,1,1,1)+0.01)  
  map("world",regions="Norway", fill=TRUE,col="gray80",ylim=c(58,71),xlim=c(5,30.5),border=FALSE)
  box()
  points(long, lat, pch=16, col=all.colors[xigr])
  legend("bottomright",legend=ltxt,fill=,col=all.colors[1:7],pch=16)
  legend("topleft","xi")
  #qplot(long, lat)  
  par(new=TRUE)
  xc <- k
  yc <- 3
  par(mfg=c(xc,yc))
  par(mar=c(11,15,8,3))
  hist(MAP[,"xi"],breaks=50,xlab="",ylab="",main="",axes=FALSE)
  axis(1,cex.axis=0.5,line=-1.8,tick=FALSE)
}

plotReturnLevel <- function(return.level = NULL,
                            return.level.local = NULL,
                            return.level.old = NULL,
                            return.level.pred = NULL,
                            floodData,floodCov,
                            T = NULL,
                            T.old = NULL,
                            writeToFile=FALSE,
                            plotDir="plot2016/",
                            fname = NULL,
                            ymin=NULL,ymax=NULL,
                            useLiter=TRUE)
{
  x.lab <- c(1,2,3,5,10,20,30,50,100,200,300,500,1000)
  if (writeToFile) {
    if(is.null(fname)) {
      fname <- paste(plotDir,"returnLevelCV_",format(Sys.time(), "%d%b%Y"),".pdf",sep="")
    }
    print(fname)
    pdf(paste(plotDir,fname),width=6.5,height=6.5)
  }
  
  stno <- names(return.level)
  # find the cooresponding data
  id <- stno
  print(id)
  tmpfloodData <- floodData[id,,drop= FALSE]
  tmpfloodCov <- floodCov[id,,drop= FALSE]
  stationnm <- tmpfloodCov[,"station.name"]
  names(stationnm) <- row.names(tmpfloodCov)
  
  delta <- 1
  ytxt <- "Return level (m3/s/km2)"
  if (useLiter) {
    delta <- 1000
    ytxt <- "Return level (l/s/km2)"
    if (!is.null(ymin))
      ymin <- ymin*1000
    if (!is.null(ymax))
      ymax_temp <- ymax*1000    
  }
  lty.old <- 3
  lty.pred <- 2 
  lty.loc <- 1
  
  nrow <- 1
  ncol <- 1
  par(mfrow=c(nrow,ncol))
  par(mar=c(2.0,1.4,1.2,0.8),oma=c(1.5,1.5,0,0),mgp=c(0,0.05,0))
  
  for (i in as.character(stno)) {
    
    tmpdata <- tmpfloodData[i,]
    tmpdata <- tmpdata[!is.na(tmpdata)]
    tmpdata <- tmpdata*delta
    N <- length(tmpdata)
    AEP <- ((1:N)-0.44)/(N+0.12)
    
    up <- return.level[[i]]$upper*delta
    lo <- return.level[[i]]$lower*delta
    mm <- return.level[[i]]$MEAN*delta
    
    if (is.null(ymax))
      ymax_temp <- max(up)
    if (is.null(ymin))      
      ymin <- min(lo)
    plot(log(T),mm,type="n",ylim=c(ymin,ymax_temp),xlim=c(0,log(max(T))),
         ylab="",xlab="",axes=FALSE)
    # set axis
    magnify <- 1.2
    mtext(ytxt,side=2,line=0.9,cex=1*magnify)
    mtext("Return period (years)",side=1,line=0.9,cex=1*magnify)
    axis(1, at=log(x.lab), label=x.lab,font=0.5,tck=-0.008,cex.axis=0.8*magnify)
    abline(v=log(x.lab),lty=3,col="gray86") #line=-0.5,font=0.5,tck=0)
    axis(2,font=0.5,tck=-0.008,cex.axis=0.8*magnify)
    box()
    
    ############# regional out-of-sample
    polygon(c(rev(log(T)), log(T)), c(rev(lo), up), col = 'grey85', border = NA)
    lines(log(T),mm,lty=lty.pred,lwd=1.2)
    
    ############# local analysis
    if(!is.null(return.level.local)){
      mm <- return.level.local[[i]]$MEAN*delta     
      lines(log(T),mm,lty=lty.loc,lwd=1.2)
    }
    
    ############# old model
    if(!is.null(return.level.old)){
      mm <- return.level.old[i,]*delta
      lines(log(T.old),mm,lty=lty.old,lwd=1.2)
      points(log(T.old),mm,pch=18)
    }
    
    ############# regional in-sample
    if(!is.null(return.level.pred)){
      mm <- return.level.pred[[i]]*delta    
      lines(log(T),mm,lty=lty.pred,lwd=1.2)
    }
    
    # plot data
    points(log(1/(1-AEP)),sort(tmpdata),pch=20)
    mtext(paste(stationnm[i],i,";",length(tmpdata),"years"),side=3,cex=1*magnify,line=0)
    
  }
  if (writeToFile)
    dev.off()  
}

