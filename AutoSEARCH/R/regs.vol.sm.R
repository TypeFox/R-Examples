regs.vol.sm <-
function(e, vc=TRUE, arch=NULL, asym=NULL,
  log.ewma=NULL, vx=NULL, p=2, zero.adj=0.1)
{
#logep:
logep <- gLog.ep(e, zero.adj=zero.adj, p=p)
mX <- cbind(logep)
colnames(mX) <- if(p==2){"loge2"}else{"logep"}

#variance constant:
e.n <- length(e)
if(!is.null(vc)){vconst <- rep(1,e.n)}else{vconst <- NULL}
mX <- cbind(mX,vconst)

#create exponential arch terms:
archlags <- NULL
if(is.null(arch)){ archadj<-NULL ; s2<-1 }else{
  archadj <- setdiff(round(abs(arch)),0)
  for(i in 1:length(archadj))
    {
    archlags <- cbind(archlags,gLag(logep,k=archadj[i]))
    }
  colnames(archlags) <- paste(c("arch"),archadj,sep="")
  s2 <- I(max(archadj)+1)
}
mX <- cbind(mX, archlags)

#create asymmetry terms:
asymlags <- NULL
if(is.null(asym)){asymadj <- NULL}else
  {
  asymadj <- setdiff(round(abs(asym)),0)
  for(i in 1:length(asymadj))
    {
    asymlags <- cbind(asymlags,gLag(logep*(e < 0),k=asymadj[i]))
    }
  colnames(asymlags) <- paste(c("asym"),asymadj,sep="")
  s2 <- max(s2,I(max(asym)+1))
  }
mX <- cbind(mX, asymlags)

#create log-EWMA term:
if(is.null(log.ewma)){logEWMA <- NULL}else{
  logEWMA <- do.call(leqwma, c(list(e),log.ewma) )
}
mX <- cbind(mX, logEWMA)

#create matrix of variance regressors vx:
if(is.null(vx)){vX <- NULL}else{
 vxnames <- colnames(vx)
 if(is.null(vxnames))
   {
   vxnames <- paste("vx",1:ncol(vx),sep="")
   }else
   {
   for(i in 1:length(vxnames))
     {
     if(vxnames[i]==""){vxnames[i] <- paste("vx",i,sep="")}
     }
   }
 vX <- cbind(vx)
 colnames(vX) <- vxnames
 }
mX <- cbind(mX, vX)

#out-matrix:
out <- mX
return(out)

}
