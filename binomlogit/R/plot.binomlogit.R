plot.binomlogit <- function(x, auto.layout = TRUE, ...){     
    oldpar=NULL
    on.exit(par(oldpar))
    if(auto.layout){
        mfcol=c(2,min(4,x$dims))
        mai=c(0.5,0.5,0.2,0.2)             # c(bottom, left, top, right)
        oldpar=par(mfcol=mfcol,mai=mai)
    }
    for(i in 1:x$dims){
      plot((x$burn+1):x$sim,x$beta[i,][(x$burn+1):x$sim],pch=20,col="red",type="l",main=paste("Draws beta",i-1,sep=" "),xlab="",ylab="")
      acf(x$beta[i,(x$burn+1):x$sim],xlab="",ylab="")
    }
}
