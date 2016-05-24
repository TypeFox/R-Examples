#
# Post.processing of MCMC simulation
# Example of fitted time series of most active voxel 
#
post.tseries <-
function(out, slicedata, ymaskdata, vreg=2)
{
    X <- slicedata$X
		nvar <- slicedata$nvar 
		nobs <- slicedata$nobs 
    yn <- ymaskdata$yn
    nreg <- ymaskdata$nreg
    stopifnot(nobs == nrow(yn))
    #-----------------------------
    # Postprocessing 
    pmeans <- pmeans.hcoef(out$betadraw) 
    px <- regpostsim(pmeans, vreg=vreg, plot=F)
    pm2 <- pmeans[,vreg]
    spma <- px$spma
    spmn <- px$spmn
    #-----------------------------
     if(length(spma)) {
         pm2rg <- range(pm2)
         cat("range pm2:",pm2rg,"\n")
         pxa <- which(pm2 == pm2rg[2]) # most active
         # pxn <- which(pm2 == pm2rg[1]) # most non-active
         betabar.a <- pmeans[pxa,]
         # betabar.n <- pmeans[pxn,]
         yfa <- X%*%betabar.a
         # yfn <- X%*%betabar.n
         #--------------
         # x11(width=7, height=3.5)
         # par(mfrow=c(2,1), mar=c(4,2,2,1)+0.1)
         # par(mfrow=c(2,1), mar=c(4,2,2,1)+0.1)
         #--------------
         ylim <- range(yn[,pxa])
         plot(yn[,pxa], ty="l", ylab="", main="Example of fitted time-series for active voxel", lty="dotted", ylim=ylim)
         points(yn[,pxa])
         lines(yfa)
         #----------
         # x11(width=7, height=3)
         # par(mfrow=c(2,1), mar=c(4,2,2,1)+0.1)
         # par(mar=c(4,2,2,1)+0.1)
         # plot(yn[,pxn], ty="l", ylab="", main="non-activated voxel", lty="dotted", ylim=ylim)
         # points(yn[,pxn])
         # lines(yfn)
     } else {
        cat("\nNo active voxels detected !\n");
    }
}
