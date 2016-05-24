ridgeCV <-
function(formula,data,lambdaopt, repl=5, 
  segments=10,segment.type = c("random", "consecutive", "interleaved"),length.seg,
  trace = FALSE, plot.opt=TRUE, ...)
{
# Repeated Cross Validation  for ridge regression
#

#require(pls)
#require(MASS)

    mf <<- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- delintercept(model.matrix(mt, mf))

ym <- mean(y)
Xsc <- apply(X,2,sd)
Xm <- apply(X,2,mean)
obsuse <- NULL

n <- nrow(X)

    pred <- matrix(NA,nrow=n,ncol=repl) # collect all predicted values here
    for (i in 1:repl){

	print(i)

       if (missing(length.seg)) {
            segment <- cvsegments(n, k = segments, type = segment.type)
        }
        else {
            segment <- cvsegments(n, length.seg = length.seg, type = segment.type)
        }
      if (trace)
        cat(paste("Replication ",i,": ",sep=""))
      for (n.seg in 1:length(segment)) {
        # if (trace) { cat(n.seg, "")}
        seg <- segment[[n.seg]]
        obsuse <- as.numeric(unlist(segment[-n.seg]))
	Xsub <- X[obsuse,]
	ysub <- y[obsuse]
	datsub <- list(ysub=ysub,Xsub=as.matrix(Xsub))	
        mod_ridge <- lm.ridge(ysub~Xsub,data=datsub,lambda=lambdaopt)
        # ypred <- ym-sum(Xm*mod_ridge$coef/Xsc)+X[seg,]%*%(mod_ridge$coef/Xsc)
        ypred <- mean(ysub)-sum(mod_ridge$xm*mod_ridge$coef/mod_ridge$scale)+
		 #X[seg,]%*%(mod_ridge$coef/mod_ridge$scale)
		 X[-obsuse,]%*%(mod_ridge$coef/mod_ridge$scale)
	#pred[seg,i] <- ypred
	pred[-obsuse,i] <- ypred
      }
    }
    resid <- pred-y
    SEP <- apply(resid,2,sd)
    SEPm <- mean(SEP)
    sMAD <- apply(resid,2,mad)
    sMADm <- mean(sMAD)
    RMSEP <- sqrt(apply(resid^2,2,mean))
    RMSEPm <- mean(RMSEP)

if (plot.opt){
par(mfrow=c(1,2))
ylimits=range(pred)
plot(y,pred[,1],xlab="Measured y",ylab="Predicted y",cex.lab=1.2,cex=0.7,
#   type="n",pch=3,col=1,ylim=ylimits,...)
   type="n",pch=3,col=1,ylim=ylimits)
for (i in 1:ncol(pred)){
  points(y,pred[,i],pch=3,cex=0.7,col=gray(0.6))
}
predm=apply(pred,1,mean)
points(y,predm,pch=3,cex=0.7,col=1)
title("Predictions from Repeated CV")
abline(c(0,1),lty=1)

ylimits=range(predm)
plot(y,predm,xlab="Measured y",ylab="Predicted y",cex.lab=1.2,cex=0.7,
#   pch=3,col=1,ylim=ylimits,...)
   pch=3,col=1,ylim=ylimits)
title("Average of predictions")
abline(c(0,1),lty=1)

legend("topleft",c(paste("SEP=",round(SEPm,2),sep=""),
  paste("sMAD=",round(sMADm,2),sep="")))
}


list(residuals=resid, predicted=pred, SEP=SEP, SEPm=SEPm, sMAD=sMAD,sMADm=sMADm,
         RMSEP=RMSEP, RMSEPm=RMSEPm)
}

