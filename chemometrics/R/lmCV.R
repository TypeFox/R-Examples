lmCV <-
function(formula,data, repl=100, 
  segments=4,segment.type = c("random", "consecutive", "interleaved"),length.seg,
  trace = FALSE, ...)
{
# Repeated Cross Validation for multiple linear regression
#    see also function "mvr" and "mvrCv" from library(pls)

#require(pls)

    mf <<- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- delintercept(model.matrix(mt, mf))

n <- nrow(X)

    pred <- matrix(NA,nrow=n,ncol=repl) # collect all predicted values here
    for (i in 1:repl){
       if (missing(length.seg)) {
            segment <- cvsegments(n, k = segments, type = segment.type)
        }
        else {
            segment <- cvsegments(n, length.seg = length.seg, type = segment.type)
        }
      if (trace)
        cat(paste("Replication ",i,": ",sep=""))
      for (n.seg in 1:length(segment)) {
        if (trace)
            cat(n.seg, "")
        seg <- segment[[n.seg]]
        obsuse <- as.numeric(unlist(segment[-n.seg]))
        res <- lm(y~X,subset=obsuse)
	pred[seg,i] <- res$coef[1] + X[seg,]%*%res$coef[-1]
      }
}
    resid <- pred-y
    SEP <- apply(resid,2,sd)
    SEPm <- mean(SEP)
    RMSEP <- sqrt(apply(resid^2,2,mean))
    RMSEPm <- mean(RMSEP)
    list(residuals=resid, predicted=pred, SEP=SEP, SEPm=SEPm, RMSEP=RMSEP, RMSEPm=RMSEPm)
}

