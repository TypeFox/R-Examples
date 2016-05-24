logLik.mlgarch <-
function(object, varma=FALSE, ...)
{
  if(varma==TRUE){
    result <- object$objective.varma
    attr(result, "df") <- length(object$par.varma)
  }else{
    mZhat <- residuals.mlgarch(object)
    mSigmaFit <- fitted.mlgarch(object)
    if(object$aux$yanyrowiszeron > 0){
      mZhat <- mZhat[-object$aux$yzerowhichrows,]
      mSigmaFit <- mSigmaFit[-object$aux$yzerowhichrows,]
    }
    mR <- cor(mZhat)
    mRlnDet <- as.numeric(determinant(mR, logarithm=TRUE)$modulus)
    mRinv <- solve(mR)
    tmZhatmRmZhat <- rowSums( (mZhat%*%mRinv)*mZhat )
    term1 <- -object$aux$ynonzerorowsn*(object$aux$m*log(2*pi) + mRlnDet)/2
    term2 <- -sum(rowSums(log(mSigmaFit)))
    term3 <- -sum(tmZhatmRmZhat)/2
    result <- term1 + term2 + term3
    attr(result, "df") <- length(object$par) - object$aux$m
  }
#  attr(result, "nobs") <- length(object$aux$n)
  attr(result, "nobs") <- length(object$aux$ynonzerorowsn)
  class(result) <- "logLik"
  return(result)
}
