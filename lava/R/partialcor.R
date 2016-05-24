##' Calculate partial correlations
##'
##' Calculate partial correlation coefficients and confidence limits via Fishers
##' z-transform
##'
##'
##' @param formula formula speciying the covariates and optionally the outcomes
##' to calculate partial correlation for
##' @param data data.frame
##' @param level Level of confidence limits
##' @return A coefficient matrix
##' @author Klaus K. Holst
##' @keywords models regression
##' @examples
##'
##' m <- lvm(c(y1,y2,y3)~x1+x2)
##' covariance(m) <- c(y1,y2,y3)~y1+y2+y3
##' d <- sim(m,500)
##' partialcor(~x1+x2,d)
##'
##' @export
partialcor <- function(formula,data,level=0.95) {
  if (attributes(terms(formula))$response==0) {
    preds <- all.vars(formula)
    yy <- setdiff(names(data),preds)
    if (length(yy)<2)
      return(NULL)
    res <- c()
    for (i in seq_len(length(yy)-1))
      for (j in seq(i+1,length(yy))) {
        f <- as.formula(paste("cbind(",yy[i],",",yy[j],")", paste(as.character(formula),collapse="")))
        res <- rbind(res, partialcor(f,data,level=level))
        rownames(res)[nrow(res)] <- paste(yy[i],yy[j],sep="~")
      }
    return(res)
  }
  l <- lm(formula,data)
  k <- ncol(model.matrix(l))
  n <- nrow(model.matrix(l))
  r <- residuals(l)
  rho <- cor(r)[1,2]
  zrho <- atanh(rho)
  var.z <- 1/(n-k-3)
  ci.z <- zrho + c(-1,1)*qnorm(1-(1-level)/2)*sqrt(var.z)
  ci.rho <- tanh(ci.z)
  z <- 1/sqrt(var.z)*zrho
  p.z <- 2*(pnorm(-abs(z))) # p-value using z-transform for H_0: rho=0.
  return(c(cor=rho,z=z,pval=p.z,lowerCI=ci.rho[1],upperCI=ci.rho[2]))
}
