FLLat.PVE <- function(Y,J.seq=seq(1,min(15,floor(ncol(Y)/2)),by=2),
                      B=c("pc","rand"),lams=c("same","diff"),thresh=10^(-4),
                      maxiter=100,maxiter.B=1,maxiter.T=1) {

  ## Error checking.
  if (!is.vector(J.seq)) {
    stop("'J.seq' must be a vector of integers")
  }
  B <- match.arg(B)
  lams <- match.arg(lams)
  CheckPars(Y=Y,B=B,thresh=thresh,maxiter=maxiter,maxiter.B=maxiter.B,
            maxiter.T=maxiter.T)

  ## Total sum of squares.
  tss <- sum((scale(Y,scale=F))^2)

  ## PVEs.
  pves <- rep(0,length(J.seq))

  ## Obtaining lam1 and lam2.
  if (lams=="same") {
    bic <- FLLat.BIC(Y,B=B,thresh=thresh,maxiter=maxiter,
                     maxiter.B=maxiter.B,maxiter.T=maxiter.T)
  }

  for (i in 1:length(J.seq)) {
    if (lams=="diff") {
      est <- FLLat.BIC(Y,J.seq[i],B,thresh,maxiter,maxiter.B,
                       maxiter.T)$opt.FLLat
    } else {
      est <- FLLat(Y,J.seq[i],B,bic$lam1,bic$lam2,thresh,maxiter,
                   maxiter.B,maxiter.T)
    }
    pves[i] <- 1-est$rss/tss
  }

  results <- list("PVEs"=pves,"J.seq"=J.seq)
  class(results) <- "PVE"

  return(results)

}
