#' @export
#' @import stats
poisLikSI <- function(dat, level, tol=.001, conf=F, B=500){
  f <- function(x){exp(sapply(x, function(y) {sum(log(dpois(dat, y)))})-sum(log(dpois(dat, mean(dat)))))-1/level}
  llim <- ifelse(sign(f(0))==sign(f(mean(dat))),0,uniroot(f, c(0, mean(dat)), tol=tol)$root)
  ulim <- uniroot(f, c(mean(dat), max(dat)+10000), tol=tol)$root
  if(conf==F) {return(c(llim, ulim))
  } else{
    ct <- 0
    for(i in 1:B){
      bsam <- rpois(length(dat), mean(dat))
      psi <- poisLikSI(bsam, level)
      if(psi[1]<mean(dat) & psi[2]>mean(dat)) {ct <- ct+1}
    }
    out <- list(c(llim, ulim), ct/B)
    names(out) <- c("si", "conf.equiv")
    return(out)
  }
}