
#' @export
#' @import stats
binLikSI <- function(dat, n, level, tol=.001, conf=F, B=500){
  f <- function(x){sapply(x, function(y) {dbinom(dat, n, y)})/dbinom(dat, n, dat/n)-1/level}
  llim <- ifelse(sign(f(0))==sign(f(dat/n)),0,uniroot(f, c(0, dat/n), tol=tol)$root)
  ulim <- ifelse(sign(f(dat/n))==sign(f(1)), 1,uniroot(f, c(dat/n, 1), tol=tol)$root)
  if(conf==F) {return(c(llim, ulim))
  } else {
    ct <- 0
    for(i in 1:B){
      p <- rbeta(1, dat+1, (n-dat)+1)
      bsam <- rbinom(1, n, p)
      bsi <- binLikSI(bsam, n, level)
      if(bsi[1]<=p & bsi[2]>=p) {ct <- ct+1}
    }
    out <- list(c(llim, ulim), ct/B)
    names(out) <- c("si", "conf.equiv")
    return(out)
  }
}
