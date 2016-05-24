#' @export
#' @import stats
normLikSI <- function(dat, level, tol=.001, conf=F){
  pl <- function(a) {sapply(a, function(y) mean((dat-y)^2)^(-length(dat)/2))/mean((dat-mean(dat))^2)^(-length(dat)/2)-1/level}
  llim <- uniroot(pl, c(min(dat)-2*sd(dat), mean(dat)), tol=tol)$root
  ulim <- uniroot(pl, c(mean(dat), max(dat)+2*sd(dat)), tol=tol)$root
  if(conf==F) {return(c(llim, ulim))
  } else{
    out <- list(c(llim, ulim), pchisq(-2*log(1/level),1))
    names(out) <- c("si", "conf.equiv")
    return(out)
  }
}
