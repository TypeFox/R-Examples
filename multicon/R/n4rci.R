n4rci <-
function(CIwidth=NULL, N=NULL, alpha=NULL) {
  if(is.null(N) & !is.null(CIwidth) & !is.null(alpha)) {
    a <- qnorm(1 - alpha/2)
    N <- (a / CIwidth)^2 + 3
    out <- cbind(N, CIwidth, alpha)
  }
  if(!is.null(N) & is.null(CIwidth) & !is.null(alpha)) {
    if(N < 4) {stop("N must be greater than 3 in this mode.")}
    a <- qnorm(1 - alpha/2)
    CIwidth <- a / sqrt(N-3)
    out <- cbind(N, CIwidth, alpha)
  }
  if(!is.null(N) & !is.null(CIwidth) & is.null(alpha)) {
    if(N < 4) {stop("N must be greater than 3 in this mode.")}
    alpha <- (1-pnorm(sqrt(N - 3) * CIwidth))*2
    out <- cbind(N, CIwidth, alpha)
  }
  out <- cbind(N, CIwidth, alpha)
  colnames(out) <- c("N", "CI Width", "alpha")
  rownames(out) <- "Results"
  return(out)
}
