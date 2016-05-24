#' @rdname EStep
#' @return \code{MStep} returns a list of parameters formatted as described in
#'   \code{\link{rtheta}}.
MStep <- function (x, kappa, meta.special.case = FALSE) {
  mus    <- lapply(1:ncol(kappa),
                   function(j) {colSums(x*kappa[,j])/sum(kappa[,j])})                 
  sigmas <- lapply(1:ncol(kappa),
                   function(j) {cov.wt(x, kappa[,j], method = "ML")$cov})
  names(mus)  <- names(sigmas) <- paste("comp", 1:ncol(kappa), sep = "")
  pie <- colMeans(kappa)
  pie <- pie/sum(pie)
  
  if (meta.special.case) {
    m      <- ncol(kappa)
    d      <- ncol(x)
    pie1    <- colMeans(kappa)[1]; names(pie1) <- NULL
    mu2    <- mean(mus[[2]])
    sigma2 <- sqrt(mean(diag(sigmas[[2]])))    
    rho2   <- (sum(sigmas[[2]])-d*sigma2^2)/(d*(d-1)*sigma2^2)
    par    <- c(pie1 = pie1, mu = mu2, sigma = sigma2, rho = rho2)
    return(meta2full(par, d = d))
  }
  return(list(m = ncol(kappa), d = ncol(x), pie = pie,
              mu = mus, sigma = sigmas))
}