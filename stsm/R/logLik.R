
logLik.stsm <- function(object, domain = c("frequency", "time"), xreg = NULL,
  td.args = list(P0cov = FALSE, t0 = 1, 
    KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version)), 
  check.td.args = TRUE, 
  barrier = list(type = c("1", "2"), mu = 0), 
  inf = 99999, ...)
{
  if(length(list(...)))
    warning("Extra arguments '...' were ignored.")

  domain <- match.arg(domain)[1]

  mll <- switch(domain,
    "frequency" = -mloglik.fd(NULL, object, barrier, inf),

    "time" = {
      if (is.null(td.args$KF.version)) {
        KF.version <- "KFKSDS"
      } else {
        KF.version <- td.args$KF.version
        td.args$KF.version <- NULL
      }
      -mloglik.td(NULL, object, KF.version, td.args, check.td.args, barrier, inf)
    })

  #attr(mll, "nobs") <- length(object@y)
  attr(mll, "df") <- length(object@pars)
  class(mll) <- "logLik"
  mll
}
