rpht <-
function(X, f, theta, cmat, conf.level, alternative, R, args)
{
  Btheta <- function(X, i, f)
    {
      XNEW <- as.data.frame(X[i, ])
      est <- theta(X = XNEW, f = f)
      return(est$estimate)
    }
  bargs <- args
  bargs$data <- as.data.frame(X)
  bargs$statistic = Btheta
  bargs$strata = f
  bargs$f <- f
  bargs$R <- R
  if (is.null(bargs$sim))
    {
      bargs$sim <- "ordinary"
    }
  if (is.null(bargs$stype))
    {
      bargs$stype <- "i"
    }
  bootout <- do.call("boot", bargs)
  diffchains <- Boutrp(x = bootout, cmat = cmat)
  out <- SCIrp(x = diffchains$chains, conf.level = conf.level, alternative = alternative)
  return(out)
}

