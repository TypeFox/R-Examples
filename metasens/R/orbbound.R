orbbound <- function(x, k.suspect=1, tau=x$tau, left=NULL,
                     backtransf=x$backtransf){
  
  ## Copas J, Jackson D. A bound for publication bias based on the
  ## fraction of unpublished studies. Biometrics 2004,
  ## Mar;60(1):146-53.
  
  
  meta:::chkclass(x, "meta")
  
  
  if (!(is.numeric(k.suspect)))
    stop("Argument 'k.suspect' must be a numeric vector")
  ##
  if (any(k.suspect < 0))
    stop("Negative values not allowed for argument 'k.suspect'")
  ##
  if (!is.numeric(tau) || length(tau)!=1 || tau < 0)
    stop("Argument 'tau' must be a positive numeric of length 1")
  
  
  if (is.null(left))
    left <- as.logical(sign(metabias(x, meth="linreg", k.min=3)$estimate[1])==1)
  else if (!(length(left)==1 & is.logical(left)))
    stop("Argument 'left' must be a logical of length 1.")
  
  
  sel <- !is.na(x$seTE)

  
  if (x$hakn)
    warning("Hartung-Knapp adjustment not considered to evaluate outcome reporting bias.")
  
  
  if (min(k.suspect)!=0)
    k.suspect <- c(0, k.suspect)
  ##
  k.suspect <- sort(k.suspect)
  
  
  if (left)
    direction.bias <- -1
  else
    direction.bias <-  1
  
  
  maxbias <- ((x$k + k.suspect)/x$k *
              dnorm(qnorm(x$k/(x$k+k.suspect))) *
              sum(sqrt(x$seTE[sel]^2 + tau^2)^(-1)) /
              sum(sqrt(x$seTE[sel]^2 + tau^2)^(-2))
              )
  ##
  maxbias <- direction.bias*maxbias
  
  
  ci.f <- ci(x$TE.fixed  + maxbias, x$seTE.fixed)
  ci.r <- ci(x$TE.random + maxbias, x$seTE.random)
  
  
  res <- list(maxbias=maxbias,
              k.suspect=k.suspect,
              tau=tau,
              fixed=ci.f,
              random=ci.r,
              left=left,
              x=x,
              call=match.call())
  
  res$backtransf <- backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  class(res) <- "orbbound"
  
  res
}
