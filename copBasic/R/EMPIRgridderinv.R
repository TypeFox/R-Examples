"EMPIRgridderinv" <-
function(empgrid=NULL, kumaraswamy=FALSE, ...) {

  if(is.null(empgrid)) {
    warning("The gridded empirical copula (say from EMPIRgrid) is NULL")
    return(NULL)
  }

  the.deriv <- EMPIRgridder(empgrid=empgrid,...)

  F <- empgrid$v
  n <- length(F)

  # use a mid-point integration method so 1/2s are needed
  # on the tails.
  delF <- diff(F)
  delF <- c(0.5*delF[1], delF)
  delF[n] <- 0.5*delF[n]

  the.inverse <- matrix(nrow=n, ncol=n)
  Alphas <- Betas <- vector(mode="numeric", length=n)
  Alphas[1] <- NA; Betas[1] <- NA
  for(i in 2:n) {
    x <- the.deriv[i,]

    # invert the CDF by linear approximation
    # we want the QDF with F (horizontal axis values on same spacing)
    inv <- approx(x, y=F, xout=F, rule=2)$y


    if(kumaraswamy) {
      beta0 <- sum(inv*delF) # first PWM (mean)
      beta1 <- sum(inv*F*delF) # second PWM (no name)

      lmr <- lmomco::vec2lmom(c(beta0, 2*beta1 - beta0)) # L-moments
      par.of.kur <- lmomco::parkur(lmr)
      X <- lmomco::quakur(F, par.of.kur) # Kumuraswamy distribution

      the.inverse[i,] <- X
      Alphas[i] <- par.of.kur$para[1]
      Betas[i]  <- par.of.kur$para[2]
    } else {
       the.inverse[i,] <- inv
       Alphas[i] <- NA
       Betas[i]  <- NA
    }
  }
  attributes(the.inverse) <- list(dim=dim(empgrid$empcop),
                                  rownames=empgrid$u,
                                  colnames=empgrid$v,
                                  kumaraswamy=list(Alpha=Alphas,
                                                   Beta=Betas),
                                  message="use the rows!, wrt U")

  return(the.inverse)
}
