### Fit multinomial logistic regression and find Delta.t and log(mu) (selection and mutation
### effects) based on vglm() in VGAM package.
###
### These functions are for one amino acid.

### Get the specific function according to the options.
get.my.fitMultinomOne <- function(model){
  if(!any(model[1] %in% .CF.CT$model)){
    stop("model is not found.")
  }
  ret <- eval(parse(text = paste("my.fitMultinomOne.", model[1], sep = "")))
  assign("my.fitMultinomOne", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.fitMultinomOne().


### Fit a ROC + NSEf model.
my.fitMultinomOne.rocnsef <- function(reu13.df.aa, phi, yaa, naa,
    phi.new = NULL, coefstart = NULL, x.arg = FALSE, y.arg = FALSE,
    qr.arg = FALSE){
  ### If phi.new is not NULL, it means the caller is in MCMC steps.
  if(!is.null(phi.new)){
    tmp.phi <- rep(phi.new, naa)
  } else{
    tmp.phi <- reu13.df.aa$phi
  }

  ### Obtain new beta (M, S_1, S_2) from vglm.
  ret <- VGAM::vglm(reu13.df.aa$Codon ~ tmp.phi + tmp.phi:reu13.df.aa$Pos,
                    VGAM::multinomial, coefstart = coefstart,
                    x.arg = x.arg, y.arg = y.arg, qr.arg = qr.arg)
  ret <- list(coefficients = ret@coefficients,
              coef.mat = matrix(ret@coefficients, nrow = 3, byrow = TRUE),
              R = ret@R)
  ret
} # End of my.fitMultinomOne.rocnsef().

### Fit a ROC model.
my.fitMultinomOne.roc <- function(reu13.df.aa, phi, yaa, naa, phi.new = NULL,
    coefstart = NULL, x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE){
  ### If phi.new is not NULL, it means the caller is in MCMC steps.
  if(!is.null(phi.new)){
    tmp.phi <- phi.new
  } else{
    tmp.phi <- phi
  }

  ### Extra care for ROC since some genes have total zero counts for some amino
  ### acid. Avoid those cases since there is no information for the parameters.
  names.phi <- names(phi)
  tmp.id <- names.phi %in% names(naa)[naa > 0]

  ### Obtain new beta (M, S_1) from vglm.
  ret <- VGAM::vglm(yaa[tmp.id,] ~ tmp.phi[tmp.id],
                    VGAM::multinomial, coefstart = coefstart,
                    x.arg = x.arg, y.arg = y.arg, qr.arg = qr.arg)
  coefficients <- ret@coefficients
  ## convert delta.t to delta.eta
  #coefficients[(length(ret@coefficients)/2+1):length(coefficients)] <- -1 * coefficients[(length(ret@coefficients)/2+1):length(coefficients)]
  coefficients <- -coefficients
  ret <- list(coefficients = coefficients,
              coef.mat = matrix(ret@coefficients, nrow = 2, byrow = TRUE),
              R = ret@R)
  ret
} # End of my.fitMultinomOne.roc().

### Fit a NSEf model.
my.fitMultinomOne.nsef <- function(reu13.df.aa, phi, yaa, naa, phi.new = NULL,
    coefstart = NULL, x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE){
  ### If phi.new is not NULL, it means the caller is in MCMC steps.
  if(!is.null(phi.new)){
    tmp.phi <- rep(phi.new, naa)
  } else{
    tmp.phi <- reu13.df.aa$phi
  }

  ### Obtain new beta (M, S_2) from vglm.
  ret <- VGAM::vglm(reu13.df.aa$Codon ~ tmp.phi:reu13.df.aa$Pos,
                    VGAM::multinomial, coefstart = coefstart,
                    x.arg = x.arg, y.arg = y.arg, qr.arg = qr.arg)
  ret <- list(coefficients = ret@coefficients,
              coef.mat = matrix(ret@coefficients, nrow = 2, byrow = TRUE),
              R = ret@R)
  ret
} # End of my.fitMultinomOne.nsef().

