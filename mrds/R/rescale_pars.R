#' Calculate the parameter rescaling for parameters associated with covariates
#'
#' This will calculate the rescaling needed when covariates to be included in the scale of the detection function are "too big". Based on code from \code{\link{optimx}}.
#'
#' Derivative-free methods like nlminb are sensitive to the parameters being poorly scaled. This can also cause problems for quasi-Newton methods too (at least, bad scaling won't _help_ the optimisation). So here we rescale the parameters if necessary (unless we already got scaling from control)
#'
#' @author David L Miller
#' @param initialvalues starting values for the optimisation
#' @param ddfobj detection function object
#' @importFrom stats sd terms
rescale_pars <- function(initialvalues, ddfobj){

  par_scaling <- rep(1, length(initialvalues))

  # from optimx:::optimx.setup, scaletol = 3 (so setting 3 here for consistency)
  # here we use a local copy in scalecheck.R
  if(scalecheck(initialvalues, NA, NA, dowarn=FALSE)$lpratio>3){
    # do the rescaling to the scale parameters only
    # this is (still) a bit hackish
    # divide by the standard deviation of the distances
    #  -- this seems appropriate due to the form of the detection function
    #  may not be appropriate in all settings?

    # match parameter vector indices to indices of model matrix
    #ind <- match(colnames(ddfobj$scale$dm),
    #             names(initialvalues))
    ind <- getpar(ddfobj, index=TRUE)
    ind <- ind[1]:(ind[2]-1)

    # find the scalings
    par_scaling[ind] <- apply(ddfobj$scale$dm, 2, sd)/sd(ddfobj$xmat$distance)

    # ensure that the intercept has scaling 1 & any zero is set back to 1
    par_scaling[abs(par_scaling)<sqrt(.Machine$double.eps)] <- 1

    # set the factor scalings to be 1
    non_factors <- !(colnames(ddfobj$scale$dm) %in%
                     rownames(attr(terms(as.formula(ddfobj$scale$formula)),
                                   "factors")))
    par_scaling[non_factors] <- 1
  }

  # return the parameter rescaling vector
  return(par_scaling)
}
