#' Iterative offset model fitting of mark-recapture with removal model
#' 
#' Detection function fitting from mark-recapture data with a removal
#' configuration in which a secondary observer knows what the primary observer
#' detects and detects objects missed by the primary observer.  The iterative
#' offset glm/gam uses an offset to compensate for the conditioning on the set
#' of objects seen by either observer (eg 00 those missed by both observers are
#' not included in the analysis.  This function is similar to
#' \code{\link{io.glm}}.
#' 
#' The only difference between this function and \code{\link{io.glm}} is the
#' offset and the data construction because there is only one detection
#' function being estimated for the primary observer. The two functions could
#' be merged.
#' 
#' @param datavec dataframe containing records seen by either observer 1 or 2
#' @param fitformula logit link formula
#' @param eps convergence criterion
#' @param iterlimit maximum number of iterations allowed
#' @param GAM uses GAM instead of GLM for fitting
#' @param gamplot set to TRUE to get a gam plot object if \code{GAM=TRUE}
#' @param datavec2 dataframe containing all records for observer 1 and observer
#'   2 as in io.glm form; this is used in case there is an observer(not
#'   platform effect)
#' @return list of class("remglm","glm","lm") or class("remglm","gam")
#'   \item{glmobj}{GLM or GAM object} \item{offsetvalue}{offsetvalues from
#'   iterative fit} \item{plotobj}{gam plot object (if GAM & gamplot==TRUE,
#'   else NULL)}
#' @note currently the code in this function for GAMs has been commented out
#'   until the remainder of the mrds package will work with GAMs.
#' @author Jeff Laake
#' @references Buckland, S.T., J.M. breiwick, K.L. Cattanach, and J.L. Laake.
#'   1993. Estimated population size of the California gray whale.  Marine
#'   Mammal Science, 9:235-249.
#' 
#' Burnham, K.P., S.T. Buckland, J.L. Laake, D.L. Borchers, T.A. Marques,
#'   J.R.B. Bishop, and L. Thomas. 2004.  Further topics in distance sampling.
#'   pp: 360-363. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'   D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'   Oxford University Press.
#' @keywords Statistical Models
rem.glm <-
function(datavec, fitformula, eps = 0.00001, iterlimit = 500, GAM = FALSE, 
  gamplot = TRUE, datavec2)
{
# ---------------------------------------------------------------
#  This is the code that uses the iterative offset glm or gam
#  approach; iteration is done until parameters are within a 
#  certain epsilon (eps) or iteration limit (iterlimit) exceeded.
#
#  Note: David used offset in formula and I've put it as an
#  argument to glm and gam functions.
#
# Input : datavec = dataframe
#         fitformula = formula
#         eps = convergence criterion - fixed
#         iterlimit = maximum number of iterations allowed - fixed
#     gamplot = TRUE if want a gam plot object
#
# Output: list with
#       Note: modified to return glm object only with 
#               class("remglm","glm","lm")
#               class("remglm","gam")
#
# $glmobj:  glm model
# $offsetvalue: final offsetvalues from iterative fit
# $plotobj: gam plot object (if GAM & gamplot==TRUE, else NULL)
# ---------------------------------------------------------------- 
# 
  datavec2$offsetvalue=0
  done <- FALSE
  i <- 1
  plotobj <- NULL
  while(i <= iterlimit & !done) {
#  fit the glm or gam
    if(GAM) {
      ioglm <- mgcv::gam(formula = fitformula, family = binomial, data = datavec)
    }else{
      ioglm <- glm(formula = fitformula, family = binomial, data = datavec)
    }

    coeff <- ioglm$coeff
    fittedp <- ioglm$fitted.values

    if(i == 1) {
      oldmodel <- ioglm
      oldcoeff <- coeff
      oldp <- fittedp
    }else{
#    calculate differences between previous and present set of model outputs
      reldiff <- max(abs(plogis(coeff) - plogis(oldcoeff))/plogis(oldcoeff))
      
      if(is.na(reldiff)) {
        print("Can't calculate regression coefficients - model has not converged")
        print(" - last fit used for estimation" )
        ioglm <- oldmodel
        done <- TRUE
      }
      
      if(reldiff < eps & !done) {
        done <- TRUE
      }else{
        oldmodel <- ioglm
        oldcoeff <- coeff
        oldp <- fittedp
      }
    }
    if(!done){
      oldoff <- datavec$offsetvalue
#     if(GAM)
        off <-  - log(plogis(predict(ioglm,newdata=datavec2)))
#       off <-  - log(plogis(ioglm$
#         additive.predictors - datavec$
#         offsetvalue))
#     else off <-  - log(plogis(ioglm$
#         linear.predictors - datavec$
#         offsetvalue))
      datavec$offsetvalue <- off
    }
    i <- i + 1
  }
# if(GAM & gamplot) {
#   assign("fitformula", fitformula, frame = 1)
#   plotobj <- plot.gam(ioglm, se = T)
#   assign("fitformula", NULL, frame = 1)
# }

  if(!done){
    datavec$offsetvalue <- oldoff
    warning("Iteration limit exceeded - last fit used for estimation")
  }

# list(glm = ioglm, offsetvalue = datavec$offsetvalue, plotobj = 
#   plotobj)
  if(GAM) ioglm$offset=off

  class(ioglm)=c("remglm",class(ioglm))

  return(ioglm)
}
