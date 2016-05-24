#' Extract the Design of a linear mixed model
#' 
#' These functions extract various elements of the design of a fitted
#' \code{lme}-, \code{mer} or \code{lmerMod}-Object.  They are called by
#' \code{exactRLRT} and \code{exactLRT}.
#' 
#' 
#' @aliases extract.lmerModDesign extract.lmeDesign
#' @param m a fitted \code{lme}- or \code{merMod}-Object
#' @return a a list with components
#' \itemize{
#' \item \code{Vr} estimated covariance of the random effects divided by the
#' estimated variance of the residuals
#' \item \code{X} design of the fixed effects
#' \item \code{Z} design of the random effects
#' \item \code{sigmasq} variance of the residuals
#' \item \code{lambda} ratios of the variances of the random effects and the
#' variance of the residuals
#' \item \code{y} response variable
#' }
#' @author Fabian Scheipl, \code{extract.lmerModDesign} by Ben Bolker.
#' Many thanks to Andrzej Galecki and Tomasz Burzykowski for bug fixes.
#' @keywords utilities
#' @examples
#' 
#' library(nlme)
#' design <- extract.lmeDesign(lme(distance ~ age + Sex, data = Orthodont, 
#'                              random = ~ 1))
#' str(design)
#' 
#' @export extract.lmeDesign
#' @importFrom stats complete.cases formula model.frame model.matrix 
#' @importFrom nlme getGroups 
#' @importFrom mgcv tensor.prod.model.matrix
extract.lmeDesign <- function(m)
{
  start.level = 1
  
  data <- if(any(!complete.cases(m$data))){
    warning("Removing incomplete cases from supplied data.") 
    m$data[complete.cases(m$data),]
  } else m$data
  grps <- getGroups(m)
  n <- length(grps)
  X <- list()
  grp.dims <- m$dims$ncol
  Zt <- model.matrix(m$modelStruct$reStruct, data)
  cov <- as.matrix(m$modelStruct$reStruct)
  i.col <- 1
  n.levels <- length(m$groups)
  Z <- matrix(0, n, 0)
  if (start.level <= n.levels) {
    for (i in 1:(n.levels - start.level + 1)) {
      if(length(levels(m$groups[[n.levels-i+1]]))!=1)
      {
        X[[1]] <- model.matrix(~m$groups[[n.levels - i +
            1]] - 1, 
          contrasts.arg = c("contr.treatment",
            "contr.treatment"))
      }
      else X[[1]]<-matrix(1, n, 1)
      X[[2]] <- as.matrix(Zt[, i.col:(i.col + grp.dims[i] -
          1)])
      i.col <- i.col + grp.dims[i]
      Z <- cbind(tensor.prod.model.matrix(X), Z)
    }
    Vr <- matrix(0, ncol(Z), ncol(Z))
    start <- 1
    for (i in 1:(n.levels - start.level + 1)) {
      k <- n.levels - i + 1
      for (j in 1:m$dims$ngrps[i]) {
        stop <- start + ncol(cov[[k]]) - 1
        Vr[ncol(Z)+1-(stop:start),ncol(Z)+1-(stop:start)] <- cov[[k]]
        start <- stop + 1
      }
    }
  }
  X <- if(class(m$call$fixed) == "name" &&  !is.null(m$data$X)){
    m$data$X
  } else 	{
    model.matrix(formula(eval(m$call$fixed)),data)
  }
  y<-as.vector(matrix(m$residuals, ncol=NCOL(m$residuals))[,NCOL(m$residuals)] + 
      matrix(m$fitted, ncol=NCOL(m$fitted))[,NCOL(m$fitted)])
  return(list(
    Vr=Vr, #Cov(RanEf)/Var(Error)
    X=X,
    Z=Z,
    sigmasq=m$sigma^2,
    lambda=unique(diag(Vr)),
    y=y,
    k=n.levels
  )
  )
}

