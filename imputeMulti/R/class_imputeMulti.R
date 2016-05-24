
#' Class "mod_imputeMulti"
#'  
#' @name mod_imputeMulti-class
#' @description A multivariate multinomial model imputed by EM or Data Augmentation is 
#' represented as a \code{\linkS4class{mod_imputeMulti}} object. A complete 
#' dataset and model is represented as an \code{\linkS4class{imputeMulti}} object.
#' Slots for \code{mod_imputeMulti} objects include: (1) the modeling method; 
#' (2) the call to the estimation function; (3) the number of iterations in estimation;
#' (4) the final log-likelihood; (5) the conjugate prior if any; (6) the MLE estimate of
#' the sufficient statistics and parameters.
#' @slot method the modeling method
#' @slot mle_call the call to the estimation function
#' @slot mle_iter the number of iterations in estimation
#' @slot mle_log_lik the final log-likelihood
#' @slot mle_cp the conjugate prior if any
#' @slot mle_x_y the MLE estimate of the sufficient statistics and parameters
#' @docType class
#' @section Objects from the class: Objects are created by calls to
#' \code{\link{multinomial_impute}}, \code{\link{multinomial_em}}, or
#' \code{\link{multinomial_data_aug}}.
#' @section Methods:
# \code{\link{show,mod_imputeMulti-method}}
#' \code{\link{summary,mod_imputeMulti-method}}
#' @seealso \code{\link{multinomial_impute}}, \code{\link{multinomial_em}}, 
#' \code{\link{multinomial_data_aug}}
#' @exportClass mod_imputeMulti
#' @export
setClass("mod_imputeMulti",
         representation= list(
           method= "character",
           mle_call= "call",
           mle_iter= "numeric",
           mle_log_lik= "numeric",
           mle_cp= "character",
           mle_x_y= "data.frame"),
         validity= function(object) {
           if (!object@method %in% c("EM", "DA", "NULL")) {
             return("Currently only EM and DA methods are defined.")
           } else if (object@mle_iter < 0) {
             return("A negative iteration was given.")
           }
           return(TRUE)
         }
)
  

###########################################################
## Methods
###########################################################
## Print

show.mod_imputeMulti <- function(object) {
  cat("\n Call: \n", paste(deparse(object@mle_call)),
      "\n Method: ", object@method,
      "\n\n Iterations: ", object@mle_iter,
      "\n\n Log-Likelihood: ", object@mle_log_lik)
}

# @title Print mod_imputeMulti class objects
# @description print method for class "mod_imputeMulti"
# @param object an object of class "mod_imputeMulti"
# @param ... further arguments passed to or from other methods.
# @exportMethod show
# @export
# setMethod("show", signature= "mod_imputeMulti",
#           def= show.mod_imputeMulti)


## Summary
setGeneric("summary")
# @export
summary.mod_imputeMulti <- function(object, ...) {
  print(object)
  
  summary.data.frame(object@mle_x_y[, c("alpha", "theta_y")])
}

#' @title Summarizing mod_imputMulti objects
#' @description summary method for class "mod_imputeMulti"
#' @param object an object of class "mod_imputeMulti"
#' @param ... further arguments passed to or from other methods.
#' @exportMethod summary
setMethod("summary", signature="mod_imputeMulti", def=summary.mod_imputeMulti)


#' Class "imputeMulti" 
#'  
#' @name imputeMulti-class
#' @description A multivariate multinomial model imputed by EM or Data Augmentation is 
#' represented as a \code{\linkS4class{mod_imputeMulti}} object. A complete 
#' dataset and model is represented as an \code{\linkS4class{imputeMulti}} object.
#' Inherits from \code{mod_imputeMulti}. Additional slots are supplied for (1) the
#' call to \code{multinomial_impute}; (2) the missing and imputed data;
#' and (3) the number of observations with missing values.
#' @slot Gcall the call to \code{multinomial_impute}
#' @slot method the modeling method
#' @slot mle_call the call to the estimation function
#' @slot mle_iter the number of iterations in estimation
#' @slot mle_log_lik the final log-likelihood
#' @slot mle_cp the conjugate prior if any
#' @slot mle_x_y the MLE estimate of the sufficient statistics and parameters
#' @slot data a \code{list} of the missing and imputed data
#' @slot nmiss the number of observations with missing data
#' @docType class
#' @section Objects from the class: Objects are created by calls to
#' \code{\link{multinomial_impute}}, \code{\link{multinomial_em}}, or
#' \code{\link{multinomial_data_aug}}.
#' @section Methods:
# \code{\link{show,imputeMulti-method}}
#' \code{\link{summary,imputeMulti-method}}
#' @seealso \code{\link{multinomial_impute}}, \code{\link{multinomial_em}}, 
#' \code{\link{multinomial_data_aug}}
#' @exportClass imputeMulti
#' @export
setClass("imputeMulti",
         representation= list(Gcall= "call",
                              method= "character",
                              mle_call= "call",
                              mle_iter= "numeric",
                              mle_log_lik= "numeric",
                              mle_cp= "character",
                              mle_x_y= "data.frame",
                              data= "list",
                              nmiss= "numeric"),
         contains= "mod_imputeMulti")



###########################################################
## Methods
###########################################################
## Print

show.imputeMulti <- function(object) {
  cat("\n Global Call: \n", paste(deparse(object@Gcall)),
      "\n Call: \n", paste(deparse(object@mle_call)),
      "\n Method: ", object@method,
      "\n\n Iterations: ", object@mle_iter,
      "\n\n Log-Likelihood: ", object@mle_log_lik,
      "\n Number Missing: ", object@nmiss)
}

# @title Print imputeMulti class objects
# @description print method for class "imputeMulti"
# @param object an object of class "imputeMulti"
# @param ... further arguments passed to or from other methods.
# @exportMethod show
# @export
# setMethod("show", signature= "imputeMulti",
#   def= show.imputeMulti)

## Summary
setGeneric("summary")
# @export
summary.imputeMulti <- function(object, ...) {
  print(object)
  
  summary(object@mle_x_y[, c("alpha", "theta_y")])
  
  lapply(object@data, dim)
}

#' @title Summarizing imputMulti objects
#' @description summary method for class "imputeMulti"
#' @param object an object of class "imputeMulti"
#' @param ... further arguments passed to or from other methods.
#' @exportMethod summary
setMethod("summary", signature= "imputeMulti",
  summary.imputeMulti)


###########################################################
# inheritance checks

#' @title Check mod_imputeMulti Class
#' @description Function that checks if the target object is a \code{mod_imputeMulti} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "mod_imputeMulti" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.mod_imputeMulti <- function(x) {
  inherits(x, "mod_imputeMulti")
}


#' @title Check imputeMulti Class
#' @description Function that checks if the target object is a \code{imputeMulti} object.
#' @param x any R object.
#' @return Returns \code{TRUE} if its argument has class "imputeMulti" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.imputeMulti <- function(x) {
  inherits(x, "imputeMulti")
}
