#' S3 class solarPolygenic.
#'
#' @name solarPolygenicClass
#' @rdname solarPolygenicClass
#'
#' @param x 
#'    An object of class \code{solarPolygenic}.
#' @param object
#'    An object of class \code{solarPolygenic}.
#' @param trait
#'    Logical argument for \code{residuals} method, 
#'    indicating whether values of trait to be returned instead of residuals.
#'    The default value is FALSE.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass solarPolygenic

#--------------------
# Print method
#--------------------

#' @rdname solarPolygenicClass
#' @export
print.solarPolygenic <- function(x, ...)
{
  cat("\nCall: ")
  print(x$call)
  
  cat("\nFile polygenic.out:\n")
  l_ply(x$solar$files$model$polygenic.out, function(x) cat(x, "\n"))
}

#' @rdname solarPolygenicClass
#' @export
summary.solarPolygenic <- function(object, ...)
{
  cat("\nCall: ")
  print(object$call)
  
  cat("\nFile polygenic.out:\n")
  l_ply(object$solar$files$model$polygenic.out, function(x) cat(x, "\n"))
  
  cat("\n Loglikelihood Table:\n")
  print(object$lf)
  
  cat("\n Covariates Table:\n")
  print(object$cf)

  cat("\n Variance Components Table:\n")
  print(object$vcf)
}

#--------------------
# Generic method
#--------------------

#' @rdname solarPolygenicClass
#' @export
residuals.solarPolygenic <- function(object, trait = FALSE, ...)
{
  stopifnot(!is.null(object$resf))
  stopifnot(nrow(object$resf) > 0)
  stopifnot(all(c("id", "residual") %in% names(object$resf)))
  
  if(!trait) {
    r <- object$resf$residual
    names(r) <- object$resf$id
  } else {
    stopifnot(length(object$traits) == 1)
    trait <- object$traits

    trait <- tolower(trait) # SOLAR naming in residual files
    stopifnot(trait %in% names(object$resf))
  
    r <- subset(object$resf, select = trait, drop = TRUE)
    names(r) <- object$resf$id
  }
  
  return(r)
}


#--------------------
# Other method
#--------------------

#' Get formula string from a solarPolygenic object
#' 
#' The function returns a character string with formula.
#' The formula is derived based on \code{traits} and \code{covlist} slots of input object.
#'
#' @param x 
#'    An object of \code{solarPolygenic} object.
#' @return
#'    A character string with formula.
#' 
#' @export
getFormulaStr <- function(x)
{
  paste(
    paste(x$traits, collapse = "+"),
    "~",
    paste(x$covlist, collapse = "+"))
}
