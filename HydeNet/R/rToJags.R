#' @name TranslateFormula
#' @importFrom stringr str_trim
#' 
#' @title Translate R Formula to JAGS
#' @description While most functions available in JAGS have equivalents in R,
#'   they don't always use the exact same names.  R formulas are converted
#'   to character strings, function names translated, and the corresponding
#'   JAGS formula is returned.
#'   
#' @param f R formula object
#' 
#' @details Only a limited subset of R functions are recognized here, but no 
#'   attempt is made to restrict the user to functions that will be recognized
#'   by JAGS.  For now, the user should remain aware of what functions
#'   are available in JAGS and only use the corresponding functions in R.
#'   The JAGS functions may be referenced in the JAGS user manual (see 
#'   References).  The corresponding R functions are listed in the 
#'   \code{jagsFunctions} data set (use \code{data(jagsFunctions)} to 
#'   review).
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' @references \url{http://people.math.aau.dk/~kkb/Undervisning/Bayes14/sorenh/docs/jags_user_manual.pdf}

rToJags <- function(f){
  f <- as.character(f)
  f <- unlist(strsplit(f, "[+]"))
  f <- stringr::str_trim(f)
  
  #* Easy translations
  f <- gsub("acos[(]", "arccos(", f)
  f <- gsub("acosh[(]", "arccosh(", f)
  f <- gsub("asin[(]", "arcsin(", f)
  f <- gsub("asinh[(]", "arcsinh(", f)
  f <- gsub("atan[(]", "arctan(", f) 
  f <- gsub("pnorm[(]", "phi(", f) 
  f <- gsub("ceiling[(]", "round(", f) 
  f <- gsub("floor[(]", "trunc(", f) 

  #* Not so easy translations
  #* convert equals
  #* 21 November 2014: Testing functionality without this.  Strictly
  #*   speaking, this isn't necessary as the == notation works fine.  
  #*   In fact, it may be preferable for dealing with factors.
#   convertEquals <- function(x){
#     if (grepl("[=][=]", x)){
#       x <- stringr::str_trim(unlist(strsplit(x, "[=][=]")))
#       for (i in length(x):2){
#         x[i-1] <- paste0('equals(', x[i-1], ", ", x[i], ")")
#         x <- x[1:(i-1)]
#       }
#       return(x)
#     }
#     else return(x)
#   }
#   f <- sapply(f, convertEquals)
  
  #* convert caret to pow()
  convertCaret <- function(x){
    if (grepl("^", x, fixed=TRUE)){
      x <- stringr::str_trim(unlist(strsplit(x, "^", fixed=TRUE)))
      for (i in length(x):2){
        x[i-1] <- paste0('pow(', x[i-1], ", ", x[i], ")")
        x <- x[1:(i-1)]
      }
      return(x)
    }
    else return(x)
  }
  f <- sapply(f, convertCaret)
  
  #* Convert logit with inverse to ilogit
  convertLogit <- function(x){
    x <- gsub(" ", "", x)
    if (grepl("logit[(]", x)){
      if (grepl("inverse[=]T", x)){
        x <- gsub("(logit[(]|qlogis[(])", "ilogit(", x)
      }
      x <- gsub(",[[:print:]]+[)]", ")", x)
      return(x)
    }
    else return(x)
  }
  f <- sapply(f, convertLogit)
  
  return(paste(f[2], f[1], paste(f[-(1:2)], collapse=" + ")) )
}
