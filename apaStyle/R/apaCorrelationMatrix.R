##
#' Generic method to generate a correlation matrix with values
#'
#' @param  data Raw dataset with variables.
#' @param  position (optional) Specify whether the correlations should be displayed in the \code{upper}, or \code{lower} diagonal of the table.
#' @return \code{apa.cor.matrix} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{data}{the data with correlation values}
#' @importFrom "stats" "cor" "pf"
#' @export
#'
#' @examples
#'
#' # Use apa.cor.matrix function
#' apa.cor.matrix(
#'   data = data.frame(
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1)
#'   ),
#'   position = "upper"
#' )
##
apa.cor.matrix = function(data=data.frame(), position="lower") UseMethod("apa.cor.matrix")

##
#' Default method to generate a correlation matrix with values
#'
#' @param  data Raw dataset with variables.
#' @param  position (optional) Specify whether the correlations should be displayed in the \code{upper}, or \code{lower} diagonal of the table.
#' @return \code{apa.cor.matrix} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{data}{the data with correlation values}
#' @importFrom "stats" "cor" "pf"
#' @export
#'
#' @examples
#'
#' # Use apa.cor.matrix function
#' apa.cor.matrix(
#'   data = data.frame(
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1)
#'   ),
#'   position = "upper"
#' )
##
apa.cor.matrix.default = function(data=data.frame(), position="lower") {

  est = apaStyleCorrelation(data, position)
  est$call = match.call()
  class(est) = "apa.cor.matrix"
  est

}

##
#' Define a print method
#'
#' @param  x A \code{apa.cor.matrix} object
#' @export
##
print.apa.cor.matrix = function(x, ...) {
  if(x$succes == TRUE) {
    cat("\n")
    cat("Correlation matrix is succesfully generated.")
    cat("\n\n")
  }
}

# The main function
apaStyleCorrelation = function(data, position) {

  # Initialize function
  options(warn = 0)

  # Check if a valid data frame is supplied
  if ((!is.data.frame(data)) || (is.data.frame(data) && nrow(data)==0)) {
    error = "Invalid data is supplied."
    warning(error)
    return(list(succes = error))
  }

  # Check if a valid correlation matrix position is specified
  if((!is.character(position)) || (length(position) > 1) || (!"upper" %in% position && !"lower" %in% position)) {
    error = "The supplied display position for the correlation matrix is not valid. Only 'upper' or 'lower' position is allowed."
    warning(error)
    return(list(succes = error))
  }

  # Internal function to calculate significance
  cor.prob = function(prob.data, position, df = nrow(prob.data) - 2) {
    correlations = stats::cor(prob.data, method = "spearman", use = "complete")
    if (position == "upper") {
      pos = row(correlations) < col(correlations)
    } else {
      pos = row(correlations) > col(correlations)
    }
    r2 = correlations[pos]^2
    Fstat = r2 * df / (1 - r2)
    correlations[pos] = round(1 - stats::pf(Fstat, 1, df), digits = 3)
    return(correlations)
  }

  cor.sig = cor.prob(data, position)
  cor.val = round(stats::cor(data, method = "spearman", use = "complete"), digits = 2)

  if (position == "upper") {
    cor.sig[lower.tri(cor.val, diag=TRUE)] = NA
    cor.val[lower.tri(cor.val, diag=TRUE)] = NA
  } else {
    cor.sig[upper.tri(cor.val, diag=TRUE)] = NA
    cor.val[upper.tri(cor.val, diag=TRUE)] = NA
  }

  cor.sig = ifelse(is.na(cor.sig), "", ifelse(cor.sig < .001, "***", ifelse(cor.sig < .01, "**", ifelse(cor.sig < .05, "*", ifelse(cor.sig < .10, "\u2020", paste(c(rep("\u00A0", 6)), collapse = ""))))))
  cor.val = ifelse(is.na(cor.val), "", sprintf("%3.2f", round(cor.val, digits = 2)))

  # Merge correlation values and significance together
  cor.odd = cor.even = 0
  cor.tmp = matrix(NA, nrow = nrow(cor.val), ncol = (ncol(cor.val)*2))
  if (ncol(cor.val) > 0) {
    for(i in 1:ncol(cor.tmp)) {
      if(i %% 2){
        cor.odd = cor.odd + 1
        cor.tmp[,i] = cor.val[,cor.odd]
      } else {
        cor.even = cor.even + 1
        cor.tmp[,i] = cor.sig[,cor.even]
      }
    }
  }

  # Add a void in empty columns so that the cell width corresponds with that of the cells with content
  if (position == "lower") {
    cor.tmp[nrow(cor.tmp),ncol(cor.tmp)-1] = paste(c(rep("\u00A0", 8)), collapse = "")
    cor.tmp[nrow(cor.tmp),ncol(cor.tmp)] = paste(c(rep("\u00A0", 6)), collapse = "")
  } else {
    cor.tmp[1,1] = paste(c(rep("\u00A0", 8)), collapse = "")
    cor.tmp[1,2] = paste(c(rep("\u00A0", 6)), collapse = "")
  }


  return(list(succes = TRUE, data = cor.tmp))

}
