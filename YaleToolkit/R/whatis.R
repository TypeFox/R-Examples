#' @title Data frame summary
#' @description Summarize the characteristics of variables (columns)
#' in a data frame.
#' @details The function \code{whatis()} provides a basic examination of some
#' characteristics of each variable (column) in a data frame.
#' @param x a data frame
#' @param var.name.truncate maximum length (in characters) for truncation
#' of variable names.  The default is 20; anything less than 12 is
#' less than the column label in the resulting data frame and is a
#' waste of information.
#' @param type.truncate maximum length (in characters) for truncation of
#' variable type; \code{14} is the full width, but \code{4} works well
#' if space is at a premium.
#' @return A list of characteristics describing the variables in the data
#' frame, \code{x}. Each component of the list has \code{length(x)}
#' values, one for each variable in the data frame \code{x}. 
#' \describe{
#'   \item{variable.name}{from the \code{names(x)} attribute, possibly
#'   truncated to \code{var.name.truncate} characters in length.}
#'   \item{type}{the possibilities include \code{"pure factor"},
#'   \code{"mixed factor"}, \code{"ordered factor"}, \code{"character"}, 
#'   and \code{"numeric"}; \code{whatis()} considers the possibility that a 
#'   factor or a vector could contain character and/or numeric values.  
#'   If both character and numeric values are present, and if the variable 
#'   is a factor, then it is called a mixed factor.  If the levels of a 
#'   factor are purely character or numeric (but not both), it is a pure 
#'   factor.  Non-factors must then be either character or numeric.}
#'   \item{missing}{the number of \code{NA}s in the variable.}
#'   \item{distinct.values}{the number of distinct values in the variable,
#'   equal to \code{length(table(variable))}.}
#'   \item{precision}{the number of decimal places of precision.}
#'   \item{min}{the minumum value (if numeric) or first value (alphabetically)
#'   as appropriate.}
#'   \item{max}{the maximum value (if numeric) or the last value
#'   (alphabetically) as appropriate.}
#' }
#' @references Special thanks to John Hartigan and the students of
#' 'Statistical Case Studies' of 2004 for their help troubleshooting
#' and developing the function \code{whatis()}.
#' @seealso See also \code{\link{str}}.
#' @author John W. Emerson, Walton Green
#' @examples
#' mydf <- data.frame(a=rnorm(100),
#'                    b=sample(c("Cat", "Dog"), 100, replace=TRUE), 
#'                    c=sample(c("Apple", "Orange", "8"), 100, replace=TRUE),
#'                    d=sample(c("Blue", "Red"), 100, replace=TRUE))
#' mydf$d <- as.character(mydf$d)
#' whatis(mydf)
#'
#' data(iris)
#' whatis(iris)
#' @export
"whatis" <- function(x, var.name.truncate = 20, type.truncate = 14) {
 
  if (!is.data.frame(x)) {
    x <- data.frame(x)
    warning("Object coerced to a data frame.\n")
  }
  if ( is.na(length(dim(x))) | is.null(length(dim(x))) )
    stop("You can not be serious!\n")

  sum.na <- function(y) sum(y, na.rm=TRUE)

  size <- dim(x)
  is.fac <- unlist(lapply(x, is.factor))
  is.char <- unlist(lapply(lapply(x, as.vector), is.character))
  is.fc <- is.fac | is.char
  is.ord <- unlist(lapply(x, is.ordered))

  sum.var.na <- apply(is.na(x), 2, sum)

  is.mixed <- rep(TRUE, size[2])
  this.type <- rep("", size[2])
  num.values <- rep(0, size[2])
  summ <- matrix(NA,size[2],2) 
  summ <- as.data.frame(summ)
  names(summ) <- c("min","max")
  precision <- rep(NA,size[2])
  var.abbrev <- rep('',size[2])
  #uniquifier <- 1

  for (i in 1:(size[2])) {

    if (nchar(names(x)[i]) > var.name.truncate) {
       #var.abbrev[i] <- paste(substr(names(x)[i], 1, 15), '...',
       #                       uniquifier, sep = '')
       #uniquifier <- uniquifier + 1
      var.abbrev[i] <- paste(substr(names(x)[i], 1, var.name.truncate), 
                             '&', sep="")
    } else var.abbrev[i] <- names(x)[i]

    num.values[i] <- length(table(x[,i]))

    z <- as.character(x[,i])
    has.nonnumeric <- (regexpr("[^0-9.-]", gsub(" ", "", z), perl=TRUE) > 0)
    has.numeric <- (regexpr("[0-9.-]", gsub(" ", "", z), perl=TRUE) > 0)
    is.mixed[i] <- sum.na(has.nonnumeric)>0 & sum.na(has.numeric)>0

    if (is.fac[i] & !is.mixed[i]) this.type[i] <- "pure factor"
    if (is.fac[i] & is.mixed[i])  this.type[i] <- "mixed factor"
    if (is.ord[i])                this.type[i] <- "ordered factor"
    if (is.char[i] & !is.fac[i])  this.type[i] <- "character"
    if (!is.char[i] & !is.fac[i]) this.type[i] <- "numeric"

    xtemp <- x[!is.na(x[,i]),i]
    if (length(xtemp)==0) {
      summ[i,] <- rep(NA, 2)
      precision[i] <- NA
    } else {

      if (this.type[i] == "numeric") {
        z <- as.character(x[,i])
        has.decimal <- as.vector(regexpr('.', z, fixed=TRUE))
        has.decimal[is.na(has.decimal)] <- 0
        has.decimal[has.decimal<0] <- 0
        has.decimal[has.decimal>0] <- nchar(z[has.decimal>0]) - has.decimal[has.decimal>0]
        precision[i] <- 10^(-1*max(has.decimal[!is.na(has.decimal)]))
      }

      if (this.type[i]=="numeric")
        summ[i,] <- c(min(x[,i], na.rm=TRUE), max(x[,i], na.rm=TRUE))
      if (this.type[i]!="numeric")
        summ[i,] <- sort(as.character(x[,i]))[c(1,length(xtemp))] 

    }
  }

  summary.by.variable <- data.frame(variable.name=var.abbrev, 
                                    type=substring(this.type, 1, type.truncate),
                                    missing=sum.var.na,
                                    distinct.values=num.values,
                                    precision=precision)
  summary.by.variable <- cbind(summary.by.variable, summ)
  row.names(summary.by.variable) <- 1:ncol(x)
 
  return(summary.by.variable) 

}

