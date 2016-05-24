#' Test logical expressions
#'
#' Check if logical expressions return what you expect with a truth table
#' 
#' @details This is a nice way to check operator precedence, see \code{\link{Syntax}}
#'
#' @return Truth table as data.frame with TRUE and FALSE (and NA) combinations
#' @author Berry Boessenkool, \email{berry-b@@gmx.de},  Mrz 2016 
#' @seealso \code{\link{logical}}
#' @keywords logic
#' @export
#' @examples
#' TFtest(!a & !b)
#' TFtest(!a & !b, a&b, !(a&b))
#' TFtest(!a & !b | c)
#' TFtest(!a & !b | c, na=FALSE)
#' TFtest(!a)
#' TFtest(a&b|c, (a&b)|c, a&(b|c), na=FALSE) # AND has precedence over OR
#'
#' @param \dots Expression(s) with logical operators to be evaluated, 
#'        with single letters for variables. Each expression is to be separated with a comma
#' @param na Logical: should NAs be included in the truth table? DEFAULT: TRUE
#' 
TFtest <- function(
...,
na=TRUE
)
{
# expressions as character strings
depsub <- as.list(substitute(list(...))[-1])
depsub <- sapply(depsub, function(x) deparse(x))

# letters in expression
lets <- gsub("[^[:alpha:]]+", "", depsub)
lets <- paste(lets, collapse="")
nn <- length(unique(strsplit(lets,"")[[1]]))

# combinations of T and F in truth table:
TFvalues <- c(TRUE, FALSE, if(na) NA)
out <- expand.grid(rep(list(TFvalues),nn))[,nn:1, drop=FALSE]
colnames(out) <- letters[1:nn]

# create objects within function:
for(i in 1:nn) assign(letters[i], out[,i])

# evaluate the expression(s):
result <- lapply(depsub, function(x) eval(parse(text=x)))
result <- as.data.frame(result)
colnames(result) <- paste("__", depsub)

# return output:
cbind(out, result)
}
