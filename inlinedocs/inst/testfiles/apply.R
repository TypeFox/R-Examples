apply <- function # Apply Functions Over Array Margins
### Returns a vector or array or list of values obtained by applying a
### function to margins of an array or matrix.
(X, ##<< an array, including a matrix.
 MARGIN,
### the margin.
 FUN, 
### the function to be applied: see \sQuote{Details}. In the case of
### functions like \code{+}, \code{%*%}, etc., the function name
### must be backquoted or quoted.
 ... ##<< optional arguments to \code{FUN}.
 ){
  5
### test apply.
}

.result <-
  list(apply=list(title="Apply Functions Over Array Margins",
         description=paste("Returns a vector or array or list of values obtained by applying a",
           "function to margins of an array or matrix.",sep="\n"),
         "item{X}"="an array, including a matrix.",
         "item{MARGIN}"="the margin.",
         "item{FUN}"=paste("the function to be applied: see \\sQuote{Details}. In the case of",
           "functions like \\code{+}, \\code{%*%}, etc., the function name",
           "must be backquoted or quoted.",sep="\n"),
         "item{\\dots}"="optional arguments to \\code{FUN}.",
         value="test apply.",format="",
         definition="apply <- function # Apply Functions Over Array Margins\n### Returns a vector or array or list of values obtained by applying a\n### function to margins of an array or matrix.\n(X, ##<< an array, including a matrix.\n MARGIN,\n### the margin.\n FUN, \n### the function to be applied: see \\sQuote{Details}. In the case of\n### functions like \\code{+}, \\code{%*%}, etc., the function name\n### must be backquoted or quoted.\n ... ##<< optional arguments to \\code{FUN}.\n ){\n  5\n### test apply.\n}"))

