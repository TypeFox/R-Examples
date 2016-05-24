#' Round numeric columns of a data frame
#' 
#' Round the numeric columns of a data frame to a specified number of significant digits.
#' 
#' @author Stephen Turner
#' 
#' @param df A data.frame.
#' @param n The number of significant digits to round off to.
#' 
#' @return A data.frame rounded to n significant digits.
#' 
#' @examples
#' dfsigfig(mtcars,1)
#' 
#' @export
dfsigfig <- function(df, n=3) {
    df[ ,sapply(df, is.numeric)] <- signif(df[ ,sapply(df, is.numeric)], n)
    df
}
