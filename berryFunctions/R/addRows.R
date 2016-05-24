#' Add n rows to a data.frame
#' 
#' simple Helper-Function to add n rows to a data.frame.
#' 
#' @param df Dataframe object
#' @param n Number of rows to add
#' @param values Values to be used in the new rows. DEFAULT: NA

#' @return A data.frame
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2014
#' @seealso \code{\link{insertRows}}, \code{\link{data.frame}}, \code{\link{matrix}}, \code{\link{rbind}}
#' @keywords misc
#' @export
#' @examples
#' 
#' MYDF <- data.frame(A=5:3, B=2:4)
#' addRows(MYDF, 3)
#' 
addRows <- function(
          df,
          n,
          values=NA)
{
dnew <- data.frame(matrix(values, nrow=n, ncol=ncol(df)))
colnames(dnew) <- colnames(df)
rbind(df, dnew)
}
