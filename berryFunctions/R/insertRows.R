#' insert rows to data.frame
#' 
#' Insert (multiple) rows to a data.frame, possibly coming from another data.frame, with value and row recycling
#' 
#' @return data.frame
#' @note Has not yet been tested with RWI (really weird input), so might not be absolutely foolproof
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Oct 2015, based on code by Ari B. Friedmann 
#'         (I added the for loop, recycling, input controls and data.framification added by)
#' @seealso \code{\link{addRows}}
#' @references \url{http://stackoverflow.com/questions/11561856/add-new-row-to-dataframe}
#' @keywords misc manip array
#' @export
#' @examples
#' 
#' existingDF <- as.data.frame(matrix(1:20, nrow=5, ncol=4))
#' existingDF
#' insertRows(existingDF, 2) # default new=NA is recycled
#' insertRows(existingDF, 2, 3:6)
#' insertRows(existingDF, 3, new=matrix(10:1,ncol=2)) # input warning
#' 
#' # Works for multiple rows as well:
#' insertRows(existingDF, r=c(2,4,5), new=NA)
#' 
#' # Also works with a data.frame for insertion:
#' insertDF <- as.data.frame(matrix(101:112, nrow=3, ncol=4))
#' insertRows(existingDF, 3, new=insertDF) # excess rows in new are ignored
#' insertRows(existingDF, c(2,4,5), new=insertDF)
#' insertRows(existingDF, c(2,4:6), new=insertDF) # rows are recycled
#' 
#' @param df data.frame
#' @param r Row number (not name!), at which the \code{new} row is to be inserted. Can be a vector
#' @param new Vector with data to be inserted, is recycled. Alternatively, a data.frame, whose rows are put into the r locations. 
#'        If it has more rows than length(r), the excess rows are ignored. DEFAULT: NA
#' 
insertRows <- function(
df,
r,
new=NA
)
{
# Input checks:
if(!is.data.frame(df)) warning("df is not a data.frame.")
if(!(is.vector(new) | is.data.frame(new))) warning("new row is not a vector or data.frame.")
# recycle new row values:
if(!is.data.frame(new))
  {
  new <- rep(new, length=ncol(df)) # recycle vector to number of columns
  new <- rbind(new) # make it into a data.frame (well, a matrix, but that works fine, too)
  }
# recycle the rows:
new <- new[rep(1:nrow(new), length=length(r)), , drop=FALSE]
# for each value in r, create a row with the i-th row of new:
for(i in 1:length(r))
  {
  SEQ <- seq(from=r[i], to=nrow(df))
  df[SEQ+1,] <- df[SEQ,]
  df[r[i],] <- new[i,]
  }
# return the final output:
df
}
