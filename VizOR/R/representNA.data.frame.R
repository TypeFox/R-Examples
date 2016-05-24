##' Apply the \code{representNA} function to all columns of a data frame
##' 
##' Replace the string "NA" with R's missing value NA, wherever it appears in
##' the given data frame.
##' 
##' @param df The data frame to be 'cleaned up'.
##' @return Returns a data frame identical to df, with "NA" recoded as NA.
##' @author David C. Norris
##' @keywords manip
##' @examples # TODO: Provide an example
##' @export representNA.data.frame
representNA.data.frame <-
function(df){
  for(col in dim(df)[2])
    df[,col] <- representNA(df[,col])
}

