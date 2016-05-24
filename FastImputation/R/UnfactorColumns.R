#' Convert columns of a dataframe from factors to character or numeric.
#'
#' @param x A dataframe
#' @return A dataframe containing the same data but any \code{factor} columns have been replaced with numeric or character columns.
#' @export
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
UnfactorColumns <- 
function(x) 
{
  # test.df <- data.frame(num.factor=as.factor(c(1,1,3,4,6,6)), char.factor=as.factor(c("this", "this", "that", "this", "this", "that")))
  if( !is.data.frame(x) ) stop("UnfactorColumns requires a data.frame for input")
  
  for(ic in 1:ncol(x)) {
    if(is.factor(x[,names(x)[ic]])) {
      # Unfactor the column
      x[,names(x)[ic]] <- levels(x[,names(x)[ic]])[as.numeric(x[,names(x)[ic]])]
      
      # Change to numeric if appropriate
      old.warn <- getOption("warn")
      options(warn = -1)
      if( 0==sum(is.na(as.numeric(x[,names(x)[ic]]))) ) x[,names(x)[ic]] <- as.numeric(x[,names(x)[ic]])
      options(warn = old.warn)
    }
  }
  return(x)
}
