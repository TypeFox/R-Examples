##' @title List to DataFrame
##' @description converts a list in a dataframe
##' @param dlist a list
##' @examples
##' a <- 1:5
##' listToDataFrame(a)
##' b = listToDataFrame(a)
##' @export
listToDataFrame <- function (dlist) {
  n = length(dlist)
  if (n>0) {
    df = data.frame(dlist[1])
    for (i in 2:n) {
      df = rbind(df, dlist[[i]])
    }
  } else {
    stop("Enter with a non empty list")
  }
  return(df)
}
