#' Function to plot the Error Comparison
#'
#' @param  dataIn as input data in list format (returned by function inpute_errors())
#' @return It returns the Error comparison for different methods
#' @importFrom  reshape2 melt
#' @export
#' @examples
#' # aa <- impute_errors()
#' # bb <- plot_errors(aa)
#' # bb
plot_errors <- function(dataIn)
{
    #melt <- NULL
    Percent_of_Missing_Values <- NULL
    Error_Values <- NULL
    Methods <- NULL

    q1 <-  as.character(dataIn[[1]])
    q <- data.frame(dataIn[-1])
    #q1 <- dataIn$Parameter

  melted <- melt(q, id.vars = "Missing_Percent")
  colnames(melted) <- c("Percent_of_Missing_Values","Methods","Error_Values")
  d <- ggplot(data=melted, aes(x=Percent_of_Missing_Values, y=Error_Values, group=Methods, color=Methods)) + labs(title = q1) + geom_line()
  return(d)
}
