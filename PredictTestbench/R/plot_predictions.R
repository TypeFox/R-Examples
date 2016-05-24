#' Function to plot the Error Comparison
#'
#' @param  dataIn as input data in list format (returned by function prediction_errors())
#' @return It returns the Error comparison for different methods
#' @importFrom  reshape2 melt
#' @export
#' @examples
#' # aa <- prediction_errors()
#' # bb <- plot_predictions(aa)
#' # bb
plot_predictions <- function(dataIn)
{
  Predicted_Time_Series <- NULL
  Predicted_Values <- NULL
  Methods <- NULL

  q1 <-  as.character(dataIn[[1]])

list.condition <- sapply(dataIn, function(x) length(x)>1)
dataIn  <- dataIn[list.condition]


  q <- data.frame(dataIn)
  Missing_Percent <- 1:length(dataIn[[1]])
  mm <- length(q)
  q[,(mm+1)] <- data.frame(Missing_Percent)

  melted <- melt(q, id.vars = "Missing_Percent")
  colnames(melted) <- c("Predicted_Time_Series","Methods","Predicted_Values")
  d <- ggplot(data=melted, aes(x=Predicted_Time_Series, y=Predicted_Values, group=Methods, color=Methods)) + labs(title = q1) + geom_line()
  return(d)
}
