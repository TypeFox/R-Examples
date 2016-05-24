#' Creating a data frame for an aoristic graph using all data
#' @param data data.frame created by aoristic.df
#' @return data data.frame with two columns (hour, freq) used to create an
#'   aoristic graph for the entire study area
#' @import lubridate
#' @export
#' @examples
#' \donttest{
#' data(aoristic)
#' data <- aoristic.df(data=arlington, 
#'          DateTimeFrom="DateTimeFrom", DateTimeTo="DateTimeTo")
#' graph <- aoristic.all.graph(data=data)
#' ggplot(graph, aes(x=hour, y=freq)) + 
#'    geom_bar(stat="identity") + 
#'    ggtitle("Aoristic Graph for the Entire Study Area")
#' 
#' # with probability labels
#' graph$prob <- paste(round(graph$freq / sum(graph$freq) * 100, 1), "%", sep="")
#' ggplot(graph, aes(x=hour, y=freq)) + 
#'    geom_bar(stat="identity") + 
#'    ggtitle("Aoristic Graph for the Entire Study Area") + 
#'    geom_text(aes(y=freq, label=prob), vjust=1.5, colour="white", size=4)
#'       
#' }
aoristic.all.graph <- function(data){
  
  #defining variables (to avoid "Note" in the package creation)
  data.temp <- melt(data, id.vars=c("id", "HourFrom", "duration"))
  data.temp$variable <- as.numeric(gsub("time", "", data.temp$variable))

  graph <- aggregate(data.temp$value, by=list(data.temp$variable), FUN=sum)
  names(graph) <- c("hour", "freq")
  graph$hour <- as.numeric(graph$hour) 

  # convert as factor and reorder for a graph
  graph$hour <- factor(graph$hour, 
                     levels=c("6", "7", "8", "9", "10", "11", "12", 
                              "13", "14", "15", "16", "17", "18", "19", "20",
                              "21", "22", "23", "0", "1", "2", "3", "4", "5"))

  graph <- graph[order(graph$hour),]
  return(graph)
}