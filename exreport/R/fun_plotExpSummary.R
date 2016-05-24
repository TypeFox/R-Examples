#' Barplot for summarizing an experiment output variable
#'
#' This function builds a barplot for a given experiment output variable,
#' summarizing its distribution according to the different methods and problems.
#' The aspect of the plot can be parametrized in several ways.
#'
#' Please notice that the plot function requires that an unique configuration
#' of parameters is present in the experiment. So the user must have processed and
#' instantiated the experiment before.
#'
#' @export
#' @import ggplot2
#' @param exp - The experiment object to take the data from
#' @param output - A string identifying the name of the output variable to be
#' plotted
#' @param columns - Integer number, 0 for a wide aspect plot and any other value
#' to include n columns of facets separating the problems.
#' @param freeScale - Boolean, if using facets sets the scale of each one 
#' independent or not
#' @param fun - A function to be applied to the selected output variables before
#' being plotted.
#' @param grayscale - Defaulted to False. True for a plot in grayscale.
#' @return an exPlot object
#' 
#' 
#' @examples
#' # This example plots the distribution of the trainingTime variable in the 
#' # wekaExperiment problem.
#' 
#' # First we create the experiment from the problem.
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' 
#' # Next we must process it to have an unique parameter configuration:
#' # We select a value for the parameter featureSelection:
#' experiment <- expSubset(experiment, list(featureSelection = "yes"))
#' # Then we reduce the fold parameter:
#' experiment <- expReduce(experiment, "fold", mean)
#' # Finally we remove unary parameters by instantiation:
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Now we can generate several plots:
#' 
#' # Default plot:
#' plotExpSummary(experiment, "accuracy")
#' 
#' # We can include faceting in the plot by dividing it into columns:
#' plotExpSummary(experiment, "accuracy", columns=3)
#' 
#' # If we want to show the independent interaction for the output variable
#' # in each experiment we can make the scales for example, remark the difference
#' # in :
#' plotExpSummary(experiment, "trainingTime", columns=3, freeScale=FALSE)
#' plotExpSummary(experiment, "trainingTime", columns=3, freeScale=TRUE)
#' 
plotExpSummary <- function(exp, output, columns=0, freeScale=FALSE, fun=identity, grayscale=FALSE) {
  
  if( !is.experiment(exp) )
    stop(.callErrorMessage("wrongParameterError", "exp", "experiment"))
  
  if( !is.character(output) )
    stop(.callErrorMessage("wrongParameterError", "output", "character"))
  
  if (!all(output %in% exp$outputs))
    stop(.callErrorMessage("variableNotPresentError", output))  
  
  # Check if instantiation of parameters is needed
  if( length(exp$parameters) != 0 )
    stop(.callErrorMessage("requireInstantiationError"))
  
  data <- exp$data
  
  if (freeScale)
    scales<-"free_y"
  else
    scales<-"fixed"
  
  if (columns > 0){
    .env <- environment()
    expr <- "p <- ggplot(data, aes(y=%s, x=%s, fill=%s), environment = .env)"
    expr <- sprintf(expr,output,exp$method,exp$method)
    eval(parse(text=expr))
  }
  else{
    .env <- environment()
    expr <- "p <- ggplot(data, aes(y=%s, x=%s, group=%s, fill=%s), environment = .env)"
    expr <- sprintf(expr,output,exp$problem,exp$method,exp$method)
    eval(parse(text=expr))
  }
  
  p <- p + geom_bar(stat="identity", colour="black", position=position_dodge())
  
  if (columns > 0){
    expr <- "p <- p + facet_wrap(~%s, ncol=columns, scales=scales)"
    expr <- sprintf(expr,exp$problem)
    eval(parse(text=expr))
  }
  
  p <- p + theme(axis.text.x  = element_text(face="bold", colour="#000000", size=12, angle= 90))
  p <- p + theme(axis.title.y = element_text(face="bold", colour="#000000", size=20, vjust=0.5))
  p <- p + ylab(output)
  p <- p + xlab("")
  
  if (grayscale)
    p <- p + scale_fill_grey()
    
  title <- sprintf("Results for output \"%s\"", output)
  
  res <- .exPlot(p, target = toString(output), title = title, alias = "ExpSummaryPlot", tags = exp$tags)
  res
}
