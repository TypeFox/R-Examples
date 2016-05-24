#' Area plot for the rank distribution from a multiple test
#'
#' This function builds an area plot from a testMultiple object displaying the
#' cumulative value for each method for all the evaluated problems. The value
#' for the rankings is obtained from the Friedman test independently of the scope
#' of the test (control or pairwise).
#'
#' @export
#' @import ggplot2
#' @param testMultiple Statistical test from which the plot is generated. The 
#' rankings are obtained from the Friedman test.
#' @param grayscale Configure the plot using a grayscale palette.
#' @return an exPlot object
#' 
#' 
#' @examples
#' # First we create an experiment from the wekaExperiment problem and prepare
#' # it to apply the test:
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' experiment <- expReduce(experiment, "fold", mean)
#' experiment <- expSubset(experiment, list(featureSelection = "no"))
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Then we perform a Friedman test included ina a testMultipleControl 
#' # test procedure
#' test <- testMultipleControl(experiment, "accuracy")
#' 
#' # Finally we obtain the plot
#' plotCumulativeRank(test)
#' cat()
plotCumulativeRank <- function(testMultiple, grayscale=FALSE) {
  
  if( !is.testMultiple(testMultiple) )
    stop(.callErrorMessage("wrongParameterError", "testMultiple", "testMultiple"))
  
  friedman <- testMultiple$friedman
  
  da <- friedman$ranks
  da <- reshape2::melt(da)
  da$Var1 <- as.factor(da$Var1)
  da$Var2 <- as.factor(da$Var2)
  # Now this magic computes the method ordering for the mean ranking and updates
  # ther order of the factor levels to sort the plot
  means <- aggregate(da$value, by=list(da$Var1), FUN=mean)
  means$rank <- rank(means$x)
  means <- means[with(means, order(rank)),]
  da$Var1 <- factor(da$Var1, means$Group.1)
  # A continuos variable is needed for the plot scale, we will change the scale further
  da$continuous <- as.numeric(da$Var1)
  problemName <- testMultiple$experiment$problem
  colnames(da)[2] <- problemName # Just to print the legend title correctly named.
  
  p <- eval(parse(text = "ggplot(da, aes(x=continuous, y=value))"))
  
  expr <- "p <- p + geom_area(aes(colour = %s, fill= %s), position = 'stack')"
  expr <- sprintf(expr,problemName,problemName)
  eval(parse(text=expr))
  
  p <- p + theme(axis.text.x  = element_text(face="bold", colour="#000000", size=12, angle= 90))
  p <- p + scale_x_discrete(breaks=da$continuous, labels=da$Var1)
  p <- p + theme(axis.text.x  = element_text(face="bold", colour="#000000", size=12, angle= 90, hjust = 1))
  p <- p + theme(axis.title.y = element_text(face="bold", colour="#000000", size=20, vjust=0.5))
  p <- p + ylab(sprintf("Cumulative Ranking for output \"%s\"", friedman$tags$target))
  p <- p + xlab("")
  
  if (grayscale)
    p <- p + scale_fill_grey()
  
  title <- sprintf("Cumulative Ranking for Var %s", friedman$tags$target)
  
  res <- .exPlot(p, title = title, alias = "CumulativeRank", tags = testMultiple$tags)
  res
}