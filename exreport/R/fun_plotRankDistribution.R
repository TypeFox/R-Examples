#' Boxplot for the ranks distribution and control hypotheses from multiple test
#'
#' This function generates a boxplot from a testMultiple statistical test
#' showing the ordered distrubution of rankings for each method computed for
#' the Friedman test. If the input test features a control multiple comparison
#' then the rejected hypotheses by the Holm methd are also indicates in the plot.
#'
#' @export
#' @import ggplot2
#' @param testMultiple The statistical test from which the plot is generated.
#' The functions accepts either control and pairwise multiple tests.
#' @return an experiment object
#' 
#' @examples
#' # First we create an experiment from the wekaExperiment problem and prepare
#' # it to apply the test:
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' experiment <- expReduce(experiment, "fold", mean)
#' experiment <- expSubset(experiment, list(featureSelection = "yes"))
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Then we perform a Friedman test included ina a testMultipleControl 
#' # test procedure
#' test <- testMultipleControl(experiment, "accuracy")
#' 
#' # Finally we obtain the plot
#' plotRankDistribution(test)
#' cat()
plotRankDistribution <- function(testMultiple) {
  
  if( !is.testMultiple(testMultiple) )
    stop(.callErrorMessage("wrongParameterError", "testMultiple", "testMultiple"))
  
  da <- data.frame(method = rownames(testMultiple$friedman$ranks), ranks = testMultiple$friedman$ranks)
  da <- reshape2::melt(da, id.vars="method")
  
  means <- aggregate(value ~ method, da, mean)
  
  # The boxes are shaded only if a control post-hoc has been peformed
  if (is.testMultipleControl(testMultiple)) {
    accepted <- testMultiple$pvalues > testMultiple$tags$alpha
    color <- accepted
    color[accepted==TRUE | is.na(accepted)]  <- "white"
    color[accepted==FALSE] <- "grey"
    # Now color says if the method i in testMultiple$names pass the test or not
    # (but the index is different from da!) --> We have to reorder
    disorderedNames <- as.character(means[["method"]]) # as characters
    orderedNames <- as.character(testMultiple$names[[1]]) # as characters
    indexNames <- sapply(orderedNames,FUN=function(x) which(disorderedNames==x))
    #Now we reorder the color array
    color <- color[indexNames]
  } else {
    color <- rep("white", nrow(means))
  }

  means <- cbind(means, color)
  means <- means[order(means$value),,drop=FALSE]
  
  p <- eval(parse(text = "ggplot(da, aes(x=reorder(method, value), value))"))
  p <- p + geom_boxplot(fill = as.character(means$color))
  p <- p + eval(parse(text = "geom_text(data = means, aes(label = sprintf(\"%.2f\",value), y = value + 0.08), colour=\"black\", fontface=\"bold\", size=4)"))
  p <- p + ylab(sprintf("Ranking Distributions for var %s", testMultiple$tags$target))
  p <- p + xlab("")
  p <- p + theme(axis.text.x  = element_text(face="bold", colour="#000000", size=12, angle= 90, hjust = 1))
  p <- p + theme(axis.title.y = element_text(face="bold", colour="#000000", size=20, vjust=0.5))
  p <- p + guides(colour=FALSE)
  
  title <- sprintf("Distribution of ranks for output \"%s\"", testMultiple$tags$target)
  
  res <- .exPlot(p, title = title, alias = "RankDistributionPlot", tags = testMultiple$tags)
  res
}