# ************************************************************************
# Main functions for plotting
# 
# Functions: 
#    ggEnsTrend:    (Line)Plot of the average test accuracies for every new classifier added in the ensemble
#    ggEnsHist:     Histogram of the ensemble results
#    ggClassPred:   Barplots of the Correctly Classified Samples
#    ggPermHist:    Histogram of permutation results
#    ggFusedHist:   Fused histogram of permutation and ensemble results
# ************************************************************************


# Plot the average test accuracies for every new classifier added in the ensemble
ggEnsTrend <- function(ensObj, xlabel = NULL, ylabel=NULL, showText=FALSE, xlims=NULL, ylims= NULL){
  .argsCheck(ensObj, "cfBuild")
  
  Ensemble <- AvgAcc <- NULL 
  ensAcc   <- getAcc(ensObj)$Test
  meanVal  <- ensAcc[1]
  ensNum   <- length(getAcc(ensObj)$Test)
  
  # Calculate the average accuracies for every added classifier 
  for (i in 2:length(ensAcc)) {
    meanVal <- c(meanVal, mean(ensAcc[1:i]))
  }
  
  # Store the ensemble num (iter num) and mean accuracies in a data frame 
  avgAcc <- data.frame(1:ensNum, meanVal)
  colnames(avgAcc) <- c("Ensemble","AvgAcc")
  
  # Create the ggplot
  ggLinePlot <- ggplot(data=avgAcc, aes(x=Ensemble, y=AvgAcc)) + 
                geom_point(aes(colour=AvgAcc), size=2.3) + 
                geom_line(linetype="dotted") + 
                theme(legend.position = "none")  
  
  
  if (is.null(xlabel)) { xlabel="Ensemble Iteration" }
  if (is.null(ylabel)) { ylabel="Average Test Accuracy" }
  ggLinePlot <- ggLinePlot + xlab(xlabel) + ylab(ylabel) 
  
  # If showText is TRUE, add the values as text in the plot 
  if (showText == TRUE) { ggLinePlot <- ggLinePlot + geom_text(data=avgAcc, aes(x=Ensemble, y=AvgAcc, label=paste(round(AvgAcc,digits=2),"%",sep="")), vjust=-0.8, size=3.5) }
  
  # Set the limits for the x and y axes if provided
  if (length(xlims) == 2){ ggLinePlot <- ggLinePlot + xlim(xlims) }
  if (length(ylims) == 2){ ggLinePlot <- ggLinePlot + ylim(ylims) }
  
  return(ggLinePlot)
}


# Histogram of the ensemble results 
ggEnsHist <- function (ensObj, density = FALSE, percentiles = FALSE, mean = FALSE, median = FALSE) {
  .argsCheck(ensObj, "cfBuild")
  
  Accuracies <- ..density.. <- ..count.. <- NULL 
  
  # Get the accuracies within the ensemble and the average accuracy
  avgVal <- getAvgAcc(ensObj)$Test
  ensAcc <- data.frame(getAcc(ensObj)$Test)
  colnames(ensAcc) <- "Accuracies"
  
  # Calculate the upper and lower percentile 
  upper <- mean(as.vector(as.matrix(ensAcc))) + 2*sd(as.vector(as.matrix(ensAcc)))
  lower <- mean(as.vector(as.matrix(ensAcc))) - 2*sd(as.vector(as.matrix(ensAcc)))
  
  if (density == TRUE) {
    # Use density
    histPlot <- ggplot(data=ensAcc, aes(x=Accuracies, y=..density..)) + 
                geom_histogram(binwidth=2, colour="#999999", fill="white") +
                geom_density(alpha=.2, fill="white", colour="#333333")  
  } else {
    # Use counts/frequency
    histPlot <- ggplot(data=ensAcc, aes(x=Accuracies, y=..count..)) + 
                geom_histogram(binwidth=2, colour="#999999", fill="white")
  }
  
  # Plot the percentiles 
  if (percentiles == TRUE) {
    histPlot <- histPlot + geom_vline(xintercept=lower, color="purple", linetype="dashed", size=0.8) + 
                           geom_vline(xintercept=upper, color="purple", linetype="dashed", size=0.8)  
  }
  
  # Plot the mean 
  if (mean == TRUE) { histPlot <- histPlot + geom_vline(xintercept=mean(getAcc(ensObj)$Test), color="red", linetype="dashed", size=0.8) }
  
  # Plot the median
  if (median == TRUE) { histPlot <- histPlot + geom_vline(xintercept=median(getAcc(ensObj)$Test), color="cyan", linetype="dashed", size=0.8) }
  
  return(histPlot)
}


# Barplots of the Correctly Classified Samples
ggClassPred <- function(ensObj, position = "stack", displayAll = FALSE, showText=FALSE, xlabel = NULL, ylabel=NULL, cbPalette = FALSE, fillBrewer = FALSE) {
  .argsCheck(ensObj, "cfBuild")
  
  InitClass <- PredClass <- Class <- Percentage <- predictions <- classPlot <- NULL 
  
  # Define a color-blind-friendly palette 
  cbColors  <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Show all classes simultaneously
  if (displayAll == FALSE) {
    predictions <- data.frame(colnames(getConfMatr(ensObj)), diag(getConfMatr(ensObj)))
    colnames(predictions) <- c("Class", "Percentage")
    
    classPlot <- ggplot(data=predictions, aes(x=Class, y=as.numeric(as.vector(Percentage)), fill=Class)) + 
                 geom_bar(width=0.3, stat="identity") + theme_bw() + 
                 xlab("\nClasses") + ylab("Percentages of Correctly Classified Samples per Class (%)\n") +
                 theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_blank()) 
    
    # If showText is TRUE, add the values as text in the plot
    if (showText == TRUE) { 
      classPlot <- classPlot + geom_text(data=predictions, aes(x = Class, y = as.numeric(as.vector(Percentage)), label=paste(as.vector(Percentage),"%",sep="")), vjust=-0.5, size=3.7) 
    } 
  } else {
    predictions <- as.data.frame(ftable(getConfMatr(ensObj), row.vars = 2:1))
    colnames(predictions) <- c("PredClass", "InitClass", "Percentage")
    
    classPlot <- ggplot(data = predictions, aes(x=InitClass, y=as.numeric(as.vector(Percentage)), fill=PredClass, ymin=0, ymax=100)) + 
                 geom_bar(width=0.3, stat="identity", position=position) + theme_bw() + 
                 theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) 
    
    # If showText is TRUE, add the values as text in the plot
    if (showText == TRUE) {
      if (position == "stack") {
        # Position stack
        classPlot <- classPlot + geom_text(data=predictions, aes(x=InitClass, y=as.numeric(as.vector(Percentage)), label=ifelse(Percentage>0, paste(as.vector(Percentage),"%",sep=""), "")), size=3.7, position="stack") 
      } else {
        # Position dodge
        classPlot <- classPlot + geom_text(data=predictions, aes(x=InitClass, y=as.numeric(as.vector(Percentage)), label=paste(as.vector(Percentage),"%",sep="")), size=3.7, position=position_dodge(width=0.5))
      }
    }
  }
  
  # If color-blind-friendly palette selected
  if (cbPalette == TRUE)  { classPlot <- classPlot + scale_fill_manual(values=cbColors) }
  
  # If brewer palette selected
  if (fillBrewer == TRUE) { classPlot <- classPlot + scale_fill_brewer() }
  
  # If no x and y axes titles provided, use the default
  if (is.null(xlabel)) { xlabel="\nClasses" }
  if (is.null(ylabel)) { ylabel="Percentages of Class Predictions (%)\n" }
  classPlot <- classPlot + xlab(xlabel) + ylab(ylabel)
  
  return(classPlot)
}

# Histogram of permutation results 
ggPermHist <- function(permObj, density = FALSE, percentiles = FALSE, mean = FALSE, median = FALSE) {
  .argsCheck(permObj, "cfPermute")
  
  Accuracies <- x <- y <- ..density.. <- ..count.. <- NULL 
  
  # Get the permutation results and store in a data frame
  permAcc <- as.data.frame(permObj$avgAcc)
  colnames(permAcc) <- "Accuracies"
  
  # Calculate the upper and lower percentile
  upper <- mean(as.vector(as.matrix(permAcc))) + 2*sd(as.vector(as.matrix(permAcc)))
  lower <- mean(as.vector(as.matrix(permAcc))) - 2*sd(as.vector(as.matrix(permAcc)))
  
  if (density == TRUE){
    # Use density
    permPlot <- ggplot(data=permAcc, aes(x=Accuracies, y=..density..))+ 
                geom_histogram(binwidth=2, colour="#999999", fill="white") + 
                geom_density(alpha=.2, fill="white", colour="#333333")  
  } else {
    # Use frequency/counts
    permPlot <- ggplot(data=permAcc, aes(x=Accuracies, y=..count..)) + 
                geom_histogram(binwidth=2, colour="#999999", fill="white")
  }
  
  # Plot the percentiles 
  if (percentiles == TRUE){
    permPlot <- permPlot + geom_vline(xintercept=lower, color="purple", linetype="dashed", size=0.8) + 
                           geom_vline(xintercept=upper, color="purple", linetype="dashed", size=0.8)  
  }
  
  # Plot the mean 
  if (mean == TRUE) { permPlot <- permPlot + geom_vline(xintercept=mean(permObj$avgAcc), color="red", linetype="dashed", size=0.8) } 
  
  # Plot the median 
  if (median == TRUE) { permPlot <- permPlot + geom_vline(xintercept=median(permObj$avgAcc), color="cyan", linetype="dashed", size=0.8) }
  
  return (permPlot)
}


ggFusedHist <- function(ensObj, permObj) {
  .argsCheck(ensObj, "cfBuild")
  .argsCheck(permObj, "cfPermute")
  
  acc <- type <- NULL
  
  fusedMatr <- rbind( cbind(rep("permutation", length(permObj$avgAcc)), permObj$avgAcc), cbind(rep("ensemble", length(ensObj$testAcc)), ensObj$testAcc))
  fusedMatr <- as.data.frame(fusedMatr)
  colnames(fusedMatr) <- c("type", "acc")
  
  ggHist <- ggplot(data=fusedMatr, aes(x=as.numeric(as.vector(acc)), fill=type, colour=type)) + 
            geom_histogram(binwidth=2, alpha=.6, colour="#999999") + theme_bw() + scale_fill_manual(values = c("red", "blue")) + 
            xlab("Overall Test Accuracies (%CC)") + ylab("Frequency\n") + theme(legend.position = "none") 
  
  return (ggHist)
}

