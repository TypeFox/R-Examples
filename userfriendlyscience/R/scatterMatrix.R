scatterMatrix <- function(dat, items=NULL, plotSize=180, sizeMultiplier = 1,
                          axisLabels = "none", powerHist=TRUE, ...) {

  if (is.null(items)) {
    items <- names(dat);
  }
  
  ### Generate object with 3 sub-objects to store input,
  ### intermediate results, and output
  res <- list(input = list(dat=dat,
                           items=items,
                           plotSize=plotSize,
                           sizeMultiplier=sizeMultiplier,
                           axisLabels=axisLabels),
              intermediate = list(),
              output = list());
  
  ### Extract dataframe and select only complete cases
  res$intermediate$dat <- dat[complete.cases(dat[, items]), items];

  ### Convert all variables to numeric vectors, if they weren't already
  res$intermediate$dat <- data.frame(lapply(res$intermediate$dat, 'as.numeric'));
  
  ### The size of each panel in the scattermatrix depends
  ### on the number of items - therefore, we need to adjust
  ### the plot sizes to the number of items.
  res$intermediate$baseSize <- baseSize <-
    (sizeMultiplier * (plotSize / length(items))) / 100;
  
  res$intermediate$plotSettings <- plotSettings <-
    theme(axis.line = element_line(size = baseSize),
          panel.grid.major = element_line(size = baseSize/2),
          line = element_line(size = baseSize/2),
          axis.ticks = element_line (size=baseSize/2)
    );  
  
  ### Visual representation of bivariate correlations
  ### First generate a normal scattermatrix with histograms
  ### on the diagonal
  res$intermediate$ggpairs.normal <-
    ggpairs(res$intermediate$dat, diag=list(continuous="barDiag", discrete="barDiag"),
            axisLabels=res$input$axisLabels);
#  lower="blank",
  
  ### Then generate one with jittered points
  res$intermediate$ggpairs.jittered <-
    ggpairs(res$intermediate$dat,
            diag=list(continuous="blankDiag"),
            upper=list(continuous=GGally::wrap("cor")),
            lower=list(continuous=GGally::wrap("points", position="jitter")),
            axisLabels=res$input$axisLabels);

  ### Copy the the one with the jittered points
  res$intermediate$ggpairs.combined <- res$intermediate$ggpairs.jittered;

  if (powerHist) {
    ### Create histograms and add them to the combined plot
    res$intermediate$powerHists <- list();
    for (currentVar in 1:length(items)) {
      res$intermediate$powerHists[[items[currentVar]]] <-
        powerHist(res$intermediate$dat[[items[currentVar]]], ...);
      res$intermediate$ggpairs.combined <-
        putPlot(res$intermediate$ggpairs.combined,
                res$intermediate$powerHists[[items[currentVar]]]$plot,
                currentVar, currentVar);
    }
  }
  else {
    ### Then place the histograms from the 'normal' one
    ### on the diagonal of the jittered scattermatrix
    for (currentVar in 1:length(items)) {
      res$intermediate$ggpairs.combined <-
        putPlot(res$intermediate$ggpairs.combined,
                getPlot(res$intermediate$ggpairs.normal, currentVar, currentVar),
                currentVar, currentVar);
    }
  }
  
  ### Copy combined matrix to the output for final adjustments
  res$output$scatterMatrix <- res$intermediate$ggpairs.combined;
  
  ### Adjust the size of the plots
  for (currentRowFromTop in 1:length(items)) {
    for (currentColumnFromLeft in 1:length(items)) {
      res$output$scatterMatrix <-
        putPlot(res$output$scatterMatrix,
                getPlot(res$output$scatterMatrix, currentRowFromTop, currentColumnFromLeft) + plotSettings,
                currentRowFromTop, currentColumnFromLeft);
    }
  }

  ### Set class and return result
  class(res) <- "scatterMatrix";
  return(res);
  
}

print.scatterMatrix <- function(x, ...) {
  ###
  print(x$output$scatterMatrix, ...);
}