### Note: this is necessary to prevent Rcmd CHECK from throwing a note;
### otherwise it think these variable weren't defined yet.
utils::globalVariables(c('distribution', '..density..', 'normalX', 'normalY'));

### Make a 'better' histogram using ggplot
powerHist <- function(vector,
                      histColor = "#0000CC",
                      distributionColor = "#0000CC",
                      normalColor = "#00CC00",
                      distributionLineSize = 1,
                      normalLineSize = 1,
                      histAlpha = .25,
                      xLabel = NULL,
                      yLabel = NULL, density=FALSE,
                      theme=dlvTheme(axis.title=element_text(colour = "black")),
                      rug=TRUE, jitteredRug=TRUE, rugSides="b",
                      rugAlpha = .2) {
  varName <- deparse(substitute(vector));
  vector <- na.omit(vector);
  if (!is.numeric(vector)) {
    tryCatch(vector <- as.numeric(vector), error = function(e) {
      stop("The vector you supplied is not numeric; I tried to convert it, ",
           "but my attempt failed. The error I got is:\n", e);
      });
  }
  res <- list(input = as.list(environment()), intermediate = list(), output = list());
  res$input$sampleSize = length(vector);
  res$intermediate$normalX <- c(seq(min(res$input$vector), max(res$input$vector),
                     by=(max(res$input$vector) -
                         min(res$input$vector)) /
                         (res$input$sampleSize-1)));
  res$intermediate$normalY <- dnorm(res$intermediate$normalX,
                                    mean=mean(res$input$vector),
                                    sd=sd(res$input$vector));
  res$intermediate$distribution <- res$input$vector;
  res$dat <- data.frame(normalX = res$intermediate$normalX,
                        normalY = res$intermediate$normalY,
                        distribution = res$intermediate$distribution);
  res$intermediate$tempBinWidth <- (max(res$input$vector) -
                                    min(res$input$vector)) / 30;
  ### Generate labels if these weren't specified
  if (is.null(xLabel)) {
    xLabel <- paste0('Value of ', varName);
  }
  if (is.null(yLabel)) {
    yLabel <- ifelse(density, "Density", "Frequency");
  }
  if (density) {
    ### Plot distribution
    res$plot <- ggplot(data=res$dat, aes(x=distribution)) + 
      xlab(xLabel) +
      ylab(yLabel) +
      geom_histogram(aes(y=..density..), color=NA, fill=histColor,
                     alpha=histAlpha, binwidth=res$intermediate$tempBinWidth) +
      geom_density(color=distributionColor, size=distributionLineSize) +
      geom_line(aes(x=normalX, y=normalY), color=normalColor, size=normalLineSize) +
      theme;
  }
  else {
    ### Get approximate maximum frequency in histogram and multiply
    ### normal curve and density curve
    res$dat$density <- density(vector, n=length(vector),
                               from=min(vector),
                               to=max(vector),
                               adjust=.8)$y;
    lineScaling <- max(summary(cut(vector, breaks=30))) / max(res$dat$density);
    res$dat$density <- lineScaling * res$dat$density;
    res$dat$normalY <- lineScaling * res$dat$normalY;

    ### Plot distribution
    res$plot <- ggplot(data=res$dat, aes(x=distribution)) + 
      xlab(xLabel) +
      ylab(yLabel) +
      geom_histogram(color=NA, fill=histColor,
                     alpha=histAlpha, binwidth=res$intermediate$tempBinWidth) +
      geom_line(aes(x=normalX, y=density), color=distributionColor, size=distributionLineSize) +
      geom_line(aes(x=normalX, y=normalY), color=normalColor, size=normalLineSize) +
      theme;
  }
  if (rug) {
    if (jitteredRug) {
      res$plot <- res$plot + geom_rug(color=distributionColor, sides=rugSides,
                                      aes(y=0), position="jitter",
                                      alpha=rugAlpha);
    } else {
      res$plot <- res$plot + geom_rug(color=distributionColor, sides=rugSides,
                                      alpha=rugAlpha);
    }
  }
  if (!is.null(res$input$xLabel) && is.logical(res$input$xLabel) &&
        !(res$input$xLabel)) {
    res$plot <- res$plot + theme(axis.title.x = element_blank());
  }
  if (!is.null(res$input$yLabel) && is.logical(res$input$yLabel) &&
        !(res$input$yLabel)) {
    res$plot <- res$plot + theme(axis.title.y = element_blank());
  }
  ### Set class and return result
  class(res) <- "powerHist";
  return(res);
}

print.powerHist <- function(x, ...) {
  print(x$plot, ...);
}
