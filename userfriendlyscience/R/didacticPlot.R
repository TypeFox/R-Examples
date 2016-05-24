
### Note: this is necessary to prevent Rcmd CHECK from throwing a note;
### otherwise it think these variables weren't defined yet.
utils::globalVariables(c("x", "y", "x.lo", "y.lo", "xintercept",
                         "textX", "textY", "textLabel", "x.hi",
                         "scale_y_continuous", "scale_x_continuous",
                         "geom_ribbon", "geom_vline", "geom_text"));

###########################################################
### Basic theme
###########################################################

### Theme used for the plots
didacticPlotTheme <- function(base_size = 14, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(colour="#000000", size = rel(0.8)),
      axis.ticks        = element_line(colour = "black"),
      legend.text       = element_text(size = rel(0.6)),
      legend.key        = element_rect(colour = "grey80"),
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.key.size   = unit(6, "mm"),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major  = element_line(colour = "grey90", size = 0.2),
      panel.grid.minor  = element_line(colour = "grey98", size = 0.5),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      strip.background  = element_rect(fill = "grey80", colour = "grey50"),
      panel.margin      = unit(c(.5), "cm")
    )
}

###########################################################
### The function itself
###########################################################

didacticPlot <- function(foundValue, statistic, df1, df2 = NULL,
                         granularity = 1000, xLim = NULL, yLab = NULL,
                         lineCol = "red", lineSize=1,
                         surfaceCol = "red", textMarginFactor = 20,
                         sided="two") {
  ###
  ### foundValue       = found value of relevant statistics;
  ### statistic        = one of "r", "t", "f" or "chisq"
  ### df1              = degrees of freedom for r, t and chi^2 test;
  ###                    degrees of freedom for denominator for F-test;
  ### df2              = degrees of freedom for numerator for F-test;
  ### granularity      = steps to use for x-axis;
  ### xLim             = vector; minimum and maximum values on x axis.
  ### yLab             = label on y axis
  ### lineCol          = colour of density line
  ### lineSize         = size of density line
  ### surfaceCol       = colour of coloured surface area
  ### textMarginFactor = used to calculate how close to the vertical line
  ###                    text labels should appear
  ###                     
  
  ### Generate object to return results
  res <- list();
  ### Store parameters
  res$statistic <- statistic;
  res$df1 <- df1;
  res$df2 <- df2;
  res$granularity <- granularity;
  res$foundValue <- foundValue;
  res$foundValue.hi <- abs(foundValue);
  
  if (is.null(yLab)) {
    if (is.null(df2)) {
      yLab <- paste0("Density for ", statistic, " when df=", df1);
    }
    else {
      yLab <- paste0("Density for ", statistic, " when df1=", df1, " and df2=", df2);
    }
  }
  
  ### Set limits for x axis if none were provided
  if (is.null(xLim)) {
    if (statistic == "r") {
      xLim <- c(-1, 1);
    }
    else if (statistic == "t") {
      xLim <- c(-3, 3);
    }
    else if (statistic == "f") {
      xLim <- c(0, 15);
    }
    else if (statistic == "chisq") {
      xLim <- c(0, 5*df1);
    }
  }
  ### Store in result object
  res$xLim <- xLim;
  
  ### Both r and t can be either positive or negative
  if ((statistic == "r") | (statistic == "t")) {
    res$foundValue.lo <- ifelse(foundValue < 0, foundValue, -1 * foundValue);
  }
  else {
    res$foundValue.lo <- xLim[1];
  }
  
  ### Generate data for X axis using xLim and granularity
  res$x <- seq(xLim[1], xLim[2], by=(xLim[2]-xLim[1])/granularity);
  
  ### Generate data for y axis, get probabilities,
  ### set the label for the x axis, and determine
  ### at which density the vertical line that
  ### denotes the found value intersects the
  ### density line (we need that for the plot)
  if (statistic == "r") {
    res$y <- dPearson(res$x, df1);
    res$foundDensity.lo <- pPearson(res$foundValue, df1);
    res$foundDensity.hi <- 1 - pPearson(res$foundValue, df1);
    res$xLab <- "Pearson r";
    res$unit <- "r";
    res$pvalUnit <- "frac(1, 2) ~~ p";
    yIntersection <- dPearson(res$foundValue, df1);
  }
  else if (statistic == "t") {
    res$y <- dt(res$x, df1);
    res$foundDensity.lo <- pt(res$foundValue, df1);
    res$foundDensity.hi <- 1 - pt(res$foundValue, df1);
    res$xLab <- "Student t";
    res$unit <- "t";
    res$pvalUnit <- "frac(1, 2) ~~ p";
    yIntersection <- dt(res$foundValue, df1);
  }
  else if (statistic == "f") {
    res$y <- df(res$x, df1, df2);
    res$foundDensity.lo <- xLim[1];
    res$foundDensity.hi <- 1 - pf(res$foundValue, df1, df2);
    res$xLab <- "F";
    res$unit <- "F";
    res$pvalUnit <- "p";
    yIntersection <- df(res$foundValue, df1, df2);
  }
  else if (statistic == "chisq") {
    res$y <- dchisq(res$x, df1);
    res$foundDensity.lo <- xLim[1];
    res$foundDensity.hi <- 1 - pchisq(res$foundValue, df1);
    res$xLab <- "chi^2";
    res$unit <- "chi^2";
    res$pvalUnit <- "p";
    yIntersection <- dchisq(res$foundValue, df1);
  }
  
  ### Create dataframe for ggplot
  res$plotData <- data.frame(x=res$x, y=res$y,
                             x.lo=ifelse(res$x <= res$foundValue.lo, res$x, res$foundValue.lo),
                             x.hi=ifelse(res$x >= res$foundValue.hi, res$x, res$foundValue.hi));
  ### Generate basic plot
  res$plot.basic <- ggplot(data=res$plotData, aes(x=x, y=y, x.lo=x.lo, x.hi=x.hi)) +
    scale_y_continuous(name=yLab) + 
    scale_x_continuous(name=res$xLab) +
    didacticPlotTheme();
  ### Generate plot with line element
  res$plot.line <- res$plot.basic +
    geom_line(color=lineCol, size=lineSize);
  ### Generate plot with coloured surface;
  ### original example used layer with geom="area",
  ### but with low granularities, this doesn't work out,
  ### so we use geom_ribbon with invisible lines.
  #res$plot.surface <- res$plot.line +
  #  layer(mapping = aes(x=x.lo, y=y),
  #        geom = "area", geom_params=list(fill=surfaceCol,alpha=.3)) +
  #  layer(mapping = aes(x=x.hi, y=y),
  #        geom = "area", geom_params=list(fill=surfaceCol,alpha=.3));
  res$plot.surface <- res$plot.line +
    geom_ribbon(aes(x=x.lo, ymin=0, ymax=y), color=NA, fill=surfaceCol, alpha=.3) +
    geom_ribbon(aes(x=x.hi, ymin=0, ymax=y), color=NA, fill=surfaceCol, alpha=.3);
  
  ### Generate plot with vertical line and labels
  
  ### First compute location for text labels and generate label texts
  
  ### Margin to prevent overlap with vertical line
  xMargin <- (xLim[2] - xLim[1]) / textMarginFactor;
  yMargin <- (max(res$y) - min(res$y)) / textMarginFactor;
  
  ### Generate dataset for labels in ggplot
  foundValueLabel <- data.frame(textX = res$foundValue + xMargin,
                                textY = yIntersection + ((max(res$y) - yIntersection) / 2),
                                textLabel = paste0(res$unit, "==", res$foundValue));
  if (res$foundValue > 0) {
    surfaceLabel <- data.frame(textX = res$foundValue + xMargin,
                               textY = yIntersection,
                               textLabel = paste0(res$pvalUnit, "==", round(res$foundDensity.hi, 3)));
  } else {
    surfaceLabel <- data.frame(textX = res$foundValue + xMargin,
                               textY = yIntersection,
                               textLabel = paste0(res$pvalUnit, "==", round(res$foundDensity.lo, 3)));
  }
  
  ### Dataset with X position
  verticalLineData <- data.frame(xintercept = res$foundValue);
  
  ### Add to plot
  res$plot.complete <- res$plot.surface +
    geom_vline(data=verticalLineData,
               aes(xintercept=xintercept), color=surfaceCol, linetype="dashed") +
    geom_text(mapping=aes(x = textX, y = textY, x.lo=textX, x.hi=textX, label = textLabel),
              data=foundValueLabel,
              color=surfaceCol, hjust=0, parse=TRUE) +
    geom_text(mapping=aes(x = textX, y = textY, x.lo=textX, x.hi=textX, label = textLabel),
              data=surfaceLabel,
              color=surfaceCol, hjust=0, parse=TRUE);
  
  ### Return result object
  class(res) <- c('didacticPlot');
  return(res);
  
}

print.didacticPlot <- function(x, ...) {
  print(x$plot.complete, ...);
}
