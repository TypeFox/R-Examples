###########################################################
###########################################################
###
### Function to generate a PDF with some diagnostics to
### assess how the elements (usually items) in a scale
### relate to each other.
###
### File created by Gjalt-Jorn Peters. Questions? You can
### contact me through http://behaviorchange.eu.
###
###########################################################
###########################################################

### This function generates a confidence level for a standard deviation
### http://www.graphpad.com/guides/prism/6/statistics/index.htm?stat_confidence_interval_of_a_stand.htm
### https://www.wolframalpha.com/input/?i=confidence+interval+for+a+standard+deviation&lk=3
sdConfInt <- function(vector=NULL, sd=NULL, n=NULL, conf.level=.95) {
  if (is.null(sd) & is.null(n)) {
    if (is.null(vector)) {
      stop("Please specify either vector, or sd and n!");
    }
    sd <- sd(vector);
    n <- length(vector);
  }
  res <- list();
  res$input <- list(vector=vector, sd=sd, n=n, conf.level=conf.level);
  res$intermediate <- list(alpha = 1-conf.level);
  res$intermediate$chisq.bound.lo <- qchisq(((1-res$intermediate$alpha)/2), n-1);
  res$intermediate$chisq.bound.hi <- qchisq(res$intermediate$alpha/2, n-1);
  ci.lo <- sqrt(((n-1)*sd^2)/res$intermediate$chisq.bound.lo);
  ci.hi <- sqrt(((n-1)*sd^2)/res$intermediate$chisq.bound.hi);
  res$output <- list(ci = c(ci.lo, ci.hi));
  class(res) <- 'sdConfInt';
  return(res);
}

print.sdConfInt <- function(x, digits=2, ...) {
  print(x$output$ci, digits=digits, ...);
}

### This function generates a confidence level for a single mean
meanConfInt <- function(vector=NULL, mean=NULL, sd=NULL, n=NULL, se=NULL, conf.level=.95) {
  if (is.null(mean) & is.null(sd) & is.null(n) & is.null(se)) {
    if (is.null(vector)) {
      stop("Please specify either vector, or a mean and then also either sd and n or se!");
    }
    mean <- mean(vector);
    sd <- sd(vector);
    n <- length(vector);
    se <- sd/sqrt(n);
  }
  else if (!is.null(mean) & !is.null(sd) & !is.null(n)) {
    se <- sd/sqrt(n);
  }
  else if (is.null(mean) | is.null(se)) {
    stop("Please specify either vector, or a mean and then also either sd and n or se!");
  }
    
  res <- list();
  res$input <- list(vector=vector, mean=mean, sd=sd, n=n, se=se, conf.level=conf.level);
  res$intermediate <- list(alpha = 1-conf.level);
  res$intermediate$t.bound.lo <- qt(res$intermediate$alpha/2, df=n-1);
  res$intermediate$t.bound.hi <- qt(1-res$intermediate$alpha/2, df=n-1);
  ci.lo <- mean + res$intermediate$t.bound.lo * se;
  ci.hi <- mean + res$intermediate$t.bound.hi * se;
  res$output <- list(ci = c(ci.lo, ci.hi));
  class(res) <- 'meanConfInt';
  return(res);
}

print.meanConfInt <- function(x, digits=2, ...) {
  print(x$output$ci, digits=digits, ...);
}

### This function actually makes the scales
makeScales <- function(dat, scales) {
  for (currentScale in 1:length(scales)) {
    if (length(unlist(scales[currentScale])) > 1) {
      dat[[names(scales[currentScale])]] <-
        rowMeans(dat[, unlist(scales[currentScale])], na.rm=TRUE);
    }
    else if (length(unlist(scales[currentScale])) == 1) {
      dat[[names(scales[currentScale])]] <- dat[[unlist(scales[currentScale])]];
    }
  }
  return(dat);
}

### This function generates a pdf file with a report
### describing the variables.
scaleInspection <- function(dat, items=NULL,
                            docTitle = "Scale inspection", docAuthor = "Author",
                            pdfLaTexPath, rnwPath=getwd(),
                            filename = "scaleInspection", convertFactors=TRUE,
                            scaleReliability.ci=FALSE, conf.level=.95, digits=2,
                            rMatrixColsLandscape = 6,
                            pboxWidthMultiplier = 1,
                            scatterPlotBaseSize = 4,
                            pageMargins=15, show = FALSE,
                            pval=TRUE) {

  if (is.null(items)) {
    items <- names(dat);
  }
  
  if (!is.list(items)) {
    items <- list('none' = items);
  }
  if (FALSE %in% sapply(items, is.vector)) {
    ### If this happens, that means there 's non-vector
    ### element in the list of items, which is a problem
    stop("There is a non-vector element in the list with items:\n",
         print(sapply(items, is.vector)));
  }
  if (FALSE %in% sapply(items, is.character)) {
    ### If this happens, that means there 's non-character 
    ### vector element in the list of items, which is a problem
    stop("There is a non-character vector element in the list with items.\n",
         "Here follows a list of booleans indicating which of the list elements ",
         "are character vectors (TRUE) and which aren't (FALSE):\n",
         paste0(names(lapply(items, is.character)), ': ', lapply(items, is.character), '\n'));
  }
  
  if (!hasLaTeX(pdfLaTexPath)) {
    stop('In path "', pdfLaTexPath, '", the file pdflatex.exe (Windows) or ',
         'pdflatex (MacOS or Ubuntu (Linux)) does not exist! Please ',
         'locate the file and provide its path (without the last ',
         'slash). See ?rnwString for more information.');
  }

  res <- list();
  res$items <- items;
  res$scales <- list();
  res$describe <- list();
  res$scaleDiagnosis <- list();
  res$scaleDiagnosis.errors <- list();
  res$normality.plots <- list();
  res$rnwBit <- list();
  res$normality.sampleDist <- list();
  res$normality.samplingDist <- list();
  res$rMatrix <- list();
  res$show <- show;
  
  res$rnw <- rnwString.initiate(docTitle, docAuthor,
                                docClassArgs='a4paper,portrait,10pt',
                                newPage=FALSE, pageMargins=pageMargins);
  
  ### Process each scale separately
  for (currentScale in names(items)) {
    ### Skip scales with only one item
    if (length(items[[currentScale]]) < 2) {
      res$rnwBit[[currentScale]] <-
        paste0('\\newpage\n',
               '\\section{SCALE: ',
               sanitizeLatexString(currentScale),
               '}\n',
               sanitizeLatexString(paste(items[[currentScale]], collapse=", ")),
               '\n\n\\vspace{1ex}\n',
               'Schale has less than 2 items - skipping this one.\n',
               '\\newline\n');
    }
    else {

      ### Check whether there is a non-numeric vector
      if (FALSE %in% lapply(dat[, items[[currentScale]]], 'is.numeric')) {
        ### Check whether there is a character vector (if so, abort)
        if (TRUE %in% lapply(dat[, items[[currentScale]]], 'is.character')) {
          stop("One of the supplied variables is a character vector ",
               "(i.e. contains text strings instead of numbers)!");
        }
        else if ((TRUE %in% lapply(dat[, items[[currentScale]]], 'is.factor')) && (!convertFactors)) {
          stop("One of the supplied items is a factor (instead ",
               "of a numeric vector), and convertFactors is FALSE.");
        }
        else if ((TRUE %in% lapply(dat[, items[[currentScale]]], 'is.factor')) && convertFactors) {
          dat[, items[[currentScale]]] <- data.frame(lapply(dat[, items[[currentScale]]], 'as.numeric'));
        }
        else {
          stop("One of the items is not numeric, not a character vector, ",
               "and not a factor. Normally, items should be numeric vectors!");
        }
      }
    
      res$scales[[currentScale]] <- rowMeans(dat[, items[[currentScale]]], na.rm=TRUE);
      ### Extract univariate descriptives to show
      res$describe[[currentScale]] <-
        describe(res$scales[[currentScale]])[, c('n', 'mean', 'sd', 'median', 'min', 'max', 'skew', 'kurtosis')];
      ### Generate and store plots for assessment of normality
      res$normality.plots[[currentScale]] <-
        normalityAssessment(res$scales[[currentScale]],
                            xLabel.sampleDist = 'Sample Distribution',
                            yLabel.sampleDist = 'Density',
                            xLabel.samplingDist = 'Sampling Distribution',
                            yLabel.samplingDist = 'Density',
                            samplingDistLineSize = .5,
                            normalLineSize = .2);
      ### Create dataframe with normality statistics for the
      ### sample distribution
      res$normality.sampleDist[[currentScale]] <-
        data.frame(sw = c(round(res$normality.plots[[currentScale]]$sw.sampleDist$statistic, 4),
                          round(res$normality.plots[[currentScale]]$sw.sampleDist$p.value, 4)),
                   ad = c(round(res$normality.plots[[currentScale]]$ad.sampleDist@test$statistic, 4),
                          round(res$normality.plots[[currentScale]]$ad.sampleDist@test$p.value, 4)),
                   ks = c(round(res$normality.plots[[currentScale]]$ks.sampleDist$statistic, 4),
                          round(res$normality.plots[[currentScale]]$ks.sampleDist$p.value, 4)));
      row.names(res$normality.sampleDist[[currentScale]]) <- c('value', 'p-val');
      ### Create dataframe with normality statistics for the
      ### sampling distribution
      res$normality.samplingDist[[currentScale]] <-
        data.frame(sw = c(round(res$normality.plots[[currentScale]]$sw.samplingDist$statistic, 4),
                          round(res$normality.plots[[currentScale]]$sw.samplingDist$p.value, 4)),
                   ad = c(round(res$normality.plots[[currentScale]]$ad.samplingDist@test$statistic, 4),
                          round(res$normality.plots[[currentScale]]$ad.samplingDist@test$p.value, 4)),
                   ks = c(round(res$normality.plots[[currentScale]]$ks.samplingDist$statistic, 4),
                          round(res$normality.plots[[currentScale]]$ks.samplingDist$p.value, 4)));
      row.names(res$normality.samplingDist[[currentScale]]) <- c('value', 'p-val');
      ### Generate scale diagnosis
      tryCatch(
        res$scaleDiagnosis[[currentScale]] <-
          scaleDiagnosis(dat, as.vector(items[[currentScale]]),
                         scaleReliability.ci=scaleReliability.ci,
                         conf.level=conf.level)
        , error = function(e) {
          res$scaleDiagnosis.errors[[currentScale]] <- e;
        }
      );
      ### Generate correlation table
      res$rMatrix[[currentScale]] <- rMatrix(dat, items[[currentScale]]);
      
      ### Generate the content:
      ###  - name of measure (scale)
      ###  - list of items in scale
      ###  - table with univariate descriptives of scale
      ###  - minipage with:
      ###     - plotted sample distribution
      ###     - normality statistics for sample distribution
      ###  - minipage with:
      ###     - plotted sampling distribution
      ###     - normality statistics for sampling distribution
      ###  - internal consistency coefficients
      ###  - principal component analysis
      ###  - ggpairs plot of scatterplots, correlations, and histograms
      
      res$rnwBit[[currentScale]] <-
        paste0('\\newpage\n',
               '\\section{SCALE: ',
               sanitizeLatexString(currentScale),
               '}\n',
               sanitizeLatexString(paste(items[[currentScale]], collapse=", ")),
               '\n\n\\vspace{1ex}\n',
               '<< echo=FALSE, results="asis" >>=\n',
               '  print(xtable(res$describe[["',
               currentScale,
               '"]], digits=c(0, 0, rep(digits, 7))), tabular.environment="tabular",
               print.rownames=FALSE, floating=FALSE);\n',
               '@\n',
               '\\vspace{1ex}\\begin{minipage}[t]{80mm}\n',
               '<< echo=FALSE, warning=FALSE, dev="pdf", fig.width=8/2.54, fig.height=8/2.54 >>=\n',
               'res$normality.plots[["', currentScale, '"]]$plot.sampleDist;\n',
               '@\n',
               '\\vspace{1ex}\n<< echo=FALSE, results="asis" >>=\n',
               '  print(xtable(res$normality.sampleDist[["',
               currentScale,
               '"]], digits=digits), tabular.environment="tabular",
               floating=FALSE);\n',
               '@\n',
               '\\end{minipage}%\n',
               '\\begin{minipage}[t]{80mm}\n',
               '<< echo=FALSE, warning=FALSE, dev="pdf", fig.width=8/2.54, fig.height=8/2.54 >>=\n',
               'res$normality.plots[["', currentScale, '"]]$plot.samplingDist;\n',
               '@\n',
               '\\vspace{1ex}\n<< echo=FALSE, results="asis" >>=\n',
               'print(xtable(res$normality.samplingDist[["',
               currentScale,
               '"]], digits=digits), tabular.environment="tabular",
               floating=FALSE);\n',
               '@\n',
               '\\end{minipage}%\n\\newline\n');
      
      if (is.null(res$scaleDiagnosis[[currentScale]]$scaleReliability)) {
        stop(paste0("No scaleReliability object present for scale ", currentScale, "!"));
      }
      
      if (res$scaleDiagnosis[[currentScale]]$scaleReliability$input$n.items > 2) {
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '\\vspace{1ex}\n<< echo=FALSE, results="asis" >>=\n',
                 'cat(paste0("\n\nOmega: ", round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$scaleReliability$output$omega, digits), "\n\nGreatest Lower Bound (GLB): ", round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$scaleReliability$output$glb, digits), "\n\\nCronbach\'s alpha: ", round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$scaleReliability$output$cronbach.alpha, digits), "\n\n',
                 'Eigen values: ", paste(round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$eigen$values, digits), collapse=", "), "\n\n',
                 'Number of factors with Eigen value over 1: ", res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$factors, "\n\n"));\n',
                 '@\n');
        ### Show principal component analysis
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '\\vspace{1ex}\n\\begin{minipage}{\\linewidth}\\begin{verbatim}\n',
                 '<< echo=FALSE, results="asis" >>=\n',
                 'print(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$pca$loadings, digits=digits);\n',
                 '@\n',
                 '\\end{verbatim}\\end{minipage}\n');
  
      } else if (res$scaleDiagnosis[[currentScale]]$scaleReliability$input$n.items == 2) {
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '\\vspace{1cm}\n<< echo=FALSE, results="asis" >>=\n',
                 '  cat(paste0("\n\nSpearman Brown coefficient: ", round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$scaleReliability$output$spearman.brown, digits), "\n\nCronbach\'s alpha: ", round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$scaleReliability$output$cronbach.alpha, digits), "\n\nPearson correlation: ", round(res$scaleDiagnosis[["',
                 currentScale,
                 '"]]$scaleReliability$intermediate$cor[1,2], digits), "\n\n"));\n',
                 '@\n\\vspace{1cm}');
      }
      
      ### Include correlation table;
      ### whether to print on a portrait page or
      ### on a landscape page depends on number of
      ### columns and rMatrixColsLandscape
      if (length(res$rMatrix[[currentScale]]$variables.cols) < rMatrixColsLandscape) {
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '\n\\begin{minipage}{\\textwidth}\n\\maxsizebox{\\textwidth}{\\textheight}{\n');
      }
      else {
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '\\begin{landscape}\n\\maxsizebox{', 297 - 2*pageMargins, 'mm}{', 210 - 2*pageMargins, 'mm}{\n');
      }
      res$rnwBit[[currentScale]] <-
        paste0(res$rnwBit[[currentScale]],
               '<< echo=FALSE, results="asis" >>=\n',
               'print(res$rMatrix[["',
               currentScale,
               '"]], digits=digits, output="LaTeX", pboxWidthMultiplier=pboxWidthMultiplier, pval=pval);\n',
               '@\n');
      if (length(res$rMatrix[[currentScale]]$variables.cols) < rMatrixColsLandscape) {
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '}\n\\end{minipage}\n');
      }
      else {
        res$rnwBit[[currentScale]] <-
          paste0(res$rnwBit[[currentScale]],
                 '}\n\\end{landscape}\n');
      }
      
      ### The size of each panel in the scattermatrix depends
      ### on the number of items - therefore, we need to adjust
      ### the plot sizes to the number of items. This is mainly
      ### necessary because in ggpairs.print (which you can
      ### view with "getAnywhere('print.ggpairs');"), the
      ### fontsize is fixed at 15.
      ### knitr wants unit for outputsize, no unit for figure draw
      ### size (but this must be specified in inches).
      if (res$scaleDiagnosis[[currentScale]]$scaleReliability$input$n.items * scatterPlotBaseSize > 18) {
        figSizeInOutput <- 18;
      }
      else {
        figSizeInOutput <- res$scaleDiagnosis[[currentScale]]$scaleReliability$input$n.items * scatterPlotBaseSize;
      }
      ### For two items on a page (i.e. plots of roughly 9x9 cm),
      ### the labels of the plots have roughly the right size,
      ### so we multiply 9 cm with the number of items.
      figSizeToDraw <- (9 / 2.54) * res$scaleDiagnosis[[currentScale]]$scaleReliability$input$n.items;
      ### If figSizeToDraw is smaller than output size, set to output size
      if (figSizeToDraw < (figSizeInOutput / 2.54)) {
        figSizeToDraw <- figSizeInOutput / 2.54;
      }
      ### Add unit to size in output
      figSizeInOutput <- paste0(figSizeInOutput, "cm");
      
      res$rnwBit[[currentScale]] <-
        paste0(res$rnwBit[[currentScale]],
               '\\begin{minipage}{180mm}\n',
               '<< echo=FALSE, warning=FALSE, dev="pdf", fig.width=', figSizeToDraw, ', fig.height=', figSizeToDraw, ', out.width="', figSizeInOutput, '", out.height="', figSizeInOutput, '" >>=\n',
               'print(res$scaleDiagnosis[["', currentScale, '"]]$scatterMatrix);\n',
               '@\n',
               '<< echo=FALSE, results="asis" >>=\n',
               '@\n',
               '\\end{minipage}%\n');
    }
  }

  ### Combine all pages generated for each scale
  for (currentScale in names(items)) {
    res$rnwPanels <- paste(res$rnwPanels, res$rnwBit[[currentScale]]);
  }
  
  overviewCountString <-
    paste0('CONTENTS: ', length(names(items)), ' measures (scales):\n\\tableofcontents');
  
  ### Combine all pieces
  res$rnw <- paste0(res$rnw, overviewCountString,
                    "\n\\newpage\n", res$rnwPanels);
  
  ### Finalize rnwString and generate the PDF
  res$rnw <- rnwString.terminate(res$rnw);
  rnwString.generate(res$rnw, rnwPath, fileName=filename, pdfLaTexPath);
  
  ### Store result for later inspection
  class(res) <- c('scaleInspection');
  return(res);
  
}

print.scaleInspection <- function (x, show=x$show, ...) {
  if (show) {
    print(x, ...);
  }
  invisible();
}
