###########################################################
###########################################################
###
### Function to generate a PDF with four panels per page,
### showing some basic item characteristics.
###
### File created by Gjalt-Jorn Peters. Questions? You can
### contact me through http://behaviorchange.eu.
###
###########################################################
###########################################################

### This function generates a pdf file with a report
### describing the variables.
itemInspection <- function(dat, items,
                           docTitle = "Scale inspection", docAuthor = "Author",
                           pdfLaTexPath, rnwPath, filename="itemInspection",
                           convertFactors = TRUE, digits=4) {

  if (!is.list(items)) {
    items <- list('none' <- items);
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
  
  res$describe <- list();
  res$plot <- list();
  res$rnwBit <- list();
  res$normality.sampleDist <- list();
  res$normality.samplingDist <- list();
  
  res$rnw <- rnwString.initiate(docTitle, docAuthor,
                                docClassArgs='a4paper,landscape,11pt',
                                newPage=FALSE);
  
  ### Process items per scale
  for (currentScale in names(items)) {
    ### Create lists for each scale to store results
    res$describe[[currentScale]] <- list();
    res$plot[[currentScale]] <- list();
    res$rnwBit[[currentScale]] <- list();
    res$normality.sampleDist[[currentScale]] <- list();
    res$normality.samplingDist[[currentScale]] <- list();
    ### Process items one by one
    for (currentItem in items[[currentScale]]) {
      ### Check whether this is a numeric vector
      if (!is.numeric(dat[, currentItem])) {
        ### Check whether it's a character vector (if so, abort)
        if (is.character(dat[, currentItem])) {
          stop("One of the supplied variables is a character vector ",
               "(i.e. contains text strings instead of numbers)!");
        }
        else if (is.factor(dat[, currentItem]) && (!convertFactors)) {
          stop("One of the supplied variables is a factor (instead ",
               "of a numeric vector), and convertFactors is FALSE.");
        }
        else if (is.factor(dat[, currentItem]) && (convertFactors)) {
          dat[, currentItem] <- as.numeric(dat[, currentItem]);
        }
        else {
          stop("One of the items is not numeric, not a character vector, ",
               "and not a factor. Normally, items should be numeric vectors!");
        }
      }
      ### Extract univariate descriptives to show
      res$describe[[currentScale]][[currentItem]] <-
        describe(dat[, currentItem])[, c('n', 'mean', 'sd', 'median', 'min', 'max', 'skew', 'kurtosis')];
      ### Generate and store plots for assessment of normality
      res$plot[[currentScale]][[currentItem]] <-
        normalityAssessment(dat[, currentItem],
                            xLabel.sampleDist = 'Sample Distribution',
                            yLabel.sampleDist = 'Density',
                            xLabel.samplingDist = 'Sampling Distribution',
                            yLabel.samplingDist = 'Density',
                            samplingDistLineSize = .5,
                            normalLineSize = .2);
      ### Create dataframe with normality statistics for the
      ### sample distribution
      res$normality.sampleDist[[currentScale]][[currentItem]] <-
        data.frame(sw = c(round(res$plot[[currentScale]][[currentItem]]$sw.sampleDist$statistic, 4),
                          round(res$plot[[currentScale]][[currentItem]]$sw.sampleDist$p.value, 4)),
                   ad = c(round(res$plot[[currentScale]][[currentItem]]$ad.sampleDist@test$statistic, 4),
                          round(res$plot[[currentScale]][[currentItem]]$ad.sampleDist@test$p.value, 4)),
                   ks = c(round(res$plot[[currentScale]][[currentItem]]$ks.sampleDist$statistic, 4),
                          round(res$plot[[currentScale]][[currentItem]]$ks.sampleDist$p.value, 4)));
      row.names(res$normality.sampleDist[[currentScale]][[currentItem]]) <- c('value', 'p-val');
      ### Create dataframe with normality statistics for the
      ### sampling distribution
      res$normality.samplingDist[[currentScale]][[currentItem]] <-
        data.frame(sw = c(round(res$plot[[currentScale]][[currentItem]]$sw.samplingDist$statistic, 4),
                          round(res$plot[[currentScale]][[currentItem]]$sw.samplingDist$p.value, 4)),
                   ad = c(round(res$plot[[currentScale]][[currentItem]]$ad.samplingDist@test$statistic, 4),
                          round(res$plot[[currentScale]][[currentItem]]$ad.samplingDist@test$p.value, 4)),
                   ks = c(round(res$plot[[currentScale]][[currentItem]]$ks.samplingDist$statistic, 4),
                          round(res$plot[[currentScale]][[currentItem]]$ks.samplingDist$p.value, 4)));
      row.names(res$normality.samplingDist[[currentScale]][[currentItem]]) <- c('value', 'p-val');
      ### Generate the minipage for this item, consisting of:
      ###  - name of measure (scale)
      ###  - name of measurement (item)
      ###  - table with univariate descriptives
      ###  - minipage with:
      ###     - plotted sample distribution
      ###     - normality statistics for sample distribution
      ###  - minipage with:
      ###     - plotted sampling distribution
      ###     - normality statistics for sampling distribution
      res$rnwBit[[currentScale]][[currentItem]] <-
        paste0('\\begin{minipage}[t][90mm][t]{133.5mm}\n',
               '\\subsection{ITEM: ',
               sanitizeLatexString(currentItem),
               '}\n',
               '<< echo=FALSE, results="asis" >>=\n',
               '  print(xtable(res$describe[["',
               currentScale, '"]][["', currentItem,
               '"]], digits=c(0, 0, rep(digits, 7))), tabular.environment="tabular",
               print.rownames=FALSE, floating=FALSE);\n',
               '@\n',
               '\\begin{minipage}[t][50mm][t]{60mm}\n',
               '<< echo=FALSE, warning=FALSE, dev="pdf", fig.width=6/2.54, fig.height=5/2.54 >>=\n',
               'res$plot[["', currentScale, '"]][["', currentItem, '"]]$plot.sampleDist;\n',
               '@\n',
               '<< echo=FALSE, results="asis" >>=\n',
               '  print(xtable(res$normality.sampleDist[["',
               currentScale, '"]][["', currentItem,
               '"]], digits=digits), tabular.environment="tabular",
               floating=FALSE);\n',
               '@\n',
               '\\end{minipage}%\n',
               '\\begin{minipage}[t][50mm][t]{60mm}\n',
               '<< echo=FALSE, warning=FALSE, dev="pdf", fig.width=6/2.54, fig.height=5/2.54 >>=\n',
               'res$plot[["', currentScale, '"]][["', currentItem, '"]]$plot.samplingDist;\n',
               '@\n',
               '<< echo=FALSE, results="asis" >>=\n',
               '  print(xtable(res$normality.samplingDist[["',
               currentScale, '"]][["', currentItem,
               '"]], digits=digits), tabular.environment="tabular",
               floating=FALSE);\n',
               '@\n',
               '\\end{minipage}%\n\\\\\n',
               '\\end{minipage}');
      
    }
  }
  
  ### Combine all minipages into one character vector.
  ### Every two minipages, go to the next line;
  ### every four minipages, go to the next page
  panelCounter <- 0;
  for (currentScale in names(items)) {
    res$rnwPanels <- paste0(res$rnwPanels,
                            '\n\\newpage\n\\section{SCALE: ',
                            sanitizeLatexString(currentScale),
                            '}\n\\newpage\n');
    panelCounter <- 0;
    for (currentItem in items[[currentScale]]) {
      panelCounter <- panelCounter + 1;
      res$rnwPanels <- paste0(res$rnwPanels,
                        res$rnwBit[[currentScale]][[currentItem]]);
      if (panelCounter %% 2 == 0) {
        res$rnwPanels <- paste0(res$rnwPanels,
                          '\n');
      }
      if (panelCounter %% 4 == 0) {
        res$rnwPanels <- paste0(res$rnwPanels,
                          '\n\\newpage');
      }
    }
  }
  
  overviewCountString <-
    paste0('CONTENTS: ', panelCounter, ' panels (measurements/items) in ',
           length(names(items)), ' measures (scales).\n\\tableofcontents\n\\newpage');
  
  ### Combine all pieces
  res$rnw <- paste0(res$rnw, overviewCountString,
                    "\n\\newpage\n", res$rnwPanels);
  
  ### Finalize rnwString and generate the PDF
  res$rnw <- rnwString.terminate(res$rnw);
  rnwString.generate(res$rnw, rnwPath, fileName=filename, pdfLaTexPath);
  
  ### Store result for later inspection
  class(res) <- c('itemInspection');
  return(res);
  
}