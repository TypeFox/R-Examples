### This function generates a pdf file with a report
### describing the variables.
scaleDiagnosisToPDF <- function(scaleDiagnosisObject,
                            docTitle = "Scale diagnosis", docAuthor = "Author",
                            pdfLatexPath, rnwPath=getwd(),
                            filename = "scaleDiagnosis",
                            digits=2,
                            rMatrixColsLandscape = 6,
                            pboxWidthMultiplier = 1,
                            scatterPlotBaseSize = 4,
                            maxScatterPlotSize = NULL,
                            pageMargins=15,
                            pval=TRUE) {

  if (!('scaleDiagnosis' %in% class(scaleDiagnosisObject))) {
    stop("Argument 'scaleDiagnosisObject' must have class 'scaleDiagnosis'!");
  }
  
  if (!hasLaTeX(pdfLatexPath)) {
    stop('In path "', pdfLatexPath, '", the file pdflatex.exe (Windows) or ',
         'pdflatex (MacOS or Ubuntu (Linux)) does not exist! Please ',
         'locate the file and provide its path (without the last ',
         'slash). See ?rnwString for more information.');
  }
  
  res <- list();
  res$scaleDiagnosis <- scaleDiagnosisObject;
  
  if (is.null(maxScatterPlotSize)) {
    maxScatterPlotSize <- 21 - ((pageMargins*2)/10);
  }
  
  res$rnw <- rnwString.initiate(docTitle, docAuthor,
                                docClassArgs='a4paper,portrait,10pt',
                                newPage=FALSE, pageMargins=pageMargins);
  
  ### Generate correlation table
  res$rMatrix <- rMatrix(res$scaleDiagnosis$dat,
                         x=names(res$scaleDiagnosis$dat));
      
  ### Include correlation table;
  ### whether to print on a portrait page or
  ### on a landscape page depends on number of
  ### columns and rMatrixColsLandscape
  if (length(res$rMatrix$variables.cols) < rMatrixColsLandscape) {
    res$rnw <-
      paste0(res$rnw,
             '\n\\begin{minipage}{\\textwidth}\n\\maxsizebox{\\textwidth}{\\textheight}{\n');
  }
  else {
    res$rnw <-
      paste0(res$rnw,
             '\\begin{landscape}\n\\maxsizebox{', 297 - 2*pageMargins, 'mm}{', 210 - 2*pageMargins, 'mm}{\n');
  }
  
  res$rnw <-
    paste0(res$rnw, '<< echo=FALSE, results="asis" >>=\n',
           'print(res$rMatrix, digits=digits, output="LaTeX", pboxWidthMultiplier=pboxWidthMultiplier, pval=pval);\n',
           '@\n');

  if (length(res$rMatrix$variables.cols) < rMatrixColsLandscape) {
    res$rnw <- paste0(res$rnw, '}\n\\end{minipage}\n');
  }
  else {
    res$rnw <- paste0(res$rnw, '}\n\\end{landscape}\n');
  }
      
  ### The size of each panel in the scattermatrix depends
  ### on the number of items - therefore, we need to adjust
  ### the plot sizes to the number of items. This is mainly
  ### necessary because in ggpairs.print (which you can
  ### view with "getAnywhere('print.ggpairs');"), the
  ### fontsize is fixed at 15.
  ### knitr wants unit for outputsize, no unit for figure draw
  ### size (but this must be specified in inches).
  if (res$scaleDiagnosis$scaleReliability$input$n.items * scatterPlotBaseSize > maxScatterPlotSize) {
    figSizeInOutput <- maxScatterPlotSize;
  }
  else {
    figSizeInOutput <- res$scaleDiagnosis$scaleReliability$input$n.items * scatterPlotBaseSize;
  }
  ### For two items on a page (i.e. plots of roughly 9x9 cm),
  ### the labels of the plots have roughly the right size,
  ### so we multiply 9 cm with the number of items.
  figSizeToDraw <- (9 / 2.54) * res$scaleDiagnosis$scaleReliability$input$n.items;
  ### If figSizeToDraw is smaller than output size, set to output size
  if (figSizeToDraw < (figSizeInOutput / 2.54)) {
    figSizeToDraw <- figSizeInOutput / 2.54;
  }
  ### Add unit to size in output
  figSizeInOutput <- paste0(figSizeInOutput, "cm");
  
  res$rnw <-
    paste0(res$rnw,
           '\\begin{minipage}{180mm}\n',
           '<< echo=FALSE, warning=FALSE, dev="pdf", fig.width=', figSizeToDraw, ', fig.height=', figSizeToDraw, ', out.width="', figSizeInOutput, '", out.height="', figSizeInOutput, '" >>=\n',
           'print(res$scaleDiagnosis$ggpairs.combined);\n',
           '@\n',
           '<< echo=FALSE, results="asis" >>=\n',
           '@\n',
           '\\end{minipage}%\n');
  
  ### Finalize rnwString and generate the PDF
  res$rnw <- rnwString.terminate(res$rnw);
  rnwString.generate(res$rnw, rnwPath, fileName=filename, pdfLatexPath);
  
  ### Store result for later inspection
  class(res) <- c('scaleDiagnosisPDF');
  return(res);
  
}