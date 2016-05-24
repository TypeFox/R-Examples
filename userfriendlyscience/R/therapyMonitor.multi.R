therapyMonitor.multi <- function(dat = NULL,
                                 variableColumn = NULL,
                                 conditionColumn = NULL,
                                 conditionMoment = NULL,
                                 minLevels = 5,
                                 outputFiles = FALSE,
                                 outputFilePath = getwd(),
                                 outputFormats = c('svg', 'png'),
                                 silent=FALSE,
                                 ...) {
  
  ### Create object for results and store input
  res <- list(input = as.list(environment()),
              intermediate = list(notice=""),
              output = list());

  ### If no dataframe was specified, load it from an SPSS file
  if (is.null(dat)) {
    dat <- getData(errorMessage=paste0("No dataframe specified, and no valid datafile selected in ",
                                       "the dialog I then showed to allow selection of a dataset.",
                                       "Original error:\n\n[defaultErrorMessage]"),
                   use.value.labels=FALSE);
    res$input$dat.name <- paste0("SPSS file imported from ", attr(dat, "filename"));
  }
  else {
    if (!is.data.frame(dat)) {
      stop("Argument 'dat' must be a dataframe or NULL! Class of ",
           "provided argument: ", class(dat));
    }
    res$input$dat.name <- deparse(substitute(dat));
  }

  ### Store number of levels for each variable
  res$intermediate$nlevels <- unlist(lapply(dat, function(x) {
    return(length(unique(x)));
  }));
  
  ### Check whether we'll need to extract variables ourselves
  if ((is.null(conditionColumn) && is.null(conditionMoment)) ||
        is.null(variableColumn)) {
    ### Remove variables with only one valid level
    if (sum(res$intermediate$nlevels < 2) > 0) {
      dat <- dat[res$intermediate$nlevels > 1];    
      
      ### Take these variables out of the list of all variable levels
      res$intermediate$nlevels.valid <-
        res$intermediate$nlevels[res$intermediate$nlevels > 1];
      
      res$intermediate$notice <-
        paste0(res$intermediate$notice, "\n",
               "The following variable(s) have been removed because they had ",
               "less then two levels: ",
               vecTxt(names(res$intermediate$nlevels)[res$intermediate$nlevels < 2]),
               ".");
    }
  }
  
  ### Verify validity of specified condition variable
  if (!is.null(conditionMoment)) {

    ### Create a condition Column to do some more checks later on
    if (is.null(conditionColumn)) {
      conditionColumn <- "condition";
    }
    dat[, conditionColumn] <- c(rep('A', conditionMoment - 1),
                                rep('B', nrow(dat) - conditionMoment + 1));
    
  } else if (!(is.null(conditionColumn))) {
    if (!(conditionColumn %in% names(dat))) {
      stop("The variable name specified in the 'conditionColumn' argument ",
           "(", conditionColumn,
           ") does not exist among the variables in the dataframe!");
    }
  } else {
    
    ### Select variables with fewest levels
    res$intermediate$varsWithMinLevels <-
      res$intermediate$nlevels.valid[res$intermediate$nlevels.valid ==
                                       min(res$intermediate$nlevels.valid)];
    
    conditionColumn <- names(res$intermediate$varsWithMinLevels[1]);
    
    res$intermediate$notice <-
      paste0(res$intermediate$notice, "\n",
             "No conditionColumn was specified; selected the first ",
             "variable with the minimum number of levels (",
             res$intermediate$varsWithMinLevels[1], "), which is ",
             conditionColumn, ".");
    
  }
  res$intermediate$conditionColumn <- conditionColumn;
  
  ### Verify validity of specified measurements variable
  if (!(is.null(variableColumn))) {
    if (any(!(variableColumn %in% names(dat)))) {
      stop("The variable name specified in the 'variableColumn' argument ",
           "(", variableColumn,
           ") does not exist among the variables in the dataframe!");
    }
  } else {
  
    ### Take all variables with at least minLevels levels
    variableColumn <- names(dat[res$intermediate$nlevels.valid > (minLevels - 1)]);
    
    ### Remove conditionColumn
    variableColumn <- variableColumn[variableColumn != conditionColumn];
    
    res$intermediate$notice <-
      paste0(res$intermediate$notice, "\n",
             "No variableColumn was specified; selected all variables ",
             "with at least as many levels as specified in the ",
             "minLevels argument (",
             minLevels, "), which are ", vecTxt(variableColumn),".");
    
  }
  
  ### Remove non-numeric and non-factor variables
  variableColumn <-
    variableColumn[sapply(dat[, variableColumn], class) %in% c('numeric', 'factor')];
  
  ### Convert all factors to numeric variables
  res$intermediate$dat <- massConvertToNumeric(dat);
  
  ### Replace the conditionColumn variable
  res$intermediate$dat[, conditionColumn] <- dat[, conditionColumn];
  
  res$intermediate$therapyMonitors <- list();
  
  res$output$results <- "";
  
  for (currentVariable in variableColumn) {
    
    if (!silent) {
      cat0("\nProcessing variable ", currentVariable, " . . .");
    }
    
    ### Remove missing cases
    tmpDat <- na.omit(res$intermediate$dat[, c(conditionColumn, currentVariable)]);
    
    ### Remove unused levels from the conditionColumn variable
    tmpDat[, conditionColumn] <- factor(tmpDat[, conditionColumn]);
    
    ### Check whether enough valid measurements remain
    if (nlevels(tmpDat[, conditionColumn]) < 2) {
      res$intermediate$therapyMonitors[[currentVariable]] <-
        paste0("After removing the missing values, no observations remained ",
               "for one of the conditions - not processing this variable!");
        if (!silent) {
        cat0(res$intermediate$therapyMonitors[[currentVariable]]);
      }
    } else {
      res$intermediate$therapyMonitors[[currentVariable]] <-
        tryCatch(therapyMonitor(res$intermediate$dat,
                                conditionColumn = conditionColumn,
                                conditionMoment = conditionMoment,
                                variableColumn = currentVariable,
                                outputFile = ifelseObj(outputFiles,
                                                       file.path(outputFilePath, currentVariable),
                                                       NULL),
                                ...),
                 error = function(e) { return(e); });
      res$output$results <- paste0(res$output$results, "\n\n",
                                   "### ", currentVariable, "\n\n",
                                   paste0(capture.output(print(res$intermediate$therapyMonitors[[currentVariable]],
                                                        printPlot=FALSE)), collapse="\n"));
      if (!silent) {
        cat0(" Done!");
      }
    }
  }
  
  if (outputFiles) {
    writeLines(res$output$results, file.path(outputFilePath, "therapyMonitor results.txt"));
    if (!silent) {
      cat0("\n\nWritten files to ", outputFilePath);
    }
  }
  
  class(res) <- 'therapyMonitor.multi';
  return(res);
  
}

print.therapyMonitor.multi <- function(x, ...) {
  
  if (!is.null(x$intermediate$notice)) {
    cat0("\n\n", x$intermediate$notice, "\n");
  }
  
  cat(x$output$results);
  
}