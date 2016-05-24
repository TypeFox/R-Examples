therapyMonitor <- function(dat = NULL,
                           design="AB",
                           statistic="|A-B|",
                           conditionColumn = NULL,
                           variableColumn = NULL,
                           timeColumn = NULL,
                           conditionMoment = NULL,
                           limit=NULL, lines=NULL,
                           ylab=NULL, xlab=NULL,
                           outputFile = NULL,
                           outputFormats = c('svg', 'png'),
                           plotTitle = "therapyMonitor results",
                           plotWidth=25,
                           plotHeight=15) {
  
#   ### For easy updating if SCRT is updated
#   validDesigns <- toupper(c("AB", "ABA", "ABAB", "CRD",
#                             "RBD", "ATD", "MBD"));
#   validStatistics <- toupper(c("A-B", "B-A", "|A-B|", "PA-PB", "PB-PA",
#                                "|PA-PB|", "AA-BB", "BB-AA", "|AA-BB|"));
#   
#   ### Convert design and statistic parameters to uppercase
#   design <- toupper(design);
#   statistic <- toupper(statistic);
  
  ### Create object for results and store input
  res <- list(input = as.list(environment()),
              intermediate = list(),
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
  
  res$intermediate$rawLines <- nrow(dat);

#   ### Verify validity of design argument
#   if (!(design %in% validDesigns)) {
#     stop("The value of parameter 'design' is invalid!");
#   }
#   
#   ### Verify validity of statistic argument
#   if (!(statistic %in% validStatistics)) {
#     stop("The value of parameter 'statistic' is invalid!");
#   }
  
  ### Verify validity of limit argument
  if (is.null(limit)) {
    if (is.null(conditionColumn)) {
      res$intermediate$limit <- max(2, min(conditionMoment-1,
                                           nrow(dat)-conditionMoment+1));
    } else {
      res$intermediate$limit <- max(2, floor(min(table(dat[, conditionColumn])) / 2));
    }
  } else if (limit < 0) {
    stop("The 'limit' value must be at least 0!");
  } else {
    res$intermediate$limit <- limit;
  }  
  
  if (is.null(lines)) {
    res$intermediate$lines <- lines <- 1:nrow(dat);
  } else {
    res$intermediate$lines <- lines;
  }
  
  dat <- data.frame(dat[lines, ]);
  
  #########################################################################
  ### Check the validity of the arguments
  #########################################################################
  
  ### Verify validity of specified condition variable
  if (is.null(conditionMoment)) {
    if (!(is.null(conditionColumn))) {
      if (!(conditionColumn %in% names(dat))) {
        stop("The variable name specified in the 'conditionColumn' argument ",
             "(", conditionColumn,
             ") does not exist among the variables in the dataframe!");
      }
    } else {
      conditionColumn <- names(dat)[1];
    }
    res$intermediate$conditionColumn <- conditionColumn;
  } else {
    ### Verify validity of conditionMoment and if it's valid, determine where
    ### the conditions change
    if (conditionMoment < res$intermediate$limit) {
      stop("The conditionMoment argument that was specified (", conditionMoment,
           ") is lower than the specified minimum number of measurements ",
           "within one condition (", res$intermediate$limit, ").");
    }
    if (conditionMoment > (nrow(dat) - res$intermediate$limit)) {
      stop("The conditionMoment argument that was specified (", conditionMoment,
           ") is higher than the number of measurements (", nrow(dat),
           ") minus the specified minimum number of measurements ",
           "within one condition (", res$intermediate$limit, "), which means that using this ",
           "conditionMoment would result in too few measurements in the ",
           "second condition (i.e. less than the specified limit).");
    }
    if (is.null(conditionColumn)) {
      conditionColumn <- "condition";
    }
    dat[, conditionColumn] <- c(rep('A', conditionMoment - 1),
                                rep('B', nrow(dat) - conditionMoment + 1));
    if (is.null(res$input$lines)) {
      res$intermediate$notice <-
        paste0("Note: the 'conditionMoment' argument was used to specify the ",
               "shift between conditions, and the 'lines' argument was not ",
               "used to select a subset of measurements from the datafile. ",
               "If the dataset in fact contains measurements from more ",
               "than two consecutive conditions, this will have resulted ",
               "in combining different conditions into one.");
    }
  }
  
  ### Verify validity of specified measurements variable
  if (!(is.null(variableColumn))) {
    if (!(variableColumn %in% names(dat))) {
      stop("The variable name specified in the 'variableColumn' argument ",
           "('", variableColumn,
           "') does not exist among the variables in the dataframe!");
    }
  } else {
    variableColumn <- names(dat)[2];
  }
  res$intermediate$variableColumn <- variableColumn;
    
  #########################################################################
  ### Select data and check whether number of conditions ('phases') is
  ### correct
  #########################################################################
  
  if (is.null(timeColumn)) {
    res$intermediate$usedData <- dat <- na.omit(dat[, c(conditionColumn, variableColumn)]);
  } else {
    res$intermediate$usedData <- dat <- na.omit(dat[, c(conditionColumn, variableColumn, timeColumn)]);
    
    ### If timeColumn is a numeric variable, check which system it used and
    ### convert it into a POSIX datetime format:
    ###  13166060400, tz=CET = 1/1/2000 in SPSS (from 1582-10-14)
    startnumber.spss <- 13166060400;
    origin.spss <- "1582-10-14";
    ###  1262300400, tz=CET = 1/1/2000 in SAS (from 1960-01-01)
    startnumber.sas <- 1262300400;
    origin.sas <- "1960-01-01";
    ###  1262300400000, tz=CET = 1/1/2000 in Stata (in milliseconds, from 1960-01-01
    startnumber.stata <- 1262300400000;
    origin.stata <- "1960-01-01";
    
    if (is.numeric(res$intermediate$usedData[, timeColumn])) {

      ### Set the regr and plot timeColumn variables
      timeColumn.plot <- paste0(timeColumn, "_plot");
      timeColumn.regr <- paste0(timeColumn, "_regr");
      
      if (max(res$intermediate$usedData[, timeColumn], na.rm=TRUE) > startnumber.stata) {
        ### This is a datetime in Stata format
        res$intermediate$usedData[, timeColumn.plot] <-
          as.POSIXct(res$intermediate$usedData[, timeColumn] / 1000,
                     origin = origin.stata,
                     tz="CET");
        res$intermediate$usedData[, timeColumn.regr] <-
          as.numeric(res$intermediate$usedData[, timeColumn.plot]);
      } else if (max(res$intermediate$usedData[, timeColumn], na.rm=TRUE) > startnumber.spss) {
        res$intermediate$usedData[, timeColumn.plot] <-
          as.POSIXct(res$intermediate$usedData[, timeColumn],
                     origin = origin.spss,
                     tz="CET");        
        res$intermediate$usedData[, timeColumn.regr] <-
          as.numeric(res$intermediate$usedData[, timeColumn.plot]);
      } else if (max(res$intermediate$usedData[, timeColumn], na.rm=TRUE) > startnumber.sas) {
        res$intermediate$usedData[, timeColumn.plot] <-
          as.POSIXct(res$intermediate$usedData[, timeColumn],
                     origin = origin.stata,
                     tz="CET");
        res$intermediate$usedData[, timeColumn.regr] <-
          as.numeric(res$intermediate$usedData[, timeColumn.plot]);
      }
    } else if (any(class(res$intermediate$usedData[, timeColumn]) %in%
                 c("POSIXct", "POSIXt"))) {
      ### Set the timeColumn variable as the 'regr' timeColumn
      timeColumn.plot <- timeColumn;
      timeColumn.regr <- paste0(timeColumn, "_regr");

      res$intermediate$usedData[, timeColumn.regr] <-
        as.numeric(res$intermediate$usedData[, timeColumn]);
      
    } else if (is.character(res$intermediate$usedData[, timeColumn])) {

      ### Set the regr and plot timeColumn variables
      timeColumn.plot <- paste0(timeColumn, "_plot");
      timeColumn.regr <- paste0(timeColumn, "_regr");
      
      if (any(grepl('\\d\\d?/\\d\\d?/\\d\\d\\d\\d? \\d\\d?:\\d\\d?:\\d\\d?',
                    res$intermediate$usedData[, timeColumn]))) {
        res$intermediate$usedData[, timeColumn.plot] <-
          as.POSIXct(strptime(res$intermediate$usedData[, timeColumn],
                   "%m/%d/%Y %H:%M:%S"));
      } else if (any(grepl('\\d\\d?-\\d\\d?-\\d\\d\\d\\d? \\d\\d?:\\d\\d?:\\d\\d?',
                           res$intermediate$usedData[, timeColumn]))) {
        res$intermediate$usedData[, timeColumn.plot] <-
          as.POSIXct(strptime(res$intermediate$usedData[, timeColumn],
                              "%d-%m-%Y %H:%M:%S"));
      } else {
        stop("I could not identify the format of the variable specified in ",
             "the 'timeColumn' argument, so I'm aborting.");
      }
      
      res$intermediate$usedData[, timeColumn.regr] <-
        as.numeric(res$intermediate$usedData[, timeColumn.plot]);
      
    } else {
      stop("I could not identify the format of the variable specified in ",
           "the 'timeColumn' argument, so I'm aborting.");
    }
    
    res$intermediate$timeColumn.plot <- timeColumn.plot;
    res$intermediate$timeColumn.regr <- timeColumn.regr;
    
    ### Also sort on the time variable
    res$intermediate$usedData <- dat <-
      res$intermediate$usedData[order(res$intermediate$usedData[, timeColumn.plot]), ];
    
    ### And then reassign conditionMoment, because the order of the data
    ### may have been wrong
    if (!is.null(conditionMoment)) {
      res$intermediate$usedData[, conditionColumn] <- dat[, conditionColumn] <-
        c(rep('A', conditionMoment - 1),
          rep('B', nrow(dat) - conditionMoment + 1));
      if (is.null(res$input$lines)) {
        res$intermediate$notice <-
          paste0("Note: the 'conditionMoment' argument was used to specify the ",
                 "shift between conditions, and the 'lines' argument was not ",
                 "used to select a subset of measurements from the datafile. ",
                 "If the dataset in fact contains measurements from more ",
                 "than two consecutive conditions, this will have resulted ",
                 "in combining different conditions into one.");
      }
    }
    
  }
  
  ### Get indices of unique conditions
  res$intermediate$conditionChanges <- (1:length(dat[, conditionColumn]))[!duplicated(dat[, conditionColumn])];
  ### Remove first one (the first condition)
  res$intermediate$conditionChanges <- res$intermediate$conditionChanges[-1];
  
  #########################################################################
  ### Deal with a situation where too many phases are specified
  #########################################################################

  if (length(res$intermediate$conditionChanges) > 1) {
    firstPhaseObservations <- res$intermediate$conditionChanges[1] - 1;
    secondPhaseObservations <- res$intermediate$conditionChanges[2] - res$intermediate$conditionChanges[1];
    if ((firstPhaseObservations >= res$intermediate$limit) &&
        (secondPhaseObservations >= res$intermediate$limit)) {
      conditionChangesPlusEnd <- c(res$intermediate$conditionChanges, length(dat[, conditionColumn]));
      res$intermediate$notice <-
        paste0("Note: More than one change between conditions specified (the datafile ",
               "contains changes between conditions at lines ",
               vecTxt(res$intermediate$conditionChanges),
               "). Because the number of valid observations for the first two conditions ",
               "(", firstPhaseObservations, " and ", secondPhaseObservations,
               ", respectively) are both higher than the limit (",
               res$intermediate$limit, "), I'll only use those, discarding all ",
               "data after line ", res$intermediate$conditionChanges[2] - 1, ". ",
               "If you want to look at the differences between the other conditions, ",
               "use argument 'lines', for example ",
               vecTxt(paste0("'lines=",
                             conditionChangesPlusEnd[1:(length(conditionChangesPlusEnd)-2)],
                             ":",
                             conditionChangesPlusEnd[3:length(conditionChangesPlusEnd)],
                             "'"), lastDelimiter = " or "), ".");
      res$intermediate$lines <- 1:(res$intermediate$conditionChanges[2] - 1);
      res$intermediate$usedData <- res$intermediate$usedData[res$intermediate$lines, ];
    } else {
      stop("More than 2 changes between conditions specified (the datafile ",
           "contains changes between conditions at lines ",
           vecTxt(res$intermediate$conditionChanges),
           ").");
    }
  }
  
  #########################################################################
  ### Prepare phases variable to fit the input the SCVA and SCRT
  ### functions expect
  #########################################################################
  
  ### Remove any spaces in the condition variable values
  res$intermediate$usedData[, conditionColumn] <- trim(res$intermediate$usedData[, conditionColumn]);
  
  ### Make a factor out of the conditions column (also drops unnused levels)
  res$intermediate$usedData[, conditionColumn] <- factor(res$intermediate$usedData[, conditionColumn]);
  
  ### Check whether removing the missing values hasn't removed all
  ### observations for one of the conditions
  if (nlevels(res$intermediate$usedData[, conditionColumn]) < 2) {
    stop("After removing the missing values, no observations remained ",
         "for one of the conditions!");
  }
  
  ### Get names of levels for the plot
  res$intermediate$conditionNames <- levels(res$intermediate$usedData[, conditionColumn]);
  
  ### Give levels the names 'A' and 'B'
  res$intermediate$usedData[, conditionColumn] <- factor(res$intermediate$usedData[, conditionColumn],
                                                         labels = c("A", "B"));
  
  ### Create a vector with conditionChanges based on the first change
  res$intermediate$generatedConditionVariable <-
    factor(c(rep('A', res$intermediate$conditionChanges[1] - 1),
             rep('B', nrow(res$intermediate$usedData) - res$intermediate$conditionChanges[1] + 1)));

  ### Compare them
  res$intermediate$conditionColumnDeviations <- 
    which(res$intermediate$usedData[, conditionColumn] != res$intermediate$generatedConditionVariable);
  
  ### Check for equality
  if (length(res$intermediate$conditionColumnDeviations) > 0) {
    res$intermediate$originalConditionVariable <-
      res$intermediate$usedData[, conditionColumn];
    res$intermediate$usedData[, conditionColumn] <-
      res$intermediate$generatedConditionVariable;
    res$intermediate$notice <-
      paste0(res$intermediate$notice, "\n",
             "Note: The specified conditionColumn variable, '",
             conditionColumn,
             "', indicates a change between conditions at measurement ",
             res$intermediate$conditionChanges[1], " - but measurements ",
             "after that moment are labelled with the condition occurring ",
             "before that measurement ('",
             res$intermediate$usedData[1, variableColumn],
             "')! Replacing the conditionColumn with a sequence of ",
             "conditions I generated, where all measurements before ",
             "measurement ", res$intermediate$conditionChanges[1],
             " are considered condition A, and measurement ",
             res$intermediate$conditionChanges[1], " and all those ",
             "following will be considered condition B.");
  }
  
  ### Store number of observations for A and B
  res$intermediate$observations <- as.vector(table(res$intermediate$usedData[, conditionColumn]));
  
  ### Check whether these are still enough 
  if (sum(res$intermediate$observations < res$intermediate$limit) > 0) {
    msg <- paste0("Less observations remain than the limit (",
                  res$intermediate$limit,
                  ")! The numbers of valid observations are ",
                  vecTxt(res$intermediate$observations),
                  ".");
    if(!is.null(res$intermediate$notice)) {
      msg <- paste0(msg, " Notices from running the function:\n",
                    res$intermediate$notice);
    }
    stop(msg);
  }
  
  ### Store total number of observations we'll use
  res$intermediate$usedLines <- nrow(res$intermediate$usedData);
  
  ### Set y label
  if (is.null(ylab)) {
    res$intermediate$ylab <- variableColumn;
  } else {
    res$intermediate$ylab <- ylab;
  }

  ### Set x label
  if (is.null(xlab)) {
    res$intermediate$xlab <- paste0("Measurement times (A = condition ",
                                    res$intermediate$conditionNames[1],
                                    ", B = condition ",
                                    res$intermediate$conditionNames[2],
                                    ")");
  } else {
    res$intermediate$xlab <- xlab;
  }
  
  ### Generate dataframe for ggplot
  res$intermediate$ggplotData <- res$intermediate$usedData;
  
  ### Create vector representing 'measurement times'
  if (is.null(timeColumn)) {
    res$intermediate$ggplotData$moment <- 1:nrow(res$intermediate$ggplotData);
    res$intermediate$usedData$moment <- 1:nrow(res$intermediate$usedData);
    timeColumn.regr <- res$intermediate$timeColumn.regr <- 'moment';
    timeColumn.plot <- res$intermediate$timeColumn.plot <- 'moment';
  } else {
    ### Convert the time version for the regression analysis to seconds that
    ### elapsed since the first measurement
    firstTimestamp <- min(res$intermediate$usedData[, timeColumn.regr], na.tm=TRUE);
    
    res$intermediate$usedData[, timeColumn.regr] <-
      res$intermediate$usedData[, timeColumn.regr] - firstTimestamp;
    res$intermediate$ggPlotData[, timeColumn.regr] <-
      res$intermediate$ggPlotData[, timeColumn.regr] - firstTimestamp;
    
    ### Divide by 7*24*60^2 to get to the number of weeks elapsed since first
    ### measurement
    res$intermediate$usedData[, timeColumn.regr] <- 
      res$intermediate$usedData[, timeColumn.regr] / (7*24*60^2);
    res$intermediate$ggPlotData[, timeColumn.regr] <-
      res$intermediate$ggPlotDataData[, timeColumn.regr] / (7*24*60^2);
  }
  
  ### Create datasets for each condition
  res$intermediate$separateDat[[1]] <-
    res$intermediate$usedData[1:(res$intermediate$conditionChanges[1] - 1), ];
  res$intermediate$separateDat[[2]] <-
    res$intermediate$usedData[res$intermediate$conditionChanges[1]: nrow(res$intermediate$usedData), ];
  
  ### Compute regression analyses and means for each condition, if there is
  ### sufficient variation
  res$output[[1]] <- list();
  if (ifelse(is.na(var(res$intermediate$separateDat[[1]][, variableColumn])),
             0, var(res$intermediate$separateDat[[1]][, variableColumn])) > 0) {
    res$output[[1]] <- list(regr =
                              regr(formula(paste0(variableColumn, " ~ ", timeColumn.regr)),
                                   dat=res$intermediate$separateDat[[1]]));
  }

  res$output[[2]] <- list();
  if (ifelse(is.na(var(res$intermediate$separateDat[[2]][, variableColumn])),
             0, var(res$intermediate$separateDat[[2]][, variableColumn])) > 0) {
    res$output[[2]] <- list(regr =
                              regr(formula(paste0(variableColumn, " ~ ", timeColumn.regr)),
                                   dat=res$intermediate$separateDat[[2]]));
  }
 
 ### Compute means for each condition
 res$output[[1]]$mean <-
   mean(res$intermediate$separateDat[[1]][, variableColumn]);
 res$output[[1]]$meanConfInt <-
   meanConfInt(res$intermediate$separateDat[[1]][, variableColumn]);
 
 res$output[[2]]$mean <-
   mean(res$intermediate$separateDat[[2]][, variableColumn]);
 res$output[[2]]$meanConfInt <-
   meanConfInt(res$intermediate$separateDat[[2]][, variableColumn]);
 
 ### Establish date (moment) of intervention
 interventionMoment <- max(as.numeric(res$intermediate$separateDat[[1]][, timeColumn.plot])) +
   ((min(as.numeric(res$intermediate$separateDat[[2]][, timeColumn.plot])) -
       max(as.numeric(res$intermediate$separateDat[[1]][, timeColumn.plot]))) / 2);
 
 res$output$meansPlot <- ggplot(res$intermediate$ggplotData,
                                 aes_string(x=timeColumn.plot,
                                            y=variableColumn,
                                            group=conditionColumn,
                                            color=conditionColumn)) +
    geom_point(size=4) + geom_line(size=1) +
    geom_vline(xintercept = interventionMoment, size=1) +
    geom_smooth(method=lm, aes_string(fill=conditionColumn), formula=y ~ 1, size=1, alpha=.25) +
    dlvTheme();
  
  res$output$regressionPlot <- ggplot(res$intermediate$ggplotData,
                                      aes_string(x=timeColumn.plot,
                                                 y=variableColumn,
                                                 group=conditionColumn,
                                                 color=conditionColumn)) +
    geom_point(size=4) +
    geom_vline(xintercept = interventionMoment, size=1) +
    geom_smooth(method=lm, aes_string(fill=conditionColumn), size=1, alpha=.25) + dlvTheme();
  
  ### Get maximum and minimum y limits to synchronize these between plots
  res$intermediate$lowerLimits <- c(
    min(ggplot_build(res$output$meansPlot)$panel$y_scales[[1]]$range$range),
    min(ggplot_build(res$output$regressionPlot)$panel$y_scales[[1]]$range$range));
  res$intermediate$upperLimits <- c(
    max(ggplot_build(res$output$meansPlot)$panel$y_scales[[1]]$range$range),
    max(ggplot_build(res$output$regressionPlot)$panel$y_scales[[1]]$range$range));
  res$intermediate$plotLimits <- c(
    min(res$intermediate$lowerLimits),
    max(res$intermediate$upperLimits));
  res$output$meansPlot <- res$output$meansPlot +
    scale_y_continuous(limits = res$intermediate$plotLimits);
  res$output$regressionPlot <- res$output$regressionPlot +
    scale_y_continuous(limits = res$intermediate$plotLimits);
  
  ### Combine plots into one
  res$output$combinedPlot <-
    arrangeGrob(
      ### First the title over the full width
      textGrob(plotTitle),
      ### Then a grob with the rest of the contents
      arrangeGrob(textGrob(variableColumn,
                           rot=90),
                  res$output$meansPlot,
                  res$output$regressionPlot,
                  ncol=3,
                  widths=c(.1,1,1)),
      ### Finally one with an empty box and the x axis label
      arrangeGrob(grid.rect(gp=gpar(col="white")),
                  textGrob("Measurement moment (time)"),
                  ncol=2, widths=c(.1, 2)),
      ncol=1, heights=c(.1, 1, .1));

  ########################################################################
  ### Compute the p value of the randomization test
  #########################################################################

  res$output$p <- pvalue.systematic(design = res$input$design,
                                    statistic = statistic,
                                    limit=res$intermediate$limit,
                                    data = res$intermediate$usedData);

  #########################################################################
  ### Write file(s) to disk if we need to
  #########################################################################
  
  if (!is.null(res$input$outputFile)) {
    if ('svg' %IN% res$input$outputFormats) {
      tryCatch({svg(filename=paste0(res$input$outputFile, ".svg"),
               width=plotWidth / 2.54,
               height=plotHeight / 2.54);
           grid.draw(res$output$combinedPlot);
      }, finally= dev.off());
    }
    if ('png' %IN% res$input$outputFormats) {
      tryCatch({png(filename=paste0(res$input$outputFile, ".png"),
               width=plotWidth,
               height=plotHeight,
               units="cm", res=300);
           grid.draw(res$output$combinedPlot);
      }, finally= dev.off());
    }
  }
  
  class(res) <- 'therapyMonitor';
  return(res);
  
}

print.therapyMonitor <- function(x, digits=2, printPlot = TRUE, ...) {
  res <- x;
  cat0(  "         Measured variable: ", res$intermediate$variableColumn,
       "\nCondition (phase) variable: ", res$intermediate$conditionColumn,
       "\n  Observations in raw data: ", res$intermediate$rawLines,
       "\n     Selected observations: lines ", min(res$intermediate$lines),
       ":", max(res$intermediate$lines), " (",
       res$intermediate$usedLines, " observations)",
       "\n                    Values: ", res$intermediate$conditionNames[1],
       " (A, ", res$intermediate$observations[1], " observations) and",
       "\n                            ", res$intermediate$conditionNames[2],
       " (B, ", res$intermediate$observations[2], " observations)");
  cat0("\n\n    Testing the difference: mean of ",
       res$intermediate$conditionNames[1], " = ",
       round(res$output[[1]]$mean, digits), " (95% CI [",
       round(res$output[[1]]$meanConfInt$output$ci[1], digits), "; ",
       round(res$output[[1]]$meanConfInt$output$ci[2], digits), "])",
       "\n                            mean of ",
       res$intermediate$conditionNames[2], " = ",
       round(res$output[[2]]$mean, digits), " (95% CI [",
       round(res$output[[2]]$meanConfInt$output$ci[1], digits), "; ",
       round(res$output[[2]]$meanConfInt$output$ci[2], digits), "])",
       "\n                            ",
       formatPvalue(res$output$p, digits=digits+1));
  cat0("\n\n   Regression coefficients: at ",
       res$intermediate$conditionNames[1], ", \u03B2 = ");
  if (ifelse(is.na(var(res$intermediate$separateDat[[1]][, res$intermediate$variableColumn])),
             0, var(res$intermediate$separateDat[[1]][, res$intermediate$variableColumn])) > 0) {
    cat0(round(res$output[[1]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, 'estimate'], digits), " (95% CI [",
         round(res$output[[1]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, '95% CI, lo'], digits), "; ",
         round(res$output[[1]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, '95% CI, hi'], digits), "]), ",
         formatPvalue(res$output[[1]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, 'p'], digits=digits+1));
  } else {
    cat0("0 - that is, there is no variance among the measurements!");
  }
  cat0("\n                            at ",
       res$intermediate$conditionNames[2], ", \u03B2 = ");
  if (ifelse(is.na(var(res$intermediate$separateDat[[2]][, res$intermediate$variableColumn])),
             0, var(res$intermediate$separateDat[[2]][, res$intermediate$variableColumn])) > 0) {
    cat0(round(res$output[[2]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, 'estimate'], digits), " (95% CI [",
         round(res$output[[2]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, '95% CI, lo'], digits), "; ",
         round(res$output[[2]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, '95% CI, hi'], digits), "]), ",
         formatPvalue(res$output[[2]]$regr$output$coef.raw[res$intermediate$timeColumn.regr, 'p'], digits=digits+1));
  } else {
    cat0("0 - that is, there is no variance among the measurements!");
  }
  if (!is.null(res$intermediate$notice)) {
    cat0("\n\n", res$intermediate$notice, "\n");
  }
  if (printPlot) {
    grid.draw(res$output$combinedPlot);
  }
  invisible();
}