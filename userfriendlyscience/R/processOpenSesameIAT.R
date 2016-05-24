###########################################################
###########################################################
###
### Function to read OpenSesame IAT data files (i.e. a
### collection of .csv files), process these, compute the
### D600, and produce a wide dataframe with the output
###
### For a description of the procedure, see e.g.
### http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3769683/
###
### File created by Frederik van Acker and Gjalt-Jorn
### Peters. Questions? You can contact us through
### http://behaviorchange.eu
###
### This file is licensed under Creative Commons BY-SA 3.0
### (Attribution-ShareAlike, which means that you can
### freely use and distribute this file, and you're
### allowed to alter it as long as you release the edited
### version using the same license (i.e. again freely
### available). This license is used to promote Open
### Science and Full Disclosure. For the complete
### license, see http://creativecommons.org/licenses/by-sa/3.0/deed.en_US
### For more information about Full Disclosure, see
### http://sciencerep.org/fulldisclosure
###
###########################################################
###########################################################

###########################################################
### Secondary function definitions
###########################################################

### This function adds a message to the log
### and displays it depending on a boolean
addToLog <- function(fullLog, ..., showLog = FALSE) {
  if (showLog) {
    cat(paste0(...))
  }
  return(paste0(fullLog, ...));
}

###########################################################
### Main function
###########################################################

processOpenSesameIAT <- function(dataPath,
                                 blocks.sizes = c(18, 36, 48, 36, 48),
                                 blocks.congruent = c(2, 3),
                                 blocks.incongruent = c(4, 5),
                                 blocks.realTrials = c(3, 5),
                                 blocks.practiceTrials = c(2, 4),
                                 responseTime.min = 400,
                                 responseTime.max = 2500,
                                 responseTime.penalty = 600,
                                 outputFile = NULL,
                                 wideOutputFile = NULL,
                                 showLog = TRUE,
                                 filenameRegEx = "subject-(\\d+)(\\w+)\\.csv",
                                 regExValues = c("subject","session"),
                                 participantVarName = "subject",
                                 taskVarName = "session",
                                 openSesameVarNames = list(correct = "correct",
                                                           response_time = "response_time"),
                                 roundOutput = 6,
                                 decimalSeparator = ".") {  
  
  if (!is.null(participantVarName) && !is.null(regExValues) &&
        !(participantVarName %in% regExValues)) {
    stop(paste0("The participantVarName ('", participantVarName,
                "') is not included in the regExValues vector!"));
  }
  if (!is.null(taskVarName) && !is.null(regExValues) &&
        !(taskVarName %in% regExValues)) {
    stop(paste0("The taskVarName ('", taskVarName,
                "') is not included in the regExValues vector! ",
                "Note that you have to set it to NULL if you only ",
                "want to specify a participantVarName (by default, ",
                "when not specified, taskVarName = 'session')."));
  }
  
  ### Generate object to return
  res <- list(input = as.list(environment()));
  
  ### Store input variables
  res$dataPath <- dataPath;
  res$blocks.sizes <- blocks.sizes;
  res$blocks.congruent <- blocks.congruent;
  res$blocks.incongruent <- blocks.incongruent;
  res$blocks.realTrials <- blocks.realTrials;
  res$blocks.practiceTrials <- blocks.practiceTrials;
  res$responseTime.min <- responseTime.min;
  res$responseTime.max <- responseTime.max;
  res$responseTime.penalty <- responseTime.penalty;
  res$outputFile <- outputFile;

  ### Store object with each participants' datafile
  res$file.raw <- list();
  res$file.clean <- list();
  
  ### Store object with each participants' dataframe
  res$dat.raw <- list();
  
  ### Create a vector with the number of the block
  ### for each trial
  res$blocks.vector <- rep(1:length(blocks.sizes), blocks.sizes);

  ### Create a vector indicating whether each trial
  ### is congruent or not.
  res$congruency.vector <- ifelse(res$blocks.vector %in% res$blocks.congruent,
                                  TRUE,
                                  ifelse(res$blocks.vector %in% res$blocks.incongruent,
                                         FALSE,
                                         NA));
  
  res$blocks.sizes <- blocks.sizes;
  res$blocks.congruent <- blocks.congruent;
  res$blocks.incongruent <- blocks.incongruent;
  res$blocks.realTrials <- blocks.realTrials;
  res$blocks.practiceTrials <- blocks.practiceTrials;
  res$responseTime.min <- responseTime.min;
  res$responseTime.max <- responseTime.max;
  res$responseTime.penalty <- responseTime.penalty;
  
  logText <- addToLog("--- processOpenSesameIAT log ---\n",
                      paste0("Specified input:\n",
                             length(blocks.sizes), " blocks were specified,",
                             " of respectively ", vecTxt(blocks.sizes),
                             " trials each.\n", 
                             "Blocks ", vecTxt(blocks.congruent),
                             " are congruent and blocks ", vecTxt(blocks.incongruent),
                             " are incongruent.\n",
                             "Blocks ", vecTxt(blocks.realTrials),
                             " are the real trials and blocks ", vecTxt(blocks.practiceTrials),
                             " are the practice trials.\n",
                             "Response times below ", responseTime.min,
                             " milliseconds will be deleted, response times over ",
                             responseTime.max, " will be replaced with ",
                             responseTime.max, " and a penalty of ",
                             responseTime.penalty, " milliseconds will be ",
                             "added to incorrect responses.\n\n",
                             "In total, this design consists of ", sum(blocks.sizes),
                             " trials. That means that every datafile will need to ",
                             "have exactly ", sum(blocks.sizes), " valid rows (i.e. ",
                             "excluding the header and empty rows); datafiles ",
                             "with a different number of rows will be excluded.",
                             "\n\n"),
                             showLog=showLog);
                      
  if (!is.null(wideOutputFile) && is.null(taskVarName)) {
    logText <- addToLog(logText,
                        paste0("Note: you are requesting a wide output file, ",
                               "but you did not specify a taskVarName, so I ",
                               "don't know how to identify different datafiles for ",
                               "each participant. No wide datafile will be produced.\n\n"),
                        showLog=showLog);
  }
  
  ### Store filenames in a character vector.
  res$inputFiles.all <- list.files(dataPath);
  
  ### Check whether a regular expression was specified to select
  ### filenames and information about each participant/sessions
  if (!is.null(filenameRegEx)) {
    ### Select filenames matching regular expression
    res$inputFiles <- grep(filenameRegEx, res$inputFiles.all, value=TRUE);
    ### Start log by indicating number of valid files
    logText <- addToLog(logText,
                        paste0("Read folder '", dataPath, "'. Out of ",
                               length(res$inputFiles.all), " files, ",
                               length(res$inputFiles), " matched regular expression '",
                               filenameRegEx, "' and will be processed.\n"),
                        showLog=showLog);
  }
  else {
    ### Process all files
    res$inputFiles <- res$inputFiles.all;
    ### Start log by indicating number of valid files
    logText <- addToLog("--- processOpenSesameIAT log ---\n",
                        paste0("Read folder '", dataPath, "'. Processing all ",
                               length(res$inputFiles.all), " files.\n"),
                        showLog=showLog);
  }
  
  ### Create dataframe for aggregated results
  res$dat <- data.frame();
  
  ### Start loop to read each participants' file
  for (currentParticipant in 1:length(res$inputFiles)) {
    
    logText <- addToLog(logText,
             paste0("Starting to process file ", res$inputFiles[currentParticipant], "\n"),
             showLog=showLog);

    ### Store filename that was used
    res$dat[currentParticipant, 'filename'] <-
      res$inputFiles[currentParticipant];
    
    ### If a regular expression was set, extract information
    ### from the filename and store this as well
    if (!is.null(filenameRegEx)) {
      for (currentVariableIndex in 1:length(regExValues)) {
        res$dat[currentParticipant, regExValues[currentVariableIndex]] <-
          sub(filenameRegEx,
              paste0("\\", currentVariableIndex),
              res$inputFiles[currentParticipant]);
      }
    }
    
    ### For some reason, some datafiles have the header
    ### repeated in between the data. Therefore, we read
    ### the file manually, and then scan for lines that
    ### are duplicates and remove them.
    res$file.raw[[currentParticipant]] <-
      readLines(file.path(dataPath, res$inputFiles[currentParticipant]));
    ### Remove duplicates
    res$file.clean[[currentParticipant]] <-
      res$file.raw[[currentParticipant]][
        !duplicated(res$file.raw[[currentParticipant]])
      ];
    
    ### Also, for some reason, some datafiles use a single
    ### double quote character to quote all fields; but
    ### sometimes, this pattern gets screwed up, causing
    ### double double quotes to appear. This leads to
    ### problems when reading the file, and since there
    ### are no text fields in the file, just remove the
    ### quotes from every line.
    res$file.clean[[currentParticipant]] <-
      gsub('"', '', res$file.clean[[currentParticipant]]);
    
    ### Store number of rows for this participant
    res$dat[currentParticipant, 'nrOfRowsInDatafile']  <-
      length(res$file.clean[[currentParticipant]]);
    
    ### Check the number of rows in the datafile. This
    ### should be equal to the number of blocks + 1;
    ### otherwise, something unknown went wrong, rendering
    ### the datafile unreliable (i.e. exclude the
    ### participant)
    if (!(length(res$file.clean[[currentParticipant]]) ==
            (sum(res$blocks.sizes) + 1))) {
      logText <- addToLog(logText,
               paste0("    File has ",
                      length(res$file.clean[[currentParticipant]]) - 1,
                      " valid rows (trials); ", sum(res$blocks.sizes),
                      " required; excluding participant.\n"),
               showLog=showLog);
      ### Store 'invalidness' in dataframe
      res$dat[currentParticipant, 'validDatafile']  <- 0;
    }
    else {
      ### Store 'validness' in dataframe
      res$dat[currentParticipant, 'validDatafile']  <- 1;
      
      nrOfElements <- unlist(lapply(res$file.clean[[currentParticipant]],
                             function(x) { return (length(strsplit(x, ",")[[1]])); }));

      if (nrOfElements[1] < max(tail(nrOfElements, -1))) {
        logText <- addToLog(logText,
                            paste0("    Error: the first line (the header, ",
                                   "containing the variable names) only contains ",
                                   nrOfElements[1],
                                   " elements, whereas the other lines contain ",
                                   "between ",
                                   paste0(range(tail(nrOfElements, -1)), collapse=" and "),
                                   " elements; excluding participant.\n"),
                            showLog=showLog);
      } else {
      
        ### Parse participants' datafile
        res$dat.raw[[currentParticipant]] <-
          read.csv(text=res$file.clean[[currentParticipant]],
                   header = TRUE);
   
        ### Add vector with block numbers for each trial
        res$dat.raw[[currentParticipant]]$blockNumber <- res$blocks.vector;
    
        ### Add vector with congruency for each trial
        res$dat.raw[[currentParticipant]]$congruent <- res$congruency.vector;
  
        ### Correct the 'correct' answers, such that anything
        ### other than 1 is replaced by 0
        res$dat.raw[[currentParticipant]]$correct.corrected <-
          ifelse(res$dat.raw[[currentParticipant]][[openSesameVarNames$correct]] == 0,
                 0, 1);
        res$dat.raw[[currentParticipant]]$correct.corrected <-
          as.numeric(res$dat.raw[[currentParticipant]]$correct.corrected);
        
        ### For each trial, store whether the response
        ### time is sufficiently high
        res$dat.raw[[currentParticipant]]$rt_highEnough <-
          res$dat.raw[[currentParticipant]][[openSesameVarNames$response_time]] >= res$responseTime.min;
        
        logText <- addToLog(logText,
                 paste0("    ",
                        sum(res$dat.raw[[currentParticipant]]$rt_highEnough),
                        " trials with response time > ",
                        res$responseTime.min, "\n"),
                 showLog=showLog);
        
        ### For each trial, store whether the response
        ### time is sufficiently low
        res$dat.raw[[currentParticipant]]$rt_lowEnough <-
          res$dat.raw[[currentParticipant]][[openSesameVarNames$response_time]] <= res$responseTime.max;
    
        logText <- addToLog(logText,
                 paste0("    ",
                        sum(res$dat.raw[[currentParticipant]]$rt_lowEnough),
                        " trials with response time < ",
                        res$responseTime.max, "\n"),
                 showLog=showLog);
    
        ### For each trial, store whether the response
        ### time is valid
        res$dat.raw[[currentParticipant]]$rt_valid <-
          res$dat.raw[[currentParticipant]]$rt_lowEnough &
          res$dat.raw[[currentParticipant]]$rt_highEnough;
  
        logText <- addToLog(logText,
                 paste0("    ",
                        sum(res$dat.raw[[currentParticipant]]$rt_valid),
                        " trials with valid response time\n"),
                 showLog=showLog);
        
        ### Generate a new variable containing the response times,
        ### but where response times exceeding the highest acceptable
        ### value are replaced with max response time
        res$dat.raw[[currentParticipant]]$response_time_clipped <- ifelse(
          res$dat.raw[[currentParticipant]]$rt_lowEnough,
          res$dat.raw[[currentParticipant]][[openSesameVarNames$response_time]],
          res$responseTime.max);
        
        ### Compute the mean response times per block,
        ### but only for correct responses, and based
        ### on the clipped response times (i.e. where
        ### unacceptably high response times are replaced
        ### with the max response time).
        ### Then repeat each mean response time with
        ### the number of sessions in that block using rep
        ### and store the result in a new variable in the
        ### dataframe, so that for every trial, we also
        ### have the mean for the relevant block available
        
        ### However, some people have no (zero) correct
        ### responses within a block, which will cause an
        ### error. Therefore, check that first.
        if (length(aggregate(response_time_clipped ~ blockNumber, data =
                           res$dat.raw[[currentParticipant]][
                             res$dat.raw[[currentParticipant]]$rt_valid &
                             (res$dat.raw[[currentParticipant]]$correct.corrected == 1)
                             , ],
                         FUN = "mean", na.rm=TRUE)$response_time_clipped) < length(blocks.sizes)) {
          ### In this case, one or more blocks has zero correct responses.
          ### That means the participant has to be excluded.
          logText <- addToLog(logText,
                   paste0("    This participant has zero (0) correct ",
                          "responses in one or more blocks; excluding participant.\n"),
                   showLog=showLog);
        }
        else {
          
          ### Compute and store means per block
          res$dat.raw[[currentParticipant]]$response_time_blockMean <- 
            rep(aggregate(response_time_clipped ~ blockNumber, data =
                          res$dat.raw[[currentParticipant]][
                            res$dat.raw[[currentParticipant]]$rt_valid &
                            (res$dat.raw[[currentParticipant]]$correct.corrected == 1)
                          , ],
                          FUN = "mean", na.rm=TRUE)$response_time_clipped, blocks.sizes);
  
          ### Replace the response times for incorrect
          ### answers with the mean response time for that
          ### block plus the penalty
          res$dat.raw[[currentParticipant]]$response_time_penalized <-
            ifelse(res$dat.raw[[currentParticipant]]$correct.corrected == 1,
                   res$dat.raw[[currentParticipant]]$response_time_clipped,
                   res$dat.raw[[currentParticipant]]$response_time_blockMean +
                     res$responseTime.penalty);

          logText <- addToLog(logText,
                              paste0("    ",
                                     sum(res$dat.raw[[currentParticipant]]$correct.corrected != 1),
                                     " trials were incorrect: setting block mean plus penalty as response time.\n"),
                              showLog=showLog);          
          
          ### Compute standard deviation, first for all blocks;
          ### then only real trials, then only practice trials.
          res$dat[currentParticipant, 'sd_all'] <-
            sd(res$dat.raw[[currentParticipant]]
               [res$dat.raw[[currentParticipant]]$rt_highEnough &
                res$dat.raw[[currentParticipant]]$rt_lowEnough,
                'response_time_penalized']);
           res$dat[currentParticipant, 'sd_real'] <-
             sd(res$dat.raw[[currentParticipant]]
                [res$dat.raw[[currentParticipant]]$rt_highEnough &
                 res$dat.raw[[currentParticipant]]$rt_lowEnough &
                 res$dat.raw[[currentParticipant]]$blockNumber %in% res$blocks.realTrials,
                 'response_time_penalized']);
           res$dat[currentParticipant, 'sd_practice'] <-
             sd(res$dat.raw[[currentParticipant]]
                [res$dat.raw[[currentParticipant]]$rt_highEnough &
                 res$dat.raw[[currentParticipant]]$rt_lowEnough &
                 res$dat.raw[[currentParticipant]]$blockNumber %in% res$blocks.practiceTrials,
                 'response_time_penalized']);
          
          ### Compute means for congruent blocks
          res$dat[currentParticipant, 'mean_congruent_all'] <-
            mean(res$dat.raw[[currentParticipant]]
                 [res$dat.raw[[currentParticipant]]$rt_highEnough &
                  res$dat.raw[[currentParticipant]]$rt_lowEnough &
                  res$dat.raw[[currentParticipant]]$congruent,
                  'response_time_penalized'], na.rm=TRUE);
           res$dat[currentParticipant, 'mean_congruent_real'] <-
             mean(res$dat.raw[[currentParticipant]]
                  [res$dat.raw[[currentParticipant]]$rt_highEnough &
                   res$dat.raw[[currentParticipant]]$rt_lowEnough &
                   res$dat.raw[[currentParticipant]]$congruent &
                     res$dat.raw[[currentParticipant]]$blockNumber %in% res$blocks.realTrials,
                   'response_time_penalized']);
           res$dat[currentParticipant, 'mean_congruent_practice'] <-
             mean(res$dat.raw[[currentParticipant]]
                  [res$dat.raw[[currentParticipant]]$rt_highEnough &
                   res$dat.raw[[currentParticipant]]$rt_lowEnough &
                   res$dat.raw[[currentParticipant]]$congruent &
                     res$dat.raw[[currentParticipant]]$blockNumber %in% res$blocks.practiceTrials,
                   'response_time_penalized']);
          
          ### And for incongruent blocks
          res$dat[currentParticipant, 'mean_incongruent_all'] <-
            mean(res$dat.raw[[currentParticipant]]
                 [res$dat.raw[[currentParticipant]]$rt_highEnough &
                  res$dat.raw[[currentParticipant]]$rt_lowEnough &
                  !res$dat.raw[[currentParticipant]]$congruent,
                  'response_time_penalized'], na.rm=TRUE);
           res$dat[currentParticipant, 'mean_incongruent_real'] <-
             mean(res$dat.raw[[currentParticipant]]
                  [res$dat.raw[[currentParticipant]]$rt_highEnough &
                   res$dat.raw[[currentParticipant]]$rt_lowEnough &
                   !res$dat.raw[[currentParticipant]]$congruent &
                     res$dat.raw[[currentParticipant]]$blockNumber %in% res$blocks.realTrials,
                   'response_time_penalized']);
           res$dat[currentParticipant, 'mean_incongruent_practice'] <-
             mean(res$dat.raw[[currentParticipant]]
                  [res$dat.raw[[currentParticipant]]$rt_highEnough &
                   res$dat.raw[[currentParticipant]]$rt_lowEnough &
                   !res$dat.raw[[currentParticipant]]$congruent &
                     res$dat.raw[[currentParticipant]]$blockNumber %in% res$blocks.practiceTrials,
                   'response_time_penalized']);
          
          ### Compute D600, which is simply the difference between
          ### the means divided by the standard deviation; again,
          ### repeat this for all three groups (all, real, and
          ### practice).
          res$dat[currentParticipant, 'd600_all'] <-
            (res$dat[currentParticipant, 'mean_congruent_all'] -
             res$dat[currentParticipant, 'mean_incongruent_all']) /
            res$dat[currentParticipant, 'sd_all'];
          res$dat[currentParticipant, 'd600_real'] <-
            (res$dat[currentParticipant, 'mean_congruent_real'] -
             res$dat[currentParticipant, 'mean_incongruent_real']) /
             res$dat[currentParticipant, 'sd_real'];
          res$dat[currentParticipant, 'd600_practice'] <-
            (res$dat[currentParticipant, 'mean_congruent_practice'] -
             res$dat[currentParticipant, 'mean_incongruent_practice']) /
             res$dat[currentParticipant, 'sd_practice'];
          
          ### Store number of trials that had response times
          ### below or above the acceptable bandwidth, as well
          ### as in between, for each block
          res$dat[currentParticipant, 'rt_tooLow'] <-
            sum(!res$dat.raw[[currentParticipant]]$rt_lowEnough);
          res$dat[currentParticipant, 'rt_tooHigh'] <-
            sum(!res$dat.raw[[currentParticipant]]$rt_highEnough);
          res$dat[currentParticipant, 'rt_valid'] <-
            sum(res$dat.raw[[currentParticipant]]$rt_valid);
      
        }
      }        
      logText <- addToLog(logText,
               paste0("Done processing file ",
                      res$inputFiles[currentParticipant], "\n"),
               showLog=showLog);
    }
  }
  
  if (sum(res$dat$validDatafile) == 0) {
    
    logText <- addToLog(logText,
                        paste0("None of the datafiles is valid! Are you ",
                               "sure you specified the correct number ",
                               "of blocks and trials (see top of this log)?\n",
                               "Note that the 'nrOfRowsInDatafile' column ",
                               "contains the number of rows in each datafile.\n"),
                        showLog=showLog);
    
  } else {

    if (!is.null(participantVarName) && !is.null(taskVarName)) {
      ### Convert data to wide format
      res$dat.wide <- reshape(res$dat,
                              timevar=taskVarName,
                              idvar=participantVarName,
                              direction="wide", sep="_");
      logText <- addToLog(logText,
                          paste0("Generated wide version of datafile\n"),
                          showLog=showLog);
      if (!is.null(wideOutputFile)) {
        ### Round numbers
        if(!is.null(roundOutput)) {
          roundedWideDat <- res$dat.wide;
          numericVars <- names(roundedWideDat)[as.vector(lapply(roundedWideDat, class)=='numeric')];
          roundedWideDat[, numericVars] <-
            round(roundedWideDat[, numericVars], roundOutput);
        }
        ### Store wide datafile
        write.table(roundedWideDat, wideOutputFile, sep="\t", dec=decimalSeparator);
        logText <- addToLog(logText,
                            paste0("Saved wide datafile to ", outputFile, "\n"),
                            showLog=showLog);
      }
    }
    
    ### This is an attempt to do the long->wide conversion
    ### manually, but I didn't manage to get it working.
    ### However, I'm keeping it around just in case there's
    ### a need for some of this code later on.
  #   res$dat.wide <- ddply(a, "participant", sessionVar = 'session',
  #                         function(dat, sessionVar) {
  #     ### Create dataframe to return.
  #     res <- data.frame();
  #     ### Create object to store every line of the dataframe
  #     ### for this participant
  #     lines <- list();
  #     for (currentSession in levels(as.factor(dat[[sessionVar]]))) {
  #       ### Add current line from dataframe
  #       lines[[currentSession]] <- dat[dat[[sessionVar]] == currentSession, ];
  #       ### Prepend session to variable names
  #       names(lines[[currentSession]]) <-
  #         paste0(currentSession, "_", names(lines[[currentSession]]));
  #     }
  #     ### Combine dataframes and return using rbind.fill,
  #     ### from the plyr package
  #     return(rbind.fill(lines));
  #   });
    
  }

  if (!is.null(outputFile)) {
    ### Round numbers
    if(!is.null(roundOutput)) {
      roundedDat <- res$dat;
      numericVars <- names(roundedDat)[as.vector(lapply(roundedDat, class)=='numeric')];
      roundedDat[, numericVars] <-
        round(roundedDat[, numericVars], roundOutput);
    }
    
    ### Store aggregated datafile
    write.table(roundedDat, outputFile, sep="\t", dec=decimalSeparator);
    logText <- addToLog(logText,
                        paste0("Saved aggregated datafile to ", outputFile, "\n"),
                        showLog=showLog);
  }
  
  ### Store log
  res$log <- logText;
  
  ### Set classes for returned object and the log
  class(res) <- c('processOpenSesameIAT');
  class(res$log) <- c("processOpenSesameIAT.log");
  
  return(res);

}

print.processOpenSesameIAT <- function (x, ...) {
  cat("Ran succesfully - parsed", nrow(x$dat), "files.\n");
  if (!is.null(x$outputFile)) {
    cat("Stored aggregated datafile in", x$outputFile, "\n");
  }
}

print.processOpenSesameIAT.log <- function(x, ...) {
  cat(x);
}
