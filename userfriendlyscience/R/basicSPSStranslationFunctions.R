getData <- function(filename=NULL,
                    errorMessage = "[defaultErrorMessage]",
                    use.value.labels=FALSE,
                    to.data.frame=TRUE,
                    stringsAsFactors=FALSE, ...) {
  
  dots <- list(...);
  fullArguments <- as.list(environment());
  matchedCall <- match.call();
  fullCall <- capture.output(print(matchedCall));
  filenameArgument <- filename;
  
  ### File formats that have been implemented
  supportedFormats <- c(".sav", ".csv", ".tsv", ".rda", ".ods", ".xls", "xlsx");
  
  ### Set error message
  errorMessage <- sub("\\[defaultErrorMessage\\]",
                      paste0("Specified file does not exist or does not have an extension identifying ",
                             "it as a readable filetype (valid extensions are: '",
                             paste(supportedFormats, collapse="', '"), "')."),
                      errorMessage);
  
  if (is.null(filename)) {
    ### If no filename is specified, request one from the user
    cat("You did not specify a file to open. Therefore, please select the",
        "file to open in the File Selection Dialog.",
        "Note that this dialog can sometimes appear behind the R window.",
        "If you do not see the file dialog now, use ALT-TAB or check the ",
        "start bar (in Windows), use COMMAND-TAB (in OSX), or check the ",
        "dock (in *nux based systems such as",
        "Ubuntu or OS X).");
    filename <- file.choose();
    slashesFilename <- gsub("\\", "/", filename, fixed=TRUE);
    
    if (length(matchedCall) == 1) {
      filenameArgument <- sub('getData(', paste0('getData(filename="',
                                                 trim(slashesFilename), '"'),
                              fullCall, fixed=TRUE);
    } else {
      filenameArgument <- sub('getData(', paste0('getData(filename="',
                                                 trim(slashesFilename), '", '),
                              fullCall, fixed=TRUE);
    }
    
    filenameArgument <- gsub(", ", ",\n        ", filenameArgument, fixed=TRUE);
    
    cat("\n\nYou have selected a file. Based on your call and the filename",
        "and directory (path) you selected, this is the",
        "command you can use to read the datafile without",
        "a dialog, for example in an R script:\n\n");
    cat(filenameArgument, ";\n\n", sep="");
  }
  
  extension <- tolower(substring(filename, nchar(filename) - 3));

  if (!file.exists(filename) |
        !(extension %in% supportedFormats)) {
    ### Show error if the file doesn't exist or has the wrong extension
    stop(errorMessage);
  }
  else {
    if (extension == ".rda") {
      dat <- load(filename);
      dat <- get(dat);
    }
    if (extension == ".sav") {
      dat <- suppressWarnings(read.spss(filename, use.value.labels=use.value.labels,
                              to.data.frame=to.data.frame, ...));
      
#       dat <- read.spss(filename, use.value.labels=use.value.labels,
#                        to.data.frame=to.data.frame, ...);
#       cat("Note that a warning like 'Unrecognized record type 7, subtype ## encountered in system file'",
#           "is no cause for concern; the file is read normally.\n");
      
#       dat <- tryCatch({
#         read.spss(filename, use.value.labels=use.value.labels,
#                   to.data.frame=to.data.frame, ...);
#       }, warning=function(w) {
#         if (grepl("Unrecognized record type 7, subtype [0123456789]+ encountered in system file", w)) {
#           return(suppressWarnings(read.spss(filename, use.value.labels=use.value.labels,
#                                   to.data.frame=to.data.frame, ...)));
#          }
#          else {
#            return(read.spss(filename, use.value.labels=use.value.labels,
#                             to.data.frame=to.data.frame, ...));
#          }
#       });
      
    }
    else if (extension == ".csv") {
      dat <- read.csv(filename, stringsAsFactors=stringsAsFactors, ...);
    }
    else if (extension == ".tsv") {
      dat <- read.delim(filename, stringsAsFactors=stringsAsFactors, ...);
    }
    else if (extension == ".ods") {
      
      stop("Sorry, I currently do not know how to import OpenOffice files. If you do, ",
           "please contact me and I'll add this as well!\nOf course, you can always export from ",
           "LibreOffice or OpenOffice to .csv (comma separated values) and load that file.");
      
#       if (!is.element('ROpenOffice', installed.packages()[, 1])) {
#          stop("To load OpenOffice or LibreOffice files, I need package 'ROpenOffice', ",
#               "which is not on CRAN. Please visit http://omegahat.org for instructions, ",
#               "or you can try to downloads and install it yourself directly using:\n\n",
#               "install.packages('ROpenOffice', repos = 'http://www.omegahat.org/R', type = 'source');\n\n",
#               "Note that you might need specific tools to compile this source package ",
#               "(see Details in the install.packages() help, displayed with:\n\n?install.packages;");
#       }
#       require('ROpenOffice');
#       dat <- read.ods(filename, ...);
    }
    else if ((extension == ".xls") || (extension == "xlsx")) {
      if (!is.element('XLConnect', installed.packages()[, 1])) {
        stop("To load Excel (.xls or .xlsx) files, I need package 'XLConnect', ",
             "which in turn requires Java. Please install it yourself if you wish to ",
             "use this. You can install it using:\n\n",
             "install.packages('XLConnect')\n\nOf course, you can always export from ",
             "Excel to .csv (comma separated values) and load that file.");
      }
      else {
        wb <- XLConnect::loadWorkbook(filename, ...);
        dat <- XLConnect::readWorksheet(wb, sheet=1);
        if (requireNamespace('XLConnect')) {
          wb <- XLConnect::loadWorkbook(filename, ...);
          dat <- XLConnect::readWorksheet(wb, sheet=1);
        } else {
          stop("To load Excel (.xls or .xlsx) files, I need package 'XLConnect', ",
               "which in turn requires Java. Please install it yourself if you wish to ",
               "use this. You can install it using:\n\n",
               "install.packages('XLConnect')\n\nOf course, you can always export from ",
               "Excel to .csv (comma separated values) and load that file.");
        }
     }
   }
    
    ### Store the file where we got this dataframe
    attr(dat, "fileName") <- filename;
    ### Store the call
    attr(dat, "getDataCall") <- filenameArgument;
    
    ### Return the resuls
    return(dat);
  }  
}

getDat <- function(dfName="dat", backup=TRUE, ...) {
  dat <- getData(...);
  if (exists(dfName, envir=sys.frame(-1)) && backup) {
    backupName <- paste0(dfName, '_backup_',
                         format(Sys.time(), "%Y%m%d_%H%M%S"));
    assign(backupName,
           value=get(dfName, envir=sys.frame(-1)),
           envir=sys.frame(-1));
    warning("NOTE: an object called '", dfName, "' already existed; I renamed ",
            "it to '", backupName, "'.");
  }
  assign(dfName, value=dat, envir=sys.frame(-1));
  cat("The data has been stored in a dataframe called '",
      dfName, "'. That means that if you want to repeat this command and ",
      "store the dataframe with the same name, you have to use:\n\n",
      dfName, " <- ", attributes(dat)$getDataCall, ";\n\n",      
      sep="");
}


exportToSPSS <- function (dat, datafile, codefile, fileEncoding = "UTF-8",
                          newLinesInString = " |n| ") {
  
  ### Convert newline characters to spaces
  if (any(charVectors <- sapply(dat, is.character))) {
    dat[, charVectors] <- data.frame(lapply(dat[, charVectors],
                                           function(x) {
                                             return(gsub('\n', newLinesInString,
                                                         x));
                                           }), stringsAsFactors=FALSE);
  }
  
  ### Export datafile
  write.table(massConvertToNumeric(dat), file = datafile,
              row.names = FALSE, col.names = TRUE, 
              sep = ",", quote = TRUE, na = "",
              fileEncoding = fileEncoding);

  codeFileConnection=file(codefile, open="w", encoding=fileEncoding);
  
  cat(paste0("GET DATA
  /TYPE = TXT
  /FILE = \"", datafile, "\"
  /DELIMITERS = \",\"
  /QUALIFIER = '\"'
  /FIRSTCASE = 2
  /VARIABLES =\n"), file=codeFileConnection);
  
  varlabels = names(dat);
  varnames = gsub("[^[:alnum:]_\\$@#]", "\\.", names(dat));
  
  cat(paste0("  ", varnames, " ",
             unlist(lapply(dat, function(x) {
               if (is.character(x)) {
                 return(paste0('A', max(nchar(x))));
               } else {
                 return("F8.2");
               }
             })), collapse="\n"), file=codeFileConnection, append=TRUE);

  cat(".\n\nVARIABLE LABELS\n", file = codeFileConnection, append = TRUE);
  
  cat(paste(varnames,
            paste("\"", varlabels, "\"", sep = ""),
            "\n"), ".\n", file = codeFileConnection, 
      append = TRUE);
  
  factors <- sapply(dat, is.factor);
  
  if (any(factors)) {
    cat("\nVALUE LABELS\n", file = codeFileConnection, append = TRUE);
    for (v in which(factors)) {
      cat("/\n", file = codeFileConnection, append = TRUE);
      cat(varnames[v], " \n", file = codeFileConnection, append = TRUE, 
          sep = "");
      levs <- levels(dat[[v]]);
      cat(paste(seq_along(levs),
                paste("\"", levs, "\"", sep = ""),
                "\n", sep = " "), 
          file = codeFileConnection, append = TRUE);
    }
    cat(".\n", file = codeFileConnection, append = TRUE);
  }
  
  cat("\nEXECUTE.\n", file = codeFileConnection, append = TRUE);
  
  close(codeFileConnection);
}

mediaan <- function(vector) {
  if (is.data.frame(vector) | is.matrix(vector)) {
    stop("The first argument is not a vector! If you need to specify ",
         "a variable from a dataframe, separate the name of the ",
         "dataframe and the variable name with a dollar sign, for ",
         "example using 'dat$gender' to extract variable 'gender' from ",
         "dataframe 'dat'.");
  }
  if (is.character(vector)) {
    stop('The first argument is a character vector; please convert it ',
         'to a factor or a numeric vector first.');
  }
  ### Store original class
  originalClass <- class(vector);
  ### Store original vector
  originalVector <- vector;
  ### Convert to numeric vector
  vector <- as.numeric(vector);
  ### If need be, convert to relevant category
  if ("factor" %in% originalClass) {
    levelIndex <- median(vector, na.rm=TRUE);
    if (round(levelIndex) == levelIndex) {
      res <- levels(originalVector)[median(vector, na.rm=TRUE)];
    }
    else {
      res <- c(levels(originalVector)[round(median(vector, na.rm=TRUE)-.5)],
               levels(originalVector)[round(median(vector, na.rm=TRUE)+.5)]);
    }
  }
  else {
    res <- median(vector, na.rm=TRUE);
  }
  return(res);
}

modus <- function(vector) {
  if (is.data.frame(vector) | is.matrix(vector)) {
    stop("The first argument is not a vector! If you need to specify ",
         "a variable from a dataframe, separate the name of the ",
         "dataframe and the variable name with a dollar sign, for ",
         "example using 'dat$gender' to extract variable 'gender' from ",
         "dataframe 'dat'.");
  }
  ### Store original class
  originalClass <- class(vector);
  ### Convert to factor
  vector <- as.factor(vector);
  ### Store frequencies
  freqs <- summary(vector);
  ### Determine highest frequency
  highestFreq <- max(freqs);
  ### Store the names of the most common category (or categories)
  categoryVector <- names(freqs[freqs==highestFreq]);
  ### Now, we need to supply this back in the same class as the original.
  if (originalClass=="factor") {
    categoryVector <- as.factor(categoryVector);
  }
  else {
    class(categoryVector) <- originalClass;
  }
   return(categoryVector);
}

filterBy <- function(dat, expression,
                     replaceOriginalDataframe=TRUE,
                     envir = parent.frame()) {
  ### Store original dataframe and current time
  originalDataframeName <- as.character(substitute(dat));
  currentTime <- Sys.time();
  timeStamp <- round(as.numeric(currentTime) * 100);
  newDataframeName <- paste0('.', originalDataframeName, "_at_", timeStamp);
  
  ### Store original dataframe with new name in parent environment
  assign(newDataframeName, value=dat, envir=envir);
  
  ### Store number of rows for reporting to user
  nrOfRows <- nrow(dat);
  
  if (!is.logical(expression)) {
    if (is.character(expression)) {
      ### Replace single 'equals' characters with the 'equals' operator 
      expression <- gsub("([^=])=([^=])", "\\1==\\2", expression);
      ### Generate logical vector
      expression <- with(dat, eval(parse(text=expression)));
    }
    else {
      stop("The argument 'expression' must be either a logical vector or a character string with a logical expression!");
    }
  }
  
  ### Create filtered dataframe
  dat <- dat[expression, ];
  
  attr(dat, "originalDataframeName") <- originalDataframeName;
  attr(dat, "lastUnfilteredDataframeName") <- newDataframeName;
  attr(dat, "lastUnfilteredDataframeEnvir") <- envir;
  attr(dat, "lastFiltering") <- currentTime;
  
  cat("Filtered ", nrOfRows - nrow(dat) ,
      " rows (records, cases, participants, or datapoints) from dataframe '",
      originalDataframeName, "'; result has ", nrow(dat), " rows.\n", sep="");
  
  if (replaceOriginalDataframe) {
    assign(originalDataframeName, value=dat, envir=sys.frame(-1));
    invisible(dat);
  }
  else {
    return(dat);
  }

}

useAll <- function(dat, replaceFilteredDataframe = TRUE) {
  ### Store name of filtered dataframe
  filteredDataframeName <- as.character(substitute(dat));
  ### Store number of rows in filtered dataframe
  nrOfRows <- nrow(dat);
  
  ### Get information required to find original dataframe
  originalDataframeName <- attr(dat, "originalDataframeName");
  lastUnfilteredDataframeName <- attr(dat, "lastUnfilteredDataframeName");
  lastUnfilteredDataframeEnvir <- attr(dat, "lastUnfilteredDataframeEnvir");
  lastFiltering <- attr(dat, "lastFiltering");
  
  ### Check whether original exists
  if (exists(lastUnfilteredDataframeName, envir=lastUnfilteredDataframeEnvir)) {
    dat <- get(lastUnfilteredDataframeName, envir=lastUnfilteredDataframeEnvir);
    rm(list=lastUnfilteredDataframeName, envir=lastUnfilteredDataframeEnvir);
  }
  else {
    stop("Could not find the original, prefiltered version of the dataframe (which was stored as '",
         lastUnfilteredDataframeName, " in environment '", lastUnfilteredDataframeEnvir,"').");
  }
  
  cat("Removed last applied filter to dataframe '", filteredDataframeName, "', which was ",
      "applied at ", format(lastFiltering), " and removed (filtered) ",
      nrow(dat) - nrOfRows, " rows (records, cases, participants, or datapoints) ",
      "from the dataframe that was originally called '", originalDataframeName,
      "'. Restored dataframe has ", nrow(dat), " rows.\n", sep="");
  
  if (replaceFilteredDataframe) {
    assign(filteredDataframeName, value=dat, envir=sys.frame(-1));
    cat("Replaced filtered dataframe '", filteredDataframeName, "'.\n", sep="");
    invisible(dat);
  }
  else {
    return(dat);
  }
  
}
