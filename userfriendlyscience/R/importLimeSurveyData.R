importLimeSurveyData <- function(datafile = NULL,
                                 scriptfile = NULL,
                                 limeSurveyRegEx.varNames =
                                   "names\\(data\\)\\[\\d*\\] <- ",
                                 limeSurveyRegEx.toChar =
                                   "data\\[, \\d*\\] <- as.character\\(data\\[, \\d*\\]\\)",
                                 limeSurveyRegEx.varLabels =
                                   "attributes\\(data\\)\\$variable.labels\\[\\d*\\] <- \".*\"",
                                 limeSurveyRegEx.toFactor =
                                   paste0("data\\[, \\d*\\] <- factor\\(data\\[, \\d*\\], ",
                                          "levels=c\\(.*\\),labels=c\\(.*\\)\\)"),
                                 limeSurveyRegEx.varNameSanitizing =
                                   list(list(pattern = "#", replacement = "."),
                                        list(pattern = "\\$", replacement = ".")),
                                 setVarNames = TRUE,
                                 setLabels = TRUE,
                                 convertToCharacter = FALSE,
                                 convertToFactor = FALSE,
                                 categoricalQuestions = NULL,
                                 massConvertToNumeric = TRUE,
                                 dataHasVarNames = TRUE,
                                 encoding='UTF-8-BOM') {
  
  ### Load datafile
  if (dataHasVarNames) {
    data <- getData(datafile, quote = "'\"", na.strings=c("", "\"\""),
                    stringsAsFactors=FALSE, encoding=encoding, header=TRUE);
  } else {
    data <- getData(datafile, quote = "'\"", na.strings=c("", "\"\""),
                    stringsAsFactors=FALSE, encoding=encoding, header=FALSE);
  }
  
  ### Load scriptfile
  if (!is.null(scriptfile)) {
    if (!file.exists(scriptfile)) {
      stop("File specified as scriptfile ('", scriptfile, "') does not exist!");
    }
    datascript <- readLines(scriptfile, encoding=encoding);
    varNamesScript <- datascript[grepl(limeSurveyRegEx.varNames,
                                       datascript)];
    varLabelsScript <- datascript[grepl(limeSurveyRegEx.varLabels,
                                        datascript)];
    toCharScript <- datascript[grepl(limeSurveyRegEx.toChar,
                                     datascript)];
    toFactorScript <- datascript[grepl(limeSurveyRegEx.toFactor,
                                       datascript)];
    if (setVarNames) {
      eval(parse(text=varNamesScript));
    }
    if (setLabels) {
      eval(parse(text=varLabelsScript));
    }
    if (convertToCharacter) {
      eval(parse(text=toCharScript));
    }
    if (convertToFactor || (!is.null(categoricalQuestions))) {
      if (massConvertToNumeric) {
        data <- massConvertToNumeric(data);
      }
      if (!is.null(categoricalQuestions)) {
        if (setVarNames) {
          varNames <- names(data);
        } else {
          stop("You can't set setVarNames to FALSE and also set ",
               "categoricalQuestions to anything else than NULL, ",
               "because the content of categoricalQuestions should ",
               "be the LimeSurvey variables names!");
        }
        toFactorScript <- unlist(lapply(as.list(categoricalQuestions),
                                        function(x, string=toFactorScript,
                                                 varNms=varNames) {
                                          return(grep(paste0("data\\[, ",
                                                             which(varNms==x),
                                                             "\\] <-"),
                                                      string, value=TRUE));
                                        }));
      }
      eval(parse(text=toFactorScript));
    }
  } else {
    if (massConvertToNumeric) {
      data <- massConvertToNumeric(data);
    }
  }
  if (length(limeSurveyRegEx.varNameSanitizing)) {
    for (currentRegexPair in limeSurveyRegEx.varNameSanitizing) {
      names(data) <- gsub(currentRegexPair$pattern,
                          currentRegexPair$replacement,
                          names(data));
    }
  }
  return(data);
}