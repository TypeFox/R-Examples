if (getRversion() >= "2.15.1") globalVariables(c("x", "y", "f"))
#' @title Table With Summary Statistics
#'
#' @description Generates a table with summary statistics.
#'
#' @param .data The \code{data.frame} names.
#' @param var The name of variable or column.
#' @param by The name of grouping variable.
#' @param statistics The name of desired statistics.
#' @keywords Exploratory
#' @examples
#' tab = tabstat(Presidents, c("winner.height", "winner.vote", "turnout"))
#' # knitr::kable(tab, digits=2)
#' @export
`tabstat` <- function(.data, var, by = NULL, statistics = c("nnmiss", "mean", "sd")) UseMethod("tabstat")

#' @rdname tabstat
#' @export
`tabstat.default` <- function(.data, var, by = NULL, statistics = c("nnmiss", "mean", "sd")) {
	if (is.character(.data)) {
		nameData <- .data
		if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
		tempData <- eval(parse(text = .data))
	} else {
		if (is.data.frame(.data)) {
			nameData <- paste(deparse(substitute(.data)))
			tempData <- .data
		} else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
	}
	if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
	if (!is.character(var)) 	{ stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
	if (!is.null(by)) 		{
		if (any(is.na(match(by, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
		tempGrp <- tempData[by]
		for (i in (1:ncol(tempGrp))) { tempGrp[,i] <- factor(tempGrp[,i]) }
	} else { tempGrp <- rep(1, each = nrow(tempData)) }
	tempVar <- tempData[var]

availableStatistics <- c("n", "nnmiss", "nmiss", "sum", "css", "ucss", "se.skew", "se.kurtosis", "range", "iqr", "var", "sd", "se.mean",
		  "cv", "mean", "median", "mode", "min", "max", "q1", "q3", "skew", "kurtosis")

	if (is.null(statistics)) { statToCompute <- c("nnmiss", "mean", "sd") } else { statToCompute <- statistics }

	if (all(is.na(match(statToCompute, availableStatistics)))) { stop("All numerical descriptive measures enumerated is invalid. Expect one of the following: 'n', 'nnmiss', 'nmiss', 'sum', 'css', 'ucss', 'se.skew', 'se.kurtosis', 'range', 'iqr', 'var', 'sd', 'se.mean', 'cv', 'mean', 'median', 'mode', 'min', 'max', 'q1', 'q3', 'skew' or 'kurtosis'") }

statToCompute <- availableStatistics[na.omit(match(statToCompute, availableStatistics))]

	summaryTable <- NULL
	outputLabel <- NULL

	# --- compute for number of rows in the data set
	if(!is.na(match("n", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) { a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, length))) }
		if (any(is.na(a$Freq))) { a$Freq <- replace(a$Freq, attr(na.omit(a$Freq), "na.action"), 0) }
		summaryTable <- a
		colnames(summaryTable)[ncol(summaryTable)] <- "N_Obs"
		outputLabel <- c(outputLabel, "No. of Obs")
	}

	# --- compute for the number of non-missing observations
	if(!is.na(match("nnmiss", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) {
			newData <- na.omit(cbind(tempVar[i], tempGrp))
			a <- rbind(a, as.data.frame.table(tapply(newData[,1], newData[2:ncol(newData)], length)))
		}
		if (any(is.na(a$Freq))) { a$Freq <- replace(a$Freq, attr(na.omit(a$Freq), "na.action"), 0) }
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "ValidObs"
		outputLabel <- c(outputLabel, "No. of Non-Missing Obs.")
	}

	# --- compute for the number of missing observations
	if(!is.na(match("nmiss", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) {
			newData <- cbind(tempVar[i], tempGrp)
			newData <- newData[!complete.cases(newData),]
			if (!is.null(by)) a <- rbind(a, as.data.frame.table(tapply(newData[,1], newData[2:ncol(newData)], length)))
			else a <- rbind(a, data.frame(tempGrp = 1, Freq = nrow(newData)))
		}
		a$Freq <- replace(a$Freq, is.na(a[,"Freq"]), 0)
		#if (any(is.na(a$Freq))) { a$Freq <- replace(a$Freq, attr(na.omit(a$Freq), "na.action"), 0) }
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "MissObs"
		outputLabel <- c(outputLabel, "No. of Missing Obs.")
	}

	# --- compute the minimun observation
	if(!is.na(match("min", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, min, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Min"
		outputLabel <- c(outputLabel, "Minimum")
	}

	# --- compute the maximum observation
	if(!is.na(match("max", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, max, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Max"
		outputLabel <- c(outputLabel,"Maximum")
	}

	# --- compute the sum of the variable
	if(!is.na(match("sum", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, sum, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Sum"
		outputLabel <- c(outputLabel,"Sum")
	}

	# --- compute the mean of the variable
	if(!is.na(match("mean", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, mean, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Mean"
		outputLabel <- c(outputLabel,"Mean")
	}

	# --- compute the median of the variable
	if(!is.na(match("median", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stats::median, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Median"
		outputLabel <- c(outputLabel,"Median")
	}
	# --- determine the modal value of the variable
	if(!is.na(match("mode", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, data.frame(Var = i, Freq = paste(tapply(tempVar[[i]], tempGrp, Mode, na.rm = TRUE)[[1]], collapse = ", ", sep = "")))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Mode"
		outputLabel <- c(outputLabel,"Mode")
	}
	if(!is.na(match("q1", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stats::quantile, probs = 0.25, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Q1"
		outputLabel <- c(outputLabel,"1st Quartile")
	}
	if(!is.na(match("q3", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stats::quantile, probs = 0.75, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Q3"
		outputLabel <- c(outputLabel,"3rd Quartile")
	}
	if(!is.na(match("range", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, max, na.rm = TRUE) - tapply(tempVar[[i]], tempGrp, min, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Range"
		outputLabel <- c(outputLabel,"Range")
	}
	if(!is.na(match("iqr", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stats::quantile,probs = 0.75, na.rm = TRUE) - tapply(tempVar[[i]], tempGrp, stats::quantile, probs = 0.25, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "IQR"
		outputLabel <- c(outputLabel,"Inter Quartile Range")
	}
	if(!is.na(match("var", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, var, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Variance"
		outputLabel <- c(outputLabel,"Variance")
	}
	if(!is.na(match("sd", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, sd, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "StdDev"
		outputLabel <- c(outputLabel,"Standard Deviation")
	}
	if(!is.na(match("se.mean", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, se, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE_Mean"
		outputLabel <- c(outputLabel,"Std. Error of the Mean")
	}
	if(!is.na(match("cv", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, cv, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "CV"
		outputLabel <- c(outputLabel,"Coefficient of Variation")
	}
	if(!is.na(match("css", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, css, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "CSS"
		outputLabel <- c(outputLabel,"Corrected Sum of Squares")
	}
	if(!is.na(match("ucss", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, ucss, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "UCSS"
		outputLabel <- c(outputLabel,"Uncorrected Sum of Squares")
	}
	if(!is.na(match("skew", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, skewness, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Skewness"
		outputLabel <- c(outputLabel,"Skewness")
	}
	if(!is.na(match("se.skew", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stdskewness, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE_Skew"
		outputLabel <- c(outputLabel,"Std. Error of Skewness")
	}
	if(!is.na(match("kurtosis", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, kurtosis, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Kurtosis"
		outputLabel <- c(outputLabel,"Kurtosis")
	}
	if(!is.na(match("se.kurtosis", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stdkurtosis, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE_Kurtosis"
		outputLabel <- c(outputLabel,"Std. Error of Kurtosis")
	}

	if (is.null(ncol(tempGrp))) { summaryTable[,1] <- names(tempVar)
	} else {
		variable <- c(rep(names(tempVar), each = nrow(as.data.frame.table(table(tempGrp)))))
		summaryTable <- data.frame(variable,summaryTable)
	}
	colnames(summaryTable)[1] <- "Variable"

	options(width = 5000)

	cat("DESCRIPTIVE STATISTICS\n")
	print_tabstat(summaryTable)
	return(invisible(summaryTable))
}##--end of tabstat
NULL









print_tabstat <- function(dataFrame, border = TRUE, digits = NULL) {

  if (!is.data.frame(dataFrame)) { stop("The argument 'dataFrame' should be of class data.frame.") }
  dataChar <- DataAttribute(dataFrame)[,1:2]
  dataChar[,1] <- as.character(dataChar[,1])

  colWidth <- NULL
  onlyDecimal <- NULL
  emptyColumn <- NULL
  for (i in (1:nrow(dataChar))) {
    if (all(is.na(dataFrame[,i])) || all(is.nan(dataFrame[,i])) || all(dataFrame[,i] == Inf) || all(dataFrame[,i] == -Inf)) {
      colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 2)
      emptyColumn <- c(emptyColumn, i)
    } else{
      tempsubdata <- subset(dataFrame[,i], !is.na(dataFrame[,i]))
      if (length(tempsubdata) == 0) {
        colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 2)
      } else {
        if (dataChar[i,2] == "numeric" || dataChar[i,2] == "integer") {
          tempsubdata <- subset(tempsubdata, !is.nan(tempsubdata))
          if (length(tempsubdata) == 0) {
            colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 2)
          } else {
            tempsubdata <- subset(tempsubdata, tempsubdata != Inf)
            if (length(tempsubdata) == 0) {
              colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 2)
            } else {
              tempsubdata <- subset(tempsubdata, tempsubdata != -Inf)
              if (length(tempsubdata) == 0) {
                colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 2)
              } else {
                if (all(tempsubdata < 1) && all(tempsubdata > -1)) {
                  if (is.null(digits)) { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 5) + 2)
                  } else { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + digits) + 2) }
                  #colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 5) + 1)
                  onlyDecimal <- c(onlyDecimal, i)
                } else {
                  if (all(is.wholenumber(tempsubdata))) {
                    dataChar[i,2] <- "integer"
                    colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), nchar(max(tempsubdata))) + 2)
                  } else {
                    if (is.null(digits)) { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 3) + 2)
                    } else { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + digits) + 2) }
                  }
                }
              }
            }
          }
        } else {
          # if variable is string
          colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 2)
        }
      }
    } ## end if-else stmt
  } ## end for stmt

  if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")

  for (j in (1:ncol(dataFrame))) {
    if (dataChar[j,2] == "numeric" || dataChar[j,2] == "integer") {
      cat(formatC(dataChar[j,1], width = colWidth[j], format = "s"), sep = "")
    } else { cat(formatC(dataChar[j,1], width = colWidth[j], format = "s", flag = "-"), sep = "") }
    if (j == ncol(dataFrame)) {
      cat("\n")
      #if (nrow(dataFrame) >= 2) cat("\n\n") else cat("\n")
    } else { cat(formatC("", width = 1, format = "s")) }
  }

  if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")

  for (i in (1:nrow(dataFrame))) {
    for (j in (1:ncol(dataFrame))) {
      if(!is.na(match(j,emptyColumn))) { cat(formatC("", width = colWidth[j], format = "s", flag = "-"), sep = "")
      } else {
        if (dataChar[j,2] == "numeric") {
          if (is.na(match(j, onlyDecimal))) {
            if (is.na(dataFrame[i,j]) || is.nan(dataFrame[i,j]) || dataFrame[i,j] == Inf || all(dataFrame[i,j] == -Inf)) { cat(formatC("", width = colWidth[j], format = "s"), sep = "")
            } else {
              if (!is.null(digits)) { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = digits), sep = "")
              } else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = 2), sep = "") }
            }
          } else {
            if (is.na(dataFrame[i,j]) || is.nan(dataFrame[i,j]) || dataFrame[i,j] == Inf || all(dataFrame[i,j] == -Inf)) { cat(formatC("", width = colWidth[j], format = "s"), sep = "")
            } else {
              if (!is.null(digits)) { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = digits), sep = "")
              } else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = 4), sep = "") }
            }
          }
        } else {
          if (dataChar[j,2] == "integer") {
            if (is.na(dataFrame[i,j]) || is.nan(dataFrame[i,j]) || dataFrame[i,j] == Inf || all(dataFrame[i,j] == -Inf)) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") }
            else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "d"), sep = "") }
          } else {
            if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s", flag = "-"), sep = "") }
            else { cat(formatC(as.character(dataFrame[i,j]), width = colWidth[j], format = "s", flag = "-"), sep = "") }
          }
        }
      }
      if (j == ncol(dataFrame)) { cat("\n") } else { cat(formatC("", width = 1, format = "s")) }
    }
  }
  if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")

}
NULL




DataAttribute <- function(.data) {
  if(is.character(.data)) {
    nameData <- .data
    if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
    tempData <- eval(parse(text = .data))
  } else {
    nameData <- paste(deparse(substitute(.data)))
    #if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
    tempData <- .data
  }
  if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
  tempTable <- data.frame(names(tempData))
  withFactor <- FALSE
  for (i in (1:nrow(tempTable))) {
    if (is.factor(tempData[,i])) {
      withFactor <- TRUE
      if (is.ordered(tempData[,i])) tempTable[i,2] <- "ordered factor"
      else tempTable[i,2] <- "factor"
      tempTable[i,3] <- as.character(nlevels(tempData[,i]))
      if (nlevels(tempData[,i]) > 5) tempTable[i,4] <- paste(c(paste(levels(tempData[,i])[1:2], collapse = ", ", sep = ""),levels(tempData[,i])[nlevels(tempData[,i])]), collapse = ", ..., ", sep = "")
      else tempTable[i,4] <- paste(levels(tempData[,i]), collapse = ", ", sep = "")
    } else {
      if (typeof(tempData[,i]) == "double") { tempTable[i,2] <- "numeric" }
      else tempTable[i,2] <- typeof(tempData[,i])
      tempTable[i,3] <- ""
      tempTable[i,4] <- ""
    }
  }
  names(tempTable) <- c("VAR NAME", "TYPE", "NLEVELS", "LEVELS")
  if (!withFactor) { tempTable <- tempTable[,1:2] }
  return(tempTable)
}
NULL
