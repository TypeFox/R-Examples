################################################################################
### toLatex-method for "sts" objects
###
### Copyright (C) 2014 Dirk Schumacher, 2014 Maelle Salmon
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################


toLatex.sts <- function(object, caption = "",label=" ", columnLabels = NULL,
                        subset = NULL, 
                        alarmPrefix = "\\textbf{\\textcolor{red}{",
                        alarmSuffix = "}}", ubColumnLabel = "UB", ...) {
  # Args:
  #   object: A single sts object; must not be NULL or empty.
  #   caption: A caption for the table. Default is the empty string.
  #   label: A label for the table. Default is the empty string.
  #   columnLabels: A list of labels for each column of the resulting table.
  #   subset: A range of values which should be displayed. If Null, then all 
  #           data in the sts objects will be displayed. Else only a subset of 
  #           data. Therefore range needs to be a numerical vector of indexes
  #           from 1 to length(@observed).
  #   alarmPrefix: A latex compatible prefix string wrapped around a table cell 
  #                iff there is an alarm;i.e. alarm = TRUE
  #   alarmSuffix: A latex compatible suffix string wrapped around a table cell 
  #                iff there is an alarm;i.e. alarm[i,j] = TRUE
  #   ubColumnLabel: The label of the upper bound column; default is "UB".
  #   ...: Variable arguments passed to toLatex.xtable
  # Returns:
  #  An object of class Latex
  
  # Error Handling
  isEmpty <- function(o) is.null(o)
  if (isEmpty(object))
    stop("object must not be null or NA.")
  
  if (is.list(object))
    stop("supplying a list of sts has been removed from the api. Sorry.")
  
  if (!isS4(object) || !is(object, "sts"))
    stop("object must be of type sts from the surveillance package.")
  
  if (!is.character(caption))
    stop("caption must be a character.")
  
  if (!isEmpty(labels) && length(labels) != length(object))
    stop("number of labels differ from the number of sts objects.")
    
  # derive default values
  
  tableLabels <- colnames(object@observed)
  if (!is.null(columnLabels) && 
        length(columnLabels) != ncol(object@observed) * 2 + 2) {
    stop("the number of labels must match the number of columns in the 
         resulting table; i.e. 2 * columns of sts + 2.")
  }
  
  tableCaption <- caption
  tableLabel <- label
  
  vectorOfDates <- epoch(object, as.Date = TRUE)
  
  yearColumn <- Map(function(d)isoWeekYear(d)$ISOYear, vectorOfDates)
  
  if (object@freq == 12 )
    monthColumn <- Map(function(d) as.POSIXlt(d)$mon, vectorOfDates)
  
  if (object@freq == 52 )
    weekColumn <- Map(function(d)isoWeekYear(d)$ISOWeek, vectorOfDates)
  
  dataTable <- data.frame(unlist(yearColumn))
  colnames(dataTable) <- "year"
  
  if (object@freq == 12 ) {
    dataTable$month <- unlist(monthColumn)  
  }
  
  if (object@freq == 52 ) {
    dataTable$week <- unlist(weekColumn)
  }
  
  if (object@freq == 365 ) {
    dataTable$day <- unlist(vectorOfDates)
    dataTable <- dataTable[c(2)]
  }
  
  noCols <- ncol(dataTable)
  j <- 1 + noCols
  tableLabelsWithUB <- c()
  
  # I know it is imperative - shame on me
  for (k in 1:(ncol(object@observed))) {
    upperbounds <- round(object@upperbound[,k], 2)
    observedValues <- object@observed[,k]
    alarms <- object@alarm[,k]
    ubCol <- c()
    for (l in 1:length(upperbounds)) {
        if (is.na(upperbounds[l])) {
            ubCol <- c(ubCol, NA)
        } else {
            ubCol <- c(ubCol, upperbounds[l])
            if (!is.na(alarms[l]) && alarms[l]) {
                observedValues[l] <- paste0(alarmPrefix, observedValues[l], alarmSuffix)
            }
        }
    }
    dataTable[,(j)] <- observedValues
    dataTable[,(j + 1)] <- ubCol
    tableLabelsWithUB <- c(tableLabelsWithUB, tableLabels[k])
    tableLabelsWithUB <- c(tableLabelsWithUB, ubColumnLabel)
    j <- j + 2
  }
  
  # remove rows which should not be displayed
  if (is.null(subset))
    subset <- 1:nrow(dataTable)
  else if (min(subset) < 1 || max(subset) > nrow(dataTable))
    stop("'subset' must be a subset of 1:nrow(observed), i.e., 1:",
         nrow(dataTable))
  
  dataTable <- dataTable[subset,]
  
  # prepare everything for xtable
  newColNames <- c(colnames(dataTable)[1:noCols], tableLabelsWithUB)
  if (!is.null(columnLabels)) {
    colnames(dataTable) <- columnLabels
  } else {
    colnames(dataTable) <- newColNames
  }
  xDataTable <- xtable(dataTable, label = tableLabel, caption = tableCaption, digits = c(0))
  toLatex(xDataTable, ...) 
}

setMethod("toLatex", "sts", toLatex.sts)
