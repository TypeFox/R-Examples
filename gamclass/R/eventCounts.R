eventCounts <-
function (data, dateCol="Date", from = NULL, to = NULL,
          by = "1 month", categoryCol=NULL, takeOnly=NULL, prefix="n_")
{
  checkCols <- c(dateCol, categoryCol) %in% names(data)
  if(!is.null(categoryCol) & !all(checkCols)){
    txt <- paste("Name(s)", c(dateCol, categoryCol)[!checkCols], "not found in", deparse(data))
    stop(txt)
  }
  if(!is.null(takeOnly)){
    subdat <- eval(parse(text=takeOnly), data)
    data <- subset(data, subdat)
  }
  date <- data[, dateCol]
  if(!is(date, "Date")){date <- try(as.Date(date), silent=TRUE)
                     if(class(date)=="try-error")
                       stop(paste("Column", dateCol, "must hold a date object"))
  }
    if (is.null(from))
        from <- min(date)
    if (is.null(to))
        to <- max(date)
    dateBreaks <- seq(from = from, to = to, by = by)
    dateBreaks <- c(dateBreaks, max(dateBreaks) + diff(dateBreaks[1:2]))
    countDF <- data.frame(Date = dateBreaks[-length(dateBreaks)])
    if(!is.null(categoryCol))categs <- names(table(data[, categoryCol])) else categs <- ""
    for(cat in categs){
      if(!is.null(categoryCol)) select <- data[, categoryCol] == cat else
        select <- rep(TRUE, nrow(countDF))
      cutDates <- cut(date[select], dateBreaks, right = FALSE)
      countNam <- paste0(prefix, gsub(" ", "", cat))
      countDF[, countNam] <- as.vector(table(cutDates))
    }
    countDF
}
