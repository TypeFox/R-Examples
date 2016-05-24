#' Generate a classical TXT point and figure plot.
#' 
#' THIS FUNCTION IS STILL UNDER DEVELOPMENT,
#' THEREFORE IT MIGHT BE SUBJECT TO CHANGE!
#' 
#' @param data a data frame object containing point and figure informations to be plotted
#' @param reversal number of boxes used in pnfprocessor
#' @param boxsize the boxsize used in pnfprocessor
#' @param log are calculations done in logarithmic mode
#' @param main a string used as a main title of the chart
#' @param sub a string used as a sub title of the chart
#' @seealso \code{\link{pnfprocessor}}
#' @seealso \code{\link{pnfplot}}
#' @references \url{http://rpnf.r-forge.r-project.org}
#' @export
#' @examples
#' library(rpnf) # Load rpnf library
#' data(DOW) # (Offline) Load free available sample data from https://www.quandl.com/data/WIKI/DOW
#' pnfdata <- pnfprocessor(
#'   high=DOW$High,
#'   low=DOW$Low,
#'   date=DOW$Date,
#'   boxsize=1L,
#'   log=FALSE)  
#' pnfplottxt(pnfdata,boxsize=1L,log=FALSE)
pnfplottxt <- function(data,reversal=3,boxsize=1,log=FALSE,main=NULL,sub=NULL) {
  # Check if warning level is set to zero, otherwise output will be messy
  old.waring.level <- getOption("warn")
  if (old.waring.level>0) {
    warning("Warning level > 0. This would produce a messy plot. Reset option(warn)=0 temporally!")
    options(warn=0)
  }
  
  ## local function definiton: plot seperation line
  plotSeperationLine <- function(numOfColumns) {
    # iterate through columns
    cat("--------+")
    for (column in 1:numOfColumns) 
      cat("-")    
    cat("\n")
  }
  ## Write main and sub title
  if (!is.null(main))
    cat(paste(main,"\n"))
  if (!is.null(sub))
    cat(paste(sub,"\n"))
  # plot seperation line
  plotSeperationLine(max(data$column)-min(data$column)+1)
  # cat to connection object
  for (mybox in max(data$boxnumber):min(data$boxnumber)) {
    ### iterate over every line
    cat(format(round(box2lower(mybox,boxsize=boxsize,log=log),2),width=8,nsmall=2))
    cat("|")
    # iterate through columns
    for (column in min(data$column):max(data$column)) {
      status <- as.character(unique(data$status.xo[data$column==column]))
      mymin <- min(data$boxnumber[data$column==column])
      mymax <- max(data$boxnumber[data$column==column])
      # correct mymin and mymax, if necessary
      if (column>min(data$column)) {
        if (status=="X") 
          mymin <- min(data$boxnumber[data$column==column-1])+1
        else
          mymax <- max(data$boxnumber[data$column==column-1])-1
      }
      # determine trendline boxes, if available
      trendline <- NA
      if ("tl.brl.boxnumber" %in% names(data) | "tl.bsl.boxnumber" %in% names(data)) {
#        trendline <- c(min(data[data$column==column, c("tl.brl.boxnumber","tl.bsl.boxnumber")],na.rm=T),
#                       max(data[data$column==column, c("tl.brl.boxnumber","tl.bsl.boxnumber")],na.rm=T))
        # FIXME this is a bottleneck
        trendline <- unique(data[data$column==column, c("tl.brl.boxnumber","tl.bsl.boxnumber")])
      }
      # FIXME this is a bottleneck
      # decide on plot
      if (mymin<=mybox & mybox<=mymax) {
        cat(status)
      } else if (mybox %in% trendline) {
        cat("+")
      } else
        cat(" ")
    } # end column loop
    # check, if current line is the latest quote
    if (data$boxnumber[nrow(data)]==mybox) {
      # write a marker for current box on rhs
      cat(" <==")
    } else {
      # write a marker for current box on rhs
      cat("    ")
    }
    # right hand side 
    # cat("|")
    #cat(format(round(box2lower(mybox,boxsize=boxsize,log=log),2),width=8,nsmall=2))
    # write line feed
    cat("\n")
  }
  # plot seperation line
  plotSeperationLine(max(data$column)-min(data$column)+1)
  ## write date lines, but vertical
  for (pos in 1:10) {
    if (pos %in% c(1,2,3,4))
      cat("       Y|")
    else if (pos %in% c(6,7))
      cat("       M|")
    else if (pos %in% c(9,10))
      cat("       D|")
    else
      cat("        |")
    for (column in min(data$column):max(data$column)) {
      if (pos == 5 | pos == 8)
        cat(" ")
      else
        cat(substr(as.character(min(data$date[data$column==column])),start=pos,stop=pos))
    }
    cat("\n")
  }

  # restore default options
  if (old.waring.level>0) {
    warning("Restoring old warning level!")
    options(warn=old.waring.level)
  }
}
