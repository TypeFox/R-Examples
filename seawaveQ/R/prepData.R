#' Prepares concentration data and continuous ancillary data
#' 
#' prepData is usually called from within \link{fitswavecav} but
#' can be invoked directly.  It performs some date calculations, removes
#' rows with missing values for concentration or continous variables, 
#' and returns the the concentration and continuous ancillary data to
#' be used by \link{fitswavecav} and its other internal functions.
#' @param cdat is the concentration data.
#' @param cavdat is the continuous (daily) ancillary data.
#' @param yrstart is the starting year of the analysis (treated as January
#' 1 of that year).  Zero means the start date will be determined by the 
#' start date of cavdat, the continuous ancillary data.
#' @param yrend is the ending year of the analysis (treated as December 31
#' of that year).  Zero means the end date will be determined by the end 
#' date of cavdat, the continuous ancillary data.
#' @param dcol is the column name for the dates, should be the same for 
#' both cdat and cavdat.
#' @param pnames are the parameters (water-quality constituents) to 
#' analyze (if using USGS parameters, omit the the starting 'P', such as 
#' "00945" for sulfate).  
#' @param iwcav is a character variable indicating which continuous
#' ancillary variables to include, if none use iwcav=c("none").
#' @param qwcols is a character vector with the beginning of the
#' column headers for remarks code (default is R), and beginning of 
#' column headers for concentration data (default is P for parameter).
#' @keywords manip
#' @return a list.  The first element is the concentration data with
#' additional date information, missing values removed, and extra columns
#' removed.  The second element is the continuous ancillary data with
#' additional date information, missing values removed, and extra columns
#' removed.  
#' @export
#' @author Aldo V. Vecchia and Karen R. Ryberg
#' @examples
#' data(swData)
#' modMoRivOmaha<-combineData(qwdat=qwMoRivOmaha, cqwdat=cqwMoRivOmaha)
#' preppedDat <- prepData(modMoRivOmaha, cqwMoRivOmaha, yrstart=1995, 
#' yrend=2003, dcol="dates", pnames=c("04035", "04037", "04041"),  
#' iwcav=c("flowa30","flowa1"), qwcols=c("R","P"))
prepData <- function(cdat, cavdat, yrstart, yrend, dcol, pnames, 
                     iwcav, qwcols) {
  cdat <- subset(cdat, year(cdat[,dcol]) >= yrstart & 
                   year(cdat[, dcol]) <= yrend)
  cdat$yrc <- year(cdat[, dcol])
  cdat$moc <- month(cdat[, dcol])
  cdat$dac <- day(cdat[, dcol])
  cdat$jdayc <- julian(cdat[, dcol], 
                       origin=as.Date(paste(yrstart-1, 
                                            "-10-01", sep="")))
  # column headings of concentration data
  mycols<-c("yrc","moc","dac","jdayc")
  rnames <- paste(qwcols[1], pnames, sep='')
  pnames <- paste(qwcols[2], pnames, sep='')
  
  for ( i in 1:length(pnames) ) {
    mycols<-c(mycols, rnames[i], pnames[i])
  }
  if ( iwcav[1]!= "none" ) {
    mycols <- c(mycols, iwcav)
  }
  cdat <- cdat[, mycols]
  # columns of cdat are year, month, day, julian day, remark code, 
  # concentration value, and the selected ancillary variables
  
  if ( iwcav[1]!= "none" ) {
    pcktmp <- !is.na(cdat[,4])
    for (j in 1:length(iwcav) ) {
      pcktmp <- pcktmp & !is.na(cdat[,iwcav[j]])
    }
    cdat <- cdat[pcktmp,]
  }
  
  cavdat <- subset(cavdat, year(cavdat[,dcol]) >= yrstart & 
                     year(cavdat[, dcol]) <= yrend)
  cavdat$yrx <- year(cavdat[, dcol])
  cavdat$mox <- month(cavdat[, dcol])
  cavdat$dax <- day(cavdat[, dcol])
  cavdat$jdayx <- julian(cavdat[, dcol], 
                         origin=as.Date(paste(yrstart-1, 
                                              "-10-01", sep="")))
  # determine column headings of ancillary data
  mycols<-c("yrx","mox","dax","jdayx")
  
  if ( iwcav[1]!= "none" ) {
    mycols <- c(mycols, iwcav)
  }
  cavdat <- cavdat[, mycols]
  # columns of cavdat are year, month, day, julian day, and the 
  # selected continuous variables 
  # remove rows with missing values 
  # for concentration or continuous variables
  
  if ( iwcav[1]!= "none" ) {
    pcktmp <- !is.na(cavdat[,5])
    if(length(cavdat[1,]) > 5)  {
      for (j in 6:length(cavdat[1,])) {
        pcktmp <- pcktmp & !is.na(cavdat[,j])
      }
    }
    cavdat <- cavdat[pcktmp,]
  }
  myData<-list(cdat, cavdat)
  myData
}  
