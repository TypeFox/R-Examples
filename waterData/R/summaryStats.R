#' Function to identify and fix common problems with hydrologic data
#' 
#' @name cleanUp
#' @title Cleans up hydrologic time series data
#' @param dataset is a data frame in format of the data frame returned by 
#' \link{importDVs}.
#' @param task is either "view" or "fix."  View will return a list containing 
#' rows with negative values and rows with missing values.  Fix will 
#' replace negative values with NA and replace zeroes with the value
#' specified by the replace argument.  
#' @param replace is the value used to replace 0 values.  The default
#' is 0.1.  For streamflow in small streams, one might want to use 0.01.
#' For daily data other than streamflow, such as turbidity, users may not want
#' to replace 0 values with a nonzero value.  In those cases, replace can be 
#' set to 0.
#' @note If calculating anomalies (see \link{compAnom}), the user may need to 
#' replace isolated missing values with with a value; however, if there are 
#' larger periods with missing values, streamflow anomalies may not be an 
#' appropriate use of the data.  The substitution of some missing data with 
#' values may be done using the function \link{fillMiss} that is part of this 
#' package.  However, care needs to be taken when filling in missing data.
#' @keywords NA ts utilities manip
#' @return A list showing days with negative values and days with 
#' values of 0 when task is "view."  When task is "fix" the fixed dataset 
#' is returned.
#' When a negative value is replaced with NA, an "N" is added to the qualcode
#' field to indicate that there had been a negative number.
#' When a zero value is replaced, an "R" is added to the qualcode field to
#' indicate that a zero value was replaced.
#' @seealso \link{fillMiss}
#' @export
#' @examples
#' data(exampleWaterData)
#' head(badDataSet)
#' cleanUp(badDataSet, task="view")
#' q05054000Fix <- cleanUp(badDataSet, task="fix")
#' # replace 0s with NA, then one could use the fillMiss function
#' # to estimate values
#' q05054000Fix2 <- cleanUp(badDataSet, task="fix", replace=NA)
#' summary(badDataSet)
#' summary(q05054000Fix)
#' summary(q05054000Fix2)
cleanUp <- function(dataset, task="view", replace=0.1) {
  if (replace < 0 & !is.na(replace) ) {
    stop("The value of replace must be greater than or equal to 0.")
  }
  if (replace > 10 & !is.na(replace) ) {
    stop("The value of replace must be less than 10.")
  }
  pck <- dataset$val < 0 & !is.na(dataset$val)
  # dataset[pck,]
  pck2 <- dataset$val == 0 & !is.na(dataset$val)
  if (task == "view") {
    list(dataset[pck,], dataset[pck2,])
  }
  else if (task == "fix") {
    dataset$qualcode<-as.character(dataset$qualcode)
    dataset$val[pck] <- NA
    # add N to qualcode to indicate that there had been a negative number
    dataset$qualcode[pck] <- paste(dataset$qualcode[pck], "N",sep=" ")
    dataset$val[pck2] <- replace
    # add R to qualcode to indicate 0s were replaced
    dataset$qualcode[pck2] <- paste(dataset$qualcode[pck2], "R",sep=" ")
    dataset$qualcode<-factor(dataset$qualcode)
    dataset
  }
  else {
    stop("Task must be view or fix.")
  }
}

#' Function to calculate summary statistics for daily hydrologic time series.
#'
#' The summary statistics returned are useful for exploratory data analysis 
#' and for describing the date set.
#' @note Hydrologic data are often skewed (Helsel and Hirsch, 2002).  Summary 
#' statistics help describe the degree of skewness and help to determine
#' the degree of applicability of hypothesis tests.  Some data, in particular
#' streamflow, may need to be transformed to produce approximately normal
#' data.
#' @name summaryStats
#' @title Calculate summary statistics
#' @param dataset is the data frame containing hydrologic data
#' @param staid is used to label the output 
#' @keywords arith
#' @return a data frame containing a number of summary statistics of the daily 
#' hydrologic data series
#' @export
#' @format The returned matrix has the following columns, which are formatted
#' for putting in a report or table. \cr
#' \tabular{lll}{
#' Name \tab Type \tab Description \cr 
#' Begin \tab character \tab The beginning date of the time series \cr
#' End \tab character \tab The ending date of the time series \cr
#' n \tab character \tab Number of rows \cr
#' NA \tab character \tab Number of missing values \cr
#' Neg \tab character \tab Number of negative values \cr
#' Min \tab character \tab The minimum value \cr
#' Q1 \tab character \tab The first quartile, 25th percentile \cr
#' Med \tab character \tab The median \cr
#' Mean \tab character \tab The mean \cr
#' Q3 \tab character \tab The third quartile, 75th percentile \cr
#' Max \tab character \tab The maximum value \cr
#' StdDev \tab character \tab The standard deviation \cr
#' IQR \tab character \tab The interquartile range \cr
#'}
#' @examples 
#' data(exampleWaterData)
#' summaryStats(pH05082500, staid="05082500")
#' @references
#' Helsel, D.R. and Hirsch, R. M., 2002, Statistical methods in water resources: 
#' U.S. Geolgical Survey Techniques of Water Resources Investigations, book 4, 
#' chap. A3, 522 p. (Also available at \url{http://pubs.usgs.gov/twri/twri4a3/}).
summaryStats<-function(dataset,staid=1) {
  sdate<-dataset$dates[1]
  edate<-dataset$dates[length(dataset$dates)]
  n<-length(dataset$val)
  pck<-is.na(dataset$val)
  missing<-sum(pck)
  pck<-dataset$val<0&!is.na(dataset$val)
  negative<-sum(pck)
  my.sum<-fivenum(dataset$val,na.rm=TRUE)
  my.min<-my.sum[1]
  my.25<-my.sum[2]
  my.med<-my.sum[3]
  qmean<-mean(dataset$val,na.rm=TRUE)
  my.75<-my.sum[4]
  my.max<-my.sum[5]
  my.sd<-sd(dataset$val,na.rm=TRUE)
  my.iqr<-IQR(dataset$val,na.rm=TRUE)
  my.dfnums<-cbind(missing, negative, my.min, my.25, my.med, qmean, my.75,
                   my.max, my.sd, my.iqr)
  n<-format(n, digits=1, big.mark=",", scientific=FALSE)
  my.dfnums<-format(my.dfnums, digits=1, big.mark=",", scientific=FALSE)
  my.df<-as.data.frame(cbind(as.character(sdate),as.character(edate),n,my.dfnums),
                       stringsAsFactors=FALSE)
  dimnames(my.df)[[2]]<-c("Begin", "End", "n", "NA", "Neg", "Min","Q1", "Med",
                          "Mean", "Q3", "Max", "StdDev", "IQR")

  row.names(my.df)<-staid
  my.df
}
