#' Creates data 
#'
#' @description R implementation of the SPSS \code{COMPUTE} argument.
#' 
#' @usage xpssCompute(x, variables = NULL, fun = NULL,...)
#' @param x input data .
#' @param variables atomic character or character vector with the names of the variables.
#' @param fun atomic character as functionname.
#' @param ... further arguments passed to or from other methods.
#' @details
#' \strong{Missing Fuctions:} 
#' \tabular{rll}{
#' \tab \code{Functionname} \tab \code{Output}
#' \cr \tab \code{\link{computeMiss}} \tab Computes missing values as logical (limited to one variable).
#' \cr \tab \code{\link{computeNmiss}} \tab Computes missing values as logical.
#' \cr \tab \code{\link{computeNvalid}} \tab Computes valid values as logical.
#' \cr \tab \code{\link{computeSysmis}} \tab Computes the system missing values as logical.
#' \cr \tab \code{\link{computeValue}} \tab Computes the user-defined values as valid values.
#' }
#' \strong{Numeric Fuctions:} 
#' \tabular{rll}{
#' \tab \code{Functionname} \tab \code{Output}
#' \cr \tab \code{\link{computeAbs}} \tab Computes the absolute value. 
#' \cr \tab \code{\link{computeArsin}} \tab Computes the arc-sine.
#' \cr \tab \code{\link{computeArtan}} \tab Computes the arc-tan.
#' \cr \tab \code{\link{computeCos}} \tab Computes the cosinus.
#' \cr \tab \code{\link{computeExp}} \tab Computes the exponential.
#' \cr \tab \code{\link{computeLn}} \tab Computes the logarithmus naturalis.
#' \cr \tab \code{\link{computeLg10}} \tab Computes the logarithmus base 10.
#' \cr \tab \code{\link{computeLngamma}} \tab Computes the logarithmus of the gamma function.
#' \cr \tab \code{\link{computeMax}} \tab Computes the maxima.
#' \cr \tab \code{\link{computeMean}} \tab Computes the atithmetic mean.
#' \cr \tab \code{\link{computeMedian}} \tab Computes the median value.
#' \cr \tab \code{\link{computeMin}} \tab Computes the minima.
#' \cr \tab \code{\link{computeMod}} \tab Computes the remainder of a division.
#' \cr \tab \code{\link{computeRnd}} \tab Computes rounded values.
#' \cr \tab \code{\link{computeSd}} \tab Computes the standard deviation.
#' \cr \tab \code{\link{computeVariance}} \tab Computes the variance.
#' }
#' 
#' \strong{Character Fuctions:} 
#' \tabular{rll}{
#' \tab \code{Functionname} \tab \code{Output}
#' \cr \tab \code{\link{computeChar_index}} \tab Returns the position of the first occurence of a pattern.
#' \cr \tab \code{\link{computeChar_length}} \tab Computes the length of a string in characters.
#' \cr \tab \code{\link{computeChar_lpad}} \tab Returns an expanded strings.
#' \cr \tab \code{\link{computeChar_mblen}} \tab Computes the byte per character or sign.
#' \cr \tab \code{\link{computeChar_rindex}} \tab Returns the position of the last occurence of a pattern.
#' \cr \tab \code{\link{computeChar_rpad}} \tab Returns an expanded strings.
#' \cr \tab \code{\link{computeConcat}} \tab  Returns a concatenated string.
#' \cr \tab \code{\link{computeLength}} \tab Computes Number of bytes in a string.
#' \cr \tab \code{\link{computeLower}} \tab Returns the input data to lower-case.
#' \cr \tab \code{\link{computeLtrim}} \tab Returns a trimmed string (left side trimmed).
#' \cr \tab \code{\link{computeNtrim}} \tab Returns the values of the input data.
#' \cr \tab \code{\link{computeReplace}} \tab Replaces a pattern in a string.
#' \cr \tab \code{\link{computeRtrim}} \tab Returns a trimmed string (right side trimmed).
#' \cr \tab \code{\link{computeStrunc}} \tab Returns a truncated string.
#' \cr \tab \code{\link{computeUpcase}} \tab Returns the input data to upper-case.
#' }
#' \strong{Lag Fuction:} 
#' \tabular{rll}{
#' \tab \code{Functionname} \tab \code{Output}
#' \cr \tab \code{\link{computeLag}} \tab Shifts the variable backward or forward.
#' }
#' \strong{Date Fuctions:} 
#' \tabular{rll}{
#' \tab \code{Functionname} \tab \code{Output}
#' \cr \tab \code{\link{computeCtime_days}} \tab Computes the difference of time between two dates in days.
#' \cr \tab \code{\link{computeCtime_hours}} \tab Computes the difference of time between two dates in hours.
#' \cr \tab \code{\link{computeCtime_minutes}} \tab Computes the difference of time between two dates in minutes.
#' \cr \tab \code{\link{computeCtime_seconds}} \tab Computes the difference of time between two dates in seconds.
#' \cr \tab \code{\link{computeDate_dmy}} \tab Computes a date with the format day-month-year.
#' \cr \tab \code{\link{computeDate_mdy}} \tab  Computes a date with the format month-day-year.
#' \cr \tab \code{\link{computeDate_moyr}} \tab Computes a date with the format month-year.
#' \cr \tab \code{\link{computeDate_qyr}} \tab Computes a date with the format year/quarter.
#' \cr \tab \code{\link{computeDate_wkyr}} \tab Computes a date with the format year/calendar week.
#' \cr \tab \code{\link{computeDate_yrday}} \tab Computes a date with the format year/yearday.
#' \cr \tab \code{\link{computeTime_days}} \tab Computes the number of passed hours on basis of the given days
#' \cr \tab \code{\link{computeTime_hms}} \tab Computes a date with the format Hour-Minute-Second
#' 
#' \cr \tab \code{\link{computeXdate_date}} \tab Extracts the date out of a date string.
#' \cr \tab \code{\link{computeXdate_hour}} \tab Extracts the hour value out of a given date.
#' \cr \tab \code{\link{computeXdate_jday}} \tab Computes the date of year on basis of a given date.
#' \cr \tab \code{\link{computeXdate_mday}} \tab Extracts the date of month on basis of a given date.
#' \cr \tab \code{\link{computeXdate_minute}} \tab Extracts the minute component on basis of a given date.
#' \cr \tab \code{\link{computeXdate_month}} \tab Extracts the month component on basis of a given date.
#' \cr \tab \code{\link{computeXdate_quarter}} \tab Computes the quarter on basis of a given date.
#' \cr \tab \code{\link{computeXdate_second}} \tab Extracts the second on basis of a given date.
#' \cr \tab \code{\link{computeXdate_tday}} \tab Computes the difference of days between the entered date and October 14, 1582.
#' \cr \tab \code{\link{computeXdate_time}} \tab Extracts the time componente on basis of a given date.
#' \cr \tab \code{\link{computeXdate_week}} \tab Calcualtes the calendar week on basis of a given date.
#' \cr \tab \code{\link{computeXdate_wkday}} \tab Calcualtes the day of week on basis of a given date.
#' \cr \tab \code{\link{computeXdate_year}} \tab Extracts the year on basis of a given date.
#' }
#' 
#' Be careful about the input format of the numeric values in your data. It is possible to specify values which are outside the of valid range. \cr Those failures are called \code{Domain Errors}.
#' \cr \cr \strong{For example}:
#' \tabular{rll}{
#' \tab \code{Function} \tab \code{Output} 
#' \cr \tab **  \tab A negative number to a noninteger power.
#' \cr \tab /  \tab A divisor of 0.
#' \cr \tab computeArsin \tab  An argument whose absolute value exceeds 1.
#' \cr \tab computeExp \tab An argument that produces a result too large to be represented on the computer.
#' \cr \tab computeLg10 \tab  A negative or 0 argument.
#' \cr \tab computeLn  \tab A negative or 0 argument.
#' \cr \tab computeMod \tab A divisor of 0.
#' \cr \tab computeSqrt \tab  A negative argument.
#' }
#' @return Output is a created vector.
#' @author Bastian Wiessner
#' @seealso Related Functions \code{\link{xpssNumeric}} \code{\link{xpssString}} 
#' @examples 
#' data(fromXPSS)
#' 
#' xpssCompute(x=fromXPSS, variables="V7_2",fun="computeValue")
#' 
#' xpssCompute(x=fromXPSS, variables=c("V5","V7_2"),fun="computeMean", na.rm=T)
#' @export
xpssCompute <- function(x, variables = NULL, fun = NULL,...){
  
  if("computeMax" %in% fun || "computeMean" %in% fun || "computeMedian" %in% fun || "computeMin" %in% fun || "computeSd" %in% fun || "computeVariance" %in% fun){
    if(length(variables)<2){
      stop("For the following operations are at least two variables needed: computeMax, computeMean, computeMedian, computeMin, computeSd, computeVariance")
    }
  }
  
  #unique date statements
  value <- c("Ctime","Date","Xdate")
  dots <-  as.list(substitute(list(...)))[-1L]
  if(is.null(dots)){
    eval(call(fun,x))
  }else{
    if(is.null(x)){
      stop("argument x is missing")
    } else {
      
      #
      #
      #
      #
      ##
      ####
      if(is.atomic(x)){
        data <- list(x)
      }  
      
      #
      #
      #
      #
      ##
      ####
      if(is.data.frame(x) | is.xpssFrame(x)) {
        if(is.null(variables)){
          stop("variables argument is missing")
        } else{
          for(i in 1:length(variables)){
            if(!(is.element(variables[i],names(x)))) {
              stop("The selected variable has to be in the dataset")
            }
          } 
          if("xpssFrame" %in% class(x)) {
            attBack <- attributesBackup(x)
            data <- subset(x, select = variables)
            temp <- paste0("attBack$local$",variables)
            
            backupatts <- attBack[[4]][1]
            for(i in 1:length(temp)){
              backupatts[[i]] <- eval(parse(text=temp[[i]]))    
            }
          } else{
            data <- subset(x, select = variables)
          }
          data <-list(as.matrix(data))
        }        
      }
    }
    
    #
    #
    #
    #
    ##
    ####
    if(is.matrix(x)){
      if(is.null(variables)){
        stop("variables argument is missing")
      } else{
        for(i in 1:length(variables)){
          if(!(is.element(variables[i],colnames(x)))) {
            stop("The selected variable has to be in the dataset")
          }
        } 
        data <- subset(x, select = variables)
      }
      data <-list(data)
    }    
    
    ###############################################################
    
    
    data <-c(data,dots)
    if(length(dots$date)>0){
      if(length(grep(pattern = "date",x = dots))>0){
        data$date <- NULL
        data$date <- get(paste(dots))
      } 
      if(length(grep(pattern = "date",x = dots))<0){
        data$date <- NULL
        data$date <- eval(dots$date)
      }
    }   
    do.call(fun,data)
  }
}
