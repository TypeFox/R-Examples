######################################################################
# Extract ISO week from Date object
#
# Details:
# Code by Gustaf Rydevik <gustaf.rydevik_at_gmail.com> , revised 2010
# http://tolstoy.newcastle.edu.au/R/e10/help/10/05/5588.html
# This is a platform independent way of doing
#  format.Date(x,"%G") or format.Date(x,"%G")
# which unfortunately does not work on windows platforms.
#
# Note: The function is vectorized.
#
# Parameters:
#  Y -- Inputs a date object (POSIX) or Year
#  M -- month (NULL if Y is a Date object)
#  D -- day (NULL if Y is a Date object)
#
# Returns:
#  ISO year and wek of the date
######################################################################


isoWeekYear<-function(Y,M=NULL,D=NULL){
  #Format the date. But whatts the difference between the two statements?
  if(!class(Y)[1]%in%c("Date","POSIXt")) {
    date.posix<-strptime(paste(Y,M,D,sep="-"),"%Y-%m-%d")
  } 
  if(class(Y)[1]%in%c("POSIXt","Date")){
    date.posix<-as.POSIXlt(Y)

    Y<-as.numeric(format(date.posix,"%Y"))
    M<-as.numeric(format(date.posix,"%m"))
    D<-as.numeric(format(date.posix,"%d"))
  }

  #LY
  LY      <- (Y%%4==0 & !(Y%%100==0))|(Y%%400==0)
  LY.prev <- ((Y-1)%%4==0 & !((Y-1)%%100==0))|((Y-1)%%400==0)

  date.yday<-date.posix$yday+1
  jan1.wday<-strptime(paste(Y,"01-01",sep="-"),"%Y-%m-%d")$wday
  jan1.wday<-ifelse(jan1.wday==0,7,jan1.wday)
  date.wday<-date.posix$wday
  date.wday<-ifelse(date.wday==0,7,date.wday)


  ####If the date is in the beginning, or end of the year,
  ### does it fall into a week of the previous or next year?
  Yn<-ifelse(date.yday<=(8-jan1.wday)&jan1.wday>4,Y-1,
        ifelse(((365+LY-date.yday)<(4-date.wday)),Y+1,Y))

  ##Set the week differently if the date is in the beginning,middle or
  ##end of the year
  Wn<-ifelse(
             Yn==Y-1,
             ifelse((jan1.wday==5|(jan1.wday==6 &LY.prev)),53,52),
             ifelse(Yn==Y+1,1,(date.yday+(7-date.wday)+(jan1.wday-1))/7-(jan1.wday>4))
             )

    return(list(ISOYear=Yn,ISOWeek=Wn)) 
}

######################################################################
# Not very beautiful function implementing a platform independent
# format.Date function. See format.Date, which for the %V and %G
# format strings does not work on windows.
# Added format string %Q for formatting of the quarter (1-4) the month
# belongs to.
#
# Params:
#  x - An object of type Date to be converted.
# format - A character string. Note that only "%V and %G" are
#          processed on Windows. Otherwise the call is sent to format.Date
######################################################################

#Small helper function - vectorized gsub, but disregarding names of x
gsub2 <- function(pattern, replacement, x) {
    len <- length(x)
    mapply(FUN = gsub,
           pattern = rep_len(as.character(pattern), len),
           replacement = rep_len(as.character(replacement), len),
           x = x,
           MoreArgs = list(fixed = TRUE),
           SIMPLIFY = TRUE, USE.NAMES = FALSE)
}

#More general version also handling a mix of several formats
formatDate <- function(x, format) {
  ##Anything to do?
  if (!grepl( "%G|%V|%Q|%OQ|%q", format)) { #nope
    return(format(x,format))
  }

    #Replicate string.
    formatStr <- rep_len(format,length(x))
        
    ##If days within quarter requested (this is kind of slow)
    if (grepl("%q",format)) {
      ##Loop over vectors of dates 
      dateOfQuarter <- sapply(x, function(date) {
        ##Month number in quarter
        modQ <- (as.numeric(format(date,"%m"))-1) %% 3 
        dateInMonth <- seq(date,length.out=2,by=paste0("-",modQ," month"))[2]
        ##Move to first of month
        return(dateInMonth - as.numeric(format(dateInMonth,"%d")) + 1)
      })
      dayInQuarter <- as.numeric(x - dateOfQuarter) + 1
      formatStr <- gsub2("%q",as.character(dayInQuarter),formatStr)
    }
    
    if (grepl("%Q|%OQ",format)) {
      Q <- (as.numeric(format(x,"%m"))-1) %/% 3 + 1 #quarter
      formatStr <- gsub2("%Q",as.character(Q),formatStr)
      formatStr <- gsub2("%OQ",as.roman(Q),formatStr)
    }
    
    if (.Platform$OS.type == "windows") {
      ##Year/week
      isoYear <- isoWeekYear(x)$ISOYear
      isoWeek <- sprintf("%.2d",isoWeekYear(x)$ISOWeek)
      formatStr <- gsub2("%G",isoYear,formatStr)
      formatStr <- gsub2("%V",isoWeek,formatStr)
    }
    
    ##The rest of the formatting - works normally as defined by strptime
    res <- character(length(x))
    for (i in 1:length(x)) { res[i] <- format(x[i],formatStr[i])}
    return(res)
}
