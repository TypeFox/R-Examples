# by Jason P Becker
age_calc <- function(dob, enddate=Sys.Date(), units=c("days","months","years"), precise=TRUE){
    
  if (!(inherits(dob, "POSIXlt") | inherits(dob, "POSIXct") | inherits(dob, "Date"))   | 
      !(inherits(enddate, "POSIXlt") | inherits(enddate, "POSIXct") | inherits(enddate, "Date"))){
    stop("Both dob and enddate must be either POSIXlt or Date class objects")
  }
  
  start <- as.POSIXlt(dob)
  end <- as.POSIXlt(enddate)
  
  if(precise){
    start_is_leap <- ifelse(start$year %% 400 == 0, TRUE, 
                        ifelse(start$year %% 100 == 0, FALSE,
                               ifelse(start$year %% 4 == 0, TRUE, FALSE)))
    end_is_leap <- ifelse(end$year %% 400 == 0, TRUE, 
                        ifelse(end$year %% 100 == 0, FALSE,
                               ifelse(end$year %% 4 == 0, TRUE, FALSE)))
  }
    
  units=match.arg(units)  
  if(units=='days'){
    result <- difftime(end, start, units='days')
  }else if(units=='months'){
    months <- sapply(mapply(seq, as.POSIXct(start), as.POSIXct(end), 
                            by='months', SIMPLIFY=FALSE), 
                     length) - 1
    # length(seq(start, end, by='month')) - 1
    if(precise){
      month_length_end <- ifelse(end$mon==1, 28,
                                 ifelse(end$mon==1 & end_is_leap, 29,
                                        ifelse(end$mon %in% c(3, 5, 8, 10), 
                                               30, 31)))
      month_length_prior <- ifelse((end$mon-1)==1, 28,
                                     ifelse((end$mon-1)==1 & start_is_leap, 29,
                                            ifelse((end$mon-1) %in% c(3, 5, 8, 
                                                                      10), 
                                                   30, 31)))
      month_frac <- ifelse(end$mday > start$mday,
                           (end$mday-start$mday)/month_length_end,
                           ifelse(end$mday < start$mday, 
                            (month_length_prior - start$mday) / 
                                month_length_prior + 
                                end$mday/month_length_end, 0.0))
      result <- months + month_frac
    }else{
      result <- months
    }
  }else if(units=='years'){
    years <- sapply(mapply(seq, as.POSIXct(start), as.POSIXct(end), 
                            by='years', SIMPLIFY=FALSE), 
                     length) - 1
    if(precise){
      start_length <- ifelse(start_is_leap, 366, 365)
      end_length <- ifelse(end_is_leap, 366, 365)
      year_frac <- ifelse(start$yday < end$yday,
                          (end$yday - start$yday)/end_length,
                          ifelse(start$yday > end$yday, 
                                 (start_length-start$yday) / start_length +
                                end$yday / end_length, 0.0))
      result <- years + year_frac
    }else{
      result <- years
    }
  }
  
  return(result)
}
