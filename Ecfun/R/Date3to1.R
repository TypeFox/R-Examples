Date3to1 <- function(data, default='Start'){
## 
## 1.  check   
##  
  nc <- ncol(data)
  if(is.null(nc)){
    stop('data is not a data.frame')
  }
  if(nc != 3){
    stop('ncol(data) = ', nc, ' != 3')
  }
##
## 2.  defStart
##
  nchd <- nchar(default)
  if(nchd<1){
    stop('nchar(default) < 1:  erroneous call')
  }
  def1 <- toupper(substring(default, 1, 1))
  defStart <- (def1 == 'S')
  defSt1 <- (1+defStart)
##
## 3.  Character vector of dates 
##
#  3.1.  is.na(YEAR) 
  Dt <- as.list(data)
  YrNA <- (is.na(Dt[[1]]) | (Dt[[1]]<1))
#  3.2.  Month <1 or >12
  MoNA <- which(is.na(Dt[[2]]) |     
        (Dt[[2]]<1) | (Dt[[2]]>12))
  Dt[[2]][MoNA] <- c(12, 1)[defSt1]
#  3.3.  Days in each month used  
  Mo1 <- Dt[[2]]+1
  Mo1[Mo1>12] <- 1 
  YM1ch <- paste(Dt[[1]], Mo1, "01", sep='-')
  YM1ch[YrNA] <- NA 
  YMend <- (as.Date(YM1ch)-1) 
  daysofmonth <- as.numeric(substring(YMend, 9, 10))  
  dayout <- which(is.na(Dt[[3]]) |    
        (Dt[[3]]<1) | (daysofmonth < Dt[[3]]))
  if(defStart){
#  start     
    Dt[[3]][MoNA] <- 1
    Dt[[3]][dayout] <- 1  
  } else {
#  end    
    Dt[[3]][MoNA] <- daysofmonth[MoNA]
    Dt[[3]][dayout] <- daysofmonth[dayout]    
  }
#  3.4.  paste -> char
  Dt$sep <- "-"
  Dte <- do.call(paste, Dt)
  Dte[YrNA] <- NA
##
## 3.  attribute 
## 
  msng <- YrNA 
  msng[MoNA] <- TRUE 
  msng[dayout] <- TRUE 
##
## 4.  Finish  
##
  DTE <- as.Date(Dte)  
  if(any(msng)){
    attr(DTE, 'missing') <- which(msng)
  }
  DTE 
}
