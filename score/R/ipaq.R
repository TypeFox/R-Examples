
ipaq <- function(ipaqdata){
  
  data <- data.frame(matrix(NA,nrow=nrow(ipaqdata)))
  
  data$ID <- ipaqdata[,1]
  data$weight <- ipaqdata[,2]
  data$VigDays <- ipaqdata[,3] 
  data$VigHours <- ipaqdata[,4] 
  data$VigMin <- ipaqdata[,5] 
  data$ModDays <- ipaqdata[,6] 
  data$ModHours <- ipaqdata[,7] 
  data$ModMin <- ipaqdata[,8]                 
  data$WalkDays <- ipaqdata[,9] 
  data$WalkHours <- ipaqdata[,10] 
  data$WalkMin <- ipaqdata[,11] 
  data$SitHours <- ipaqdata[,12] 
  data$SitMin <- ipaqdata[,13]
  data <- data[,-1]
  
  vminday <- data$VigHours*60 + data$VigMin 
  mminday <- data$ModHours*60 + data$ModMin 
  wminday <- data$WalkHours*60 + data$WalkMin 
  sminday <- data$SitHours*60 + data$SitMin
  
  vminwk <- data$VigDays*vminday
  mminwk <- data$ModDays*mminday
  wminwk <- data$WalkDays*wminday
  sminwk <- sminday*7
  
  sumpa <- vminday + mminday + wminday
  sumday <- data$VigDays + data$ModDays + data$WalkDays
  
  data$excMin <- ifelse(sumpa > 960, 1,0)
  data$excDay <- ifelse(sumday > 9, 1,0)
  
  vminday <- ifelse(vminday < 10,0,vminday)
  mminday <- ifelse(mminday < 10,0,mminday)
  wminday <- ifelse(wminday < 10,0,wminday)
  
  vminwkMET =  8*vminwk
  mminwkMET =  4*mminwk
  wminwkMET = 3.3*wminwk
  
  MET = vminwkMET + mminwkMET + wminwkMET
  kilocalories = MET*(data$weight/60)
  
  data$MET <- MET
  data$kilocalories <- kilocalories
  
  data$pacat[data$VigDays >= 3 & vminday >= 20] <- 'Moderate'
  data$pacat[data$ModDays >= 5] <- 'Moderate'
  data$pacat[wminday >= 30] <- 'Moderate'
  data$pacat[(data$VigDays + data$ModDays + data$WalkDays)>= 5 & MET>= 600] <- 'Moderate'
  
  data$pacat[data$VigDays >= 3 & MET >= 1500] <- 'High'
  data$pacat[(data$VigDays + data$ModDays + data$WalkDays)>= 7 & MET>= 3000] <- 'High'
  
  data$pacat[is.na(data$pacat)] <- 'Low';
  data$pacat[is.na(data$VigDays) & is.na(data$VigHours) & is.na(data$VigMin) &
               is.na(data$ModDays) & is.na(data$ModHours) & is.na(data$ModMin) &
               is.na(data$WalkDays) & is.na(data$WalkHours) & is.na(data$WalkMin) ] <- NA
  
  data <- data[ which(data$excMin!=1 & data$excDay!=1), ]
  data <- data[,c(-14,-15)]
  rownames(data) <- NULL
  data
}