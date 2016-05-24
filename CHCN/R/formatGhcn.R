formatGhcn <- function(data, dataColumn = 7){
  
  colSelect <- c(1,3,dataColumn)
  data <- data[,colSelect]  
  Ids  <- unique(data$Id)
  CHCN <- matrix(ncol = 14)
  firstStation <- Ids[1]
  for(stationId in Ids){ 
       stationData <-  data[which(data$Id == stationId),]    
       
       keep  <- which(!is.na(stationData[,2]))
       stationData <- stationData[keep,]
       years <- unique(stationData[,2])
       startYear <- min(years)
       endYear   <- max(years)
       cat(stationId, startYear, endYear, "\n")
       if (sum(diff(years)) != length(years) -1 ) {
         
         warning("gaps in years")
         print(stationId)
         } else { 
       temps <- rep(NA,((endYear-startYear)+1)*12)
       temps[1:nrow(stationData)] <- stationData[,3]  
       tempMat <- matrix(temps,ncol = 12, byrow = TRUE)   
       thisStation <- cbind(stationId,startYear:endYear,tempMat)
       if (stationId == firstStation) {
          CHCN <- thisStation 
          } else {
             CHCN <-rbind(CHCN,thisStation)
           }
         }
  }
   
  colnames(CHCN) <- c("Id","Year", month.abb)
  return(CHCN)
}