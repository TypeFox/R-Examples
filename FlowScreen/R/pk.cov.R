#' Center of Volume
#' 
#' This function calculates center of volume metrics, including the day of the 
#' hydrologic year that 25 percent, 50 percent, and 75 percent of the total annual 
#' streamflow is reached. 
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @return Returns a data.frame with the following columns:
#'   \itemize{
#'     \item hYear - Hydrologic Years
#'     \item Q25 - day of hydrologic year for 25 percent of the total annual streamflow
#'     \item Q50 - day of hydrologic year for 50 percent of the total annual streamflow, i.e. Center of Volume
#'     \item Q75 - day of hydrologic year for 75 percent of the total annual streamflow
#'     \item Dur - duration of between the 25 percent and 75 percent day of year, in days
#'   }
#' @author Jennifer Dierauer
#' @export
#' @examples
#' data(cania.sub.ts)
#' res1 <- pk.cov(cania.sub.ts)
#' 
#' # trend and changepoint plot for baseflow peak start doy
#' res2 <- screen.metric(res1[,2], "Day of Year")

pk.cov <- function(TS) {
    
    Year <- as.factor(TS$hyear)
    NumRecords<-tapply(TS$Flow, Year, length)
    YearList <- unique(Year)
    YearList.sub <- YearList[NumRecords >= 60]
    TS <- TS[Year %in% YearList.sub,]
    
    #create factors for year
    Year <- as.factor(TS$hyear)
    
    year_list <- unique(Year)
    
    out <- data.frame(hYear=year_list, Q25=NA, Q50=NA, Q75=NA, Dur=NA)
    
    YearStack <- split(TS, Year)  #split into dataframes by year for looping
    
    for (i in 1:length(YearStack)) { #loop through years
        
        Have25 <- FALSE
        Have50 <- FALSE
        Have75 <- FALSE
        mysum <- 0
        temp <- YearStack[[i]]
        j <- 0
        
        
        while (Have75 == FALSE) { #loop through days until 75% is reached
            j <- j+1
            
            total.flow <- sum(temp$Flow)
            
            q25 <- total.flow * 0.25
            q50 <- total.flow * 0.5
            q75 <- total.flow * 0.75
            
            mysum <- mysum + temp$Flow[j]
            
            if (mysum > q25 & Have25 == FALSE){
                Have25 <- TRUE
                out[i,2] <- as.numeric(temp$hdoy[j])
            }
            
            else if (mysum > q50 & Have50 == FALSE){
                Have50 <- TRUE
                out[i,3] <- as.numeric(temp$hdoy[j])
            }
            
            else if (mysum > q75 & Have75 == FALSE){
                Have75 <- TRUE
                out[i,4] <- as.numeric(temp$hdoy[j])
            }
        }
    }
    
    out$Dur <- out[,4] - out[,2]
    
    for (i in 2:5) {
        attr(out[,i], "times") <- as.character(year_list)
    }
    
    return(out)
    
}