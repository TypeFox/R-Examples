#' Calculate baseflow peak statistics
#' 
#' This function finds the start, middle, end, and duration of the baseflow
#' peak based on percent of the total annual baseflow volume.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param bfpct numeric vector of percentages used to define the 
#'   start, middle, and end of the baseflow peak. Default is c(25, 50, 75)
#' @details This function calculates metrics intended to focus on snowmelt-related
#'   streamflow occuring in spring and summer. For catchments in cold climates, 
#'   the baseflow peak can be interpreted as snowmelt-induced. Baseflow is estimated with 
#'   \code{\link{bf_eckhardt}}.
#' @return Returns a data.frame with the following columns:
#'   \itemize{
#'     \item Start - day of year defining the start of the baseflow peak
#'     \item Mid - day of year defining the middle of the baseflow peak
#'     \item End - day of year defining the end of the baseflow peak
#'     \item Dur - duration of the baseflow peak, in days
#'   }
#' @author Jennifer Dierauer
#' @export
#' @examples
#' data(cania.sub.ts)
#' res1 <- pk.bf.stats(cania.sub.ts)
#' 
#' # trend and changepoint plot for baseflow peak start doy
#' res2 <- screen.metric(res1[,1], "Day of Year")


pk.bf.stats <- function(TS, bfpct=c(25,50,75)) {
    
    #create factors for year
    Year <- as.factor(TS$hyear)
    
    ### Set parameter values for Eckhardt RDF
    BFindex <- 0.8
    alpha <- 0.970 ##based on values suggested by Eckhardt 2012 for perrennial stream
    
    ## calculate daily BF and BFI
    temp <- TS
    temp$ro <- temp$Flow - bf_eckhardt(temp$Flow, alpha, BFindex)
    
    ## calculate 10%, 50%, and 90% of annual baseflow volume
    BFVy <- array(data=NA, c(max(as.numeric(Year)),length(bfpct)+1))
    colnames(BFVy) <- c("Sum", bfpct)
    BFVy[,1] <- tapply(temp$ro, Year, sum)
    Output <- data.frame(1:length(unique(Year)))
    for (i in 1:length(bfpct)) {
        BFVy[,i+1]<-BFVy[,1]*(bfpct[i]/100)
        Output[,i] <- NA 
    }
    
    colnames(Output) <- bfpct
    
    YearStack <- split(temp, Year)  #split into dataframes by year for looping
    
    for (i in 1:length(YearStack)) { #loop through years
        
        Have10 <- FALSE
        Have50 <- FALSE
        Have90 <- FALSE
        mysum <- 0
        temp <- YearStack[[i]]
        j <- 0
        
        while (Have90 == FALSE) { #loop through days until Q90 value is reached
            j <- j+1
            
            mysum <- mysum + temp$ro[j]
            
            if (mysum > BFVy[i,2] & Have10 == FALSE){
                Have10 <- TRUE
                Output[i,1] <- as.character(temp$Date[j])}
            
            else if (mysum > BFVy[i,3] & Have50 == FALSE){
                Have50 <- TRUE
                Output[i,2] <- as.character(temp$Date[j])}
            
            else if (mysum > BFVy[i,4] & Have90 == FALSE){
                Have90 <- TRUE
                Output[i,3] <- as.character(temp$Date[j])}
        }
    }
    
    
    #change character dates to doy for plotting
    for (i in 1:3) {
        temp <- as.Date(Output[,i],format="%Y-%m-%d")
        Output[,i] <- as.numeric(format(temp, "%j"))
    }
    
    #calculate spring flood duration
    Output$Dur<-Output[,3] - Output[,1]
    
    for (i in 1:4) {
        attr(Output[,i], "times") <- as.character(unique(Year))
    }
    
    colnames(Output) <- c("Start", "Mid", "End", "Dur")
    return(Output)
    
}