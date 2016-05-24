#' Partial Duration Series and Event Statistics for streamflow droughts
#' 
#' This function extracts the partial duration series for all streamflow droughts based
#' on a moving window quantile threshold.  Also returns summary information (start date,
#' end date, duration, deficit volume) for each drought event.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param Qdr Numeric value of the drought threshold quantile.  Default is 0.2.
#' @param WinSize Numeric value specifying the size of the moving window in
#'   days.  Default is 30. 
#' @param IntEventDur Numeric value for the minimum inter-event duration in days. Drought
#'   events with less than the specified number of days between will be pooled 
#'   and considered as one event.
#' @param EventDur Numeric value for the minimum drought duration in days. Default
#'   is 15.
#' @return Returns a list with the following elements:
#' 
#'   DroughtEvents: A data.frame with the following columns:
#'   \itemize{
#'     \item Event - Integer indicating the original event number assigned before
#'       minor drought events were removed.
#'     \item Start - Date of the start of the drought event.
#'     \item End - Date of the end of the drought event
#'     \item maxDef - Numeric value of the maximum streamflow deficit.
#'     \item Severity - Numeric value indicating the drought severity, 
#'       calculated as the cumulative daily streamflow deficit in m3/s.
#'     \item Duration - Numeric value of the drought duration in days.
#'     \item Magnitude - Numeric value indicating the drought magnitude, which 
#'       is calculated as the mean daily streamflow deficit in m3/s.
#'     \item stdtotDef - Numeric value indicating the standardized cumulative
#'       streamflow deficit, calculated as the drought severity divided by
#'       the mean annual daily streamflow.
#'   }
#'   
#'   DroughtPDS: A data.frame of the original input TS that has been subset to 
#'   include only the days on which the streamflow was below the drought threshold.
#'   The data.frame also has the following columns appended:
#'   \itemize{
#'     \item Thresh - Numeric value indicating the streamflow drought threshold,
#'       as calculated by \code{\link{mqt}}
#'     \item BelowThresh - Logical indicating whether the observed streamflow 
#'       was below the streamflow drought threshold.
#'     \item Def - Numeric value of the streamflow defict, calculated as the 
#'       streamflow drought threshold (m3/s) minus the observed streamflow (m3/s).
#'   }
#' @author Jennifer Dierauer
#' @seealso See \code{\link{dr.seas}} to calculate metrics for droughts 
#'   occurring in a user-defined season.
#'   
#'   This function calls \code{\link{dr.pds}} which calls \code{\link{mqt}}.
#' @export
#' @examples
#' data(cania.sub.ts)
#' res1 <- dr.events(cania.sub.ts)
#' events <- res1$DroughtEvents
#' plot(events$Start, events$Duration, pch=19, ylab="Drought Duration (days)", xlab="")

dr.events <- function(TS, Qdr=0.2, WinSize=30, IntEventDur=10, EventDur=15) {
        
    ## create PDS from original data
    myPDS <- dr.pds(TS, Qdr, WinSize)
    
    ### subset to keep only dates when streamflow was below threshold
    myPDS <- subset(myPDS, myPDS$BelowThresh==TRUE)

    if (length(myPDS[,1]) > 0) {
    
        ### pool drought events based on IntEventDur
        eventNo <- 1
        myPDS$Event <- 1
        
        for (i in 2:length(myPDS$Flow)) {
            
            myPDS$Event[i] <- eventNo
            
            # if duration between two dates is greater than or equal to the
            # minimum event duration, define as a new event
            meventdur <- myPDS$Date[i] - myPDS$Date[i-1]
            meventdur <- as.numeric(meventdur)
                        
            if (meventdur >= IntEventDur) {
                eventNo <- eventNo + 1
            }
            
            myPDS$Event[i] <- eventNo
            
        }
        
        myPDS$Def <- as.numeric(myPDS$Thresh - myPDS$Flow)
        
        ### set up output data.frame 
        DroughtEvents <- data.frame(Event=numeric(), 
                                    Start=as.Date(character()),
                                    End=as.Date(character()),
                                    maxDef=numeric(),
                                    totDef=numeric(),
                                    Duration=numeric())

        ### calculate event stats (start, end, duration, maxDef, totDef) 
        ### ignore drought events with duration < EventDur
        for (j in 1:eventNo) {
            PDS.sub <- subset(myPDS, myPDS$Event==j)
            dur <- length(PDS.sub[,1])
            if (dur >= EventDur) {
                
                dstats <- data.frame(Event=j, 
                                     Start=as.Date(PDS.sub$Date[1]),
                                     End=as.Date(PDS.sub$Date[dur]),
                                     maxDef=max(PDS.sub$Def), 
                                     Severity=sum(PDS.sub$Def),
                                     Duration=as.numeric(dur),
                                     Magnitude=mean(PDS.sub$Def),
                                     stdtotDef=(sum(PDS.sub$Def)/mean(TS$Flow)))
                
                DroughtEvents <- rbind(DroughtEvents, dstats)
            }
        }
        
        ### subset PDS to remove minor droughts
        Events <- unique(DroughtEvents[,1])
        myPDS <- subset(myPDS, myPDS$Event %in% Events)
        
        output <- list(DroughtEvents=DroughtEvents, DroughtPDS=myPDS)
        
    } else {output <- list(DroughtEvents=NA, DroughtPDS=NA)}

return(output)

}
