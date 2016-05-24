#' Returns plot titles and labels based on plot type and language preference
#' 
#' @param type character indicating the type of summary plot
#' @param language "English" or "French"
#' @param Qmax the flow quantile used to define peaks of threshold, e.g. 0.95
#' @author Jennifer Dierauer


get.titles.internal <- function(type, language="English", Qmax) {
    
    tQmax <- strsplit(as.character(Qmax),"\\.")
    tQmax <- tQmax[[1]]
    
    Titles.High.Eng <- c("NA", "NA", "NA", "Annual Maximum Series",
                         "Day of Annual Maximum",
                         paste("Peaks Over Threshold (Q", tQmax[2], ")", sep=""),
                         "Inter-Event Duration", "Q80", "Q90", 
                         "Day of Year 25 pct Annual Flow",
                         "Center of Volume", "Day of Year 75 pct Annual Flow", 
                         "25 pct - 75 pct Annual Flow Duration")
    
    Titles.Low.Eng <- c("NA", "NA", "NA", "Q10", "Q25", "Drought Start",
                        "Drought Center", "Drought End", "Drought Duration",
                        "Drought Severity", "Annual Minimum Flow",
                        "Mean Annual Minimum 7-day Flow",
                        "Mean Annual Minimum 10-day Flow")
    
    Titles.BF.Eng <- c("NA", "NA", "NA", "Annual Mean Daily Discharge", "Annual Baseflow Volume",
                       "Annual Mean Baseflow", "Annual Maximum Baseflow", "Annual Minimum Baseflow",
                       "Mean Baseflow Index", "Day of Year 25 pct Annual Baseflow",
                       "Baseflow Center of Volume", "Day of Year 75 pct Annual Baseflow",
                       "25 pct - 75 pct Baseflow Duration")
    
    Titles.High.Fr <- c("NA", "NA", "NA", "maximum annuel due debit",
                        "jour du maximum",
                        paste("pic au-dessus du seuil (Q", tQmax[2], ")", sep=""),
                        "duree inter-evenements", "Q80", "Q90", "debut de la crue",
                        "centre de la crue", "fin de la crue",
                        "Duree de la crue")
    
    Titles.Low.Fr <- c("NA", "NA", "NA", "Q10", "Q25", "debut de la secheresse",
                       "centre de la secheresse", "fin de la secheresse", 
                       "duree de la secheresse",
                       "gravite de la secheresse", "debit annuel minimum",
                       "debit moyen annuel de 7 jours minimum",
                       "debit moyen annuel de 10 jours minimum")
    
    Titles.BF.Fr <- c("NA", "NA", "NA", "debit moyen annuel", "Volume SDB",
                      "moyenne journaliere SDB", "maximum annuel du SDB", 
                      "minimum annuel du SDB",
                      "moyen BFI", "debut de la crue du SDB",
                      "centre de la crue du SDB", "fin de la crue du SDB",
                      "duree de la crue du SDB")
    
    Xlabs.Eng <- c(1:13)
    Xlabs.Eng[1] <- NA
    Xlabs.Eng[2] <- "Day of Year"
    Xlabs.Eng[3:13] <- "Year"
    
    Xlabs.Fr <- c(1:13)
    Xlabs.Fr[1] <- NA
    Xlabs.Fr[2] <- "Jour du An"
    Xlabs.Fr[3:13] <- "An"
    
    y1 <- expression(paste("Discharge (m" ^{3}, "/s)"))
    y2 <- "No. of Days"
    y3 <- "Day of Year"
    y4 <- expression(paste("BFS (m" ^{3}, "/s)"))
    y5 <- expression(paste("km" ^{3}))
    y6 <- "BFI"
    y7 <- "Percent"
    
    y2f <- "jour de l'annee"
    y3f <- "Jour du An"
    y7f <- "pour cent"
    
    Ylabs.High.Eng <-c(NA, y1, y1, y1, y3, y1, y2, y1, y1, y3, y3, y3, y2)
    Ylabs.Low.Eng <- c(NA, y1, y1, y1, y1, y3, y3, y3, y2, y1, y1, y1, y1)
    Ylabs.BF.Eng <- c(NA, y1, y1, y1, y5, y1, y1, y1, y6, y3, y3, y3, y2)

    Ylabs.High.Fr <-c(NA, y1, y1, y1, y3f, y1, y2f, y1, y1, y3f, y3f, y3f, y2f)
    Ylabs.Low.Fr <- c(NA, y1, y1, y1, y1, y3, y3f, y3f, y2f, y1, y1, y1, y1)
    Ylabs.BF.Fr <- c(NA, y1, y1, y1, y5, y1, y1, y1, y6, y3f, y3f, y3f, y2f)

    
    if (type=="h") {
        
        if (language=="English") {
            output <- list(Xlabs=Xlabs.Eng, Ylabs=Ylabs.High.Eng, Titles=Titles.High.Eng)
        }
        if (language=="French") {
            output <- list(Xlabs=Xlabs.Fr, Ylabs=Ylabs.High.Fr, Titles=Titles.High.Fr)
        }
    }
    
    if (type=="l") {
        
        if (language=="English") {
            output <- list(Xlabs=Xlabs.Eng, Ylabs=Ylabs.Low.Eng, Titles=Titles.Low.Eng)
        }
        if (language=="French") {
            output <- list(Xlabs=Xlabs.Fr, Ylabs=Ylabs.Low.Fr, Titles=Titles.Low.Fr)
        }
    }
    
    if (type=="b") {
        
        if (language=="English") {
            output <- list(Xlabs=Xlabs.Eng, Ylabs=Ylabs.BF.Eng, Titles=Titles.BF.Eng)
        }
        if (language=="French") {
            output <- list(Xlabs=Xlabs.Fr, Ylabs=Ylabs.BF.Fr, Titles=Titles.BF.Eng)
        }
    }
    
    return(output)
    
}