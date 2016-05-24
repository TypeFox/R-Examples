hydro.year <- function(flow.ts, hydro.year = "hydro", year.only = FALSE) {
    
    begin <- head(flow.ts[, "Date"], 1)
    finish <- tail(flow.ts[, "Date"], 1)
    
    if (hydro.year == "hydro") {
        month.runs <- c()
        month.means <- aggregate(flow.ts[, "Q"], by = list(month = strftime(flow.ts[, "Date"], format = "%m")), sum, na.rm = T)
        month.runs[1] <- sum(month.means[1:6, "x"])
        month.runs[2] <- sum(month.means[2:7, "x"])
        month.runs[3] <- sum(month.means[3:8, "x"])
        month.runs[4] <- sum(month.means[4:9, "x"])
        month.runs[5] <- sum(month.means[5:10, "x"])
        month.runs[6] <- sum(month.means[6:11, "x"])
        month.runs[7] <- sum(month.means[7:12, "x"])
        month.runs[8] <- sum(month.means[c(8:12, 1), "x"])
        month.runs[9] <- sum(month.means[c(9:12, 1:2), "x"])
        month.runs[10] <- sum(month.means[c(10:12, 1:3), "x"])
        month.runs[11] <- sum(month.means[c(11:12, 1:4), "x"])
        month.runs[12] <- sum(month.means[c(12, 1:5), "x"])
        
        alt.month <- which.min(month.runs)
        hydro.year <- c()
        hydro.year <- ifelse(as.numeric(strftime(flow.ts[, "Date"], format = "%m")) < alt.month, as.numeric(strftime(flow.ts[, "Date"], format = "%Y")) - 1, as.numeric(strftime(flow.ts[, "Date"], 
            format = "%Y")))
    }
    if (year.only == TRUE) {
        return(hydro.year)
    } else {
        flow.ts <- data.frame(flow.ts, hydro.year = hydro.year)
        return(flow.ts)
    }
    
} 
