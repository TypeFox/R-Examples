seasonality <- function(flow.ts, monthly.range = FALSE) {
    
    
    month.runs <- c()
    month.means <- aggregate(flow.ts[, "Q"], by = list(month = strftime(flow.ts[["Date"]], format = "%m")), mean, na.rm = T)
    month.runs[1] <- sum(month.means[1:6, "Q"], na.rm = T)
    month.runs[2] <- sum(month.means[2:7, "Q"], na.rm = T)
    month.runs[3] <- sum(month.means[3:8, "Q"], na.rm = T)
    month.runs[4] <- sum(month.means[4:9, "Q"], na.rm = T)
    month.runs[5] <- sum(month.means[5:10, "Q"], na.rm = T)
    month.runs[6] <- sum(month.means[6:11, "Q"], na.rm = T)
    month.runs[7] <- sum(month.means[7:12, "Q"], na.rm = T)
    month.runs[8] <- sum(month.means[c(8:12, 1), "Q"], na.rm = T)
    month.runs[9] <- sum(month.means[c(9:12, 1:2), "Q"], na.rm = T)
    month.runs[10] <- sum(month.means[c(10:12, 1:3), "Q"], na.rm = T)
    month.runs[11] <- sum(month.means[c(11:12, 1:4), "Q"], na.rm = T)
    month.runs[12] <- sum(month.means[c(12, 1:5), "Q"], na.rm = T)
    
    month.runs.sort <- sort(month.runs)
    low.6.months <- sum(month.runs.sort[1:6])
    out <- low.6.months/sum(month.runs) * 100
    
    
    if (monthly.range == TRUE) {
        flow.ts$Year <- as.factor(strftime(flow.ts$Date, format = "%Y"))
        flow.ts$month <- as.factor(strftime(flow.ts$Date, format = "%m"))
        month.year.means <- tapply(flow.ts[, "Q"], flow.ts[c("Year", "month")], sum)
        month.year.means <- na.omit(month.year.means)
        monthly.means <- apply(month.year.means, 2, mean, na.rm = T)
        month.range <- apply(month.year.means, 1, range, na.rm = T)
        min.months <- apply(month.year.means, 1, which.min)
        max.months <- apply(month.year.means, 1, which.max)
        av.ann.month.range <- mean(month.range[2, ] - month.range[1, ])
        month.difs <- factor(as.character(abs(max.months - min.months)))
        month.difs <- as.factor(recode(month.difs, oldvalue = c("11", "10", "9", "8", "7"), newvalue = c("1", "2", "3", "4", "5")))
        
        return(list(seasonality = out, monthly.means, avg.ann.month.range = av.ann.month.range, max.min.time.dif = round(mean(as.numeric(levels(month.difs)[month.difs])), 0)))
        
    } else {
        return(out)
    }
} 
