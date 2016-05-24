utilities$ccfpday <- utilities$ccf / utilities$billingDays
plot1 <- xyplot( ccfpday ~ (year + month/12), utilities, groups=month )
plot2 <- bwplot( ccfpday ~ factor(month), utilities )
