plot1 <- xyplot( ccf ~ (year + month/12), utilities, groups=month )
plot2 <- bwplot( ccf ~ factor(month), utilities )
