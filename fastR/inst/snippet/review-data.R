names(batting)
require(Hmisc)
summary(HR~team,data=batting,fun=max,
        subset=(year==2005&league=="AL"),nmin=1)
plot1 <- histogram(~AB, data=batting,subset=year==2005)
plot2 <- xyplot(HR~H, subset=(team=="DET" & year==2005), data=batting)
plot3 <- bwplot(HR~league, data=batting,subset=year==2005)
