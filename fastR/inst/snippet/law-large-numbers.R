x <- rexp(1000)
runningMean <- cumsum(x) / 1:length(x)
expPlot <- 
    xyplot(runningMean~1:1000,ylab="running mean", xlab="n",type="l",
        panel=function(...){ panel.abline(h=1,col='gray70');
            panel.xyplot(...); });
