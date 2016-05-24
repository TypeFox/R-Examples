summaryfunction = function(ftsdata, plot.type = c("summarystats", "quantilestats"),
                       quantilepercent = seq(0.1,0.9,by=0.1), plot.legend = FALSE,
                       legendpos = "topright", cex = 0.9, lwd = 1, lty = 1, ncol = 2)
{
    plot.type=match.arg(plot.type)
    if(all(class(ftsdata) != "fds"))
    {
        stop("Data are not functional data class.")
    }
    p = dim(ftsdata$y)[1]
    xname = ftsdata$xname
    yname = ftsdata$yname
    if(plot.type == "summarystats")
    {
        minvec = firstvec = medvec = meanvec = thirdvec = maxvec = vector(,p)
        for(i in 1:p)
        {
            dummy = summary(ftsdata$y[i,])
            minvec[i] = dummy[[1]]
            firstvec[i] = dummy[[2]]
            medvec[i] = dummy[[3]]
            meanvec[i] = dummy[[4]]
            thirdvec[i] = dummy[[5]]
            maxvec[i] = dummy[[6]]
        }
        summa = cbind(minvec, firstvec, medvec, meanvec, thirdvec, maxvec)
        matplot(ftsdata$x, summa, type = "l", lty = 1:6, xlab = xname, ylab = yname)
        if(plot.legend==TRUE)
        {
            legend(legendpos, c("Minimum", "1st quart", "Median", "Mean", "3rd quarter", "Maximum"), col = 1:6,
                   cex = cex, lwd = lwd, lty = lty, ncol=ncol)
        }
    }
    else
    {
        quantvec=matrix(,p,length(quantilepercent))
        for(i in 1:p)
        {
            quantvec[i,] = quantile(ftsdata$y[i,], prob = quantilepercent)
        }
        matplot(ftsdata$x, quantvec, type="l", col = rainbow(length(quantilepercent)),
                lty = 1:length(quantilepercent), xlab = xname, ylab = yname)
        if(plot.legend == TRUE)
        {
            legend(legendpos, as.character(quantilepercent), col = rainbow(length(quantilepercent)),
                   cex = cex, lwd = lwd, lty = lty, ncol = ncol)
        }
    }
}  
