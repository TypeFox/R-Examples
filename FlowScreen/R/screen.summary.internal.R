#' Internal wrapper for creating trend and change-point summary plots
#' 
#' @param input metric time series
#' @param mparam trend and change point info
#' @param type character indicating type of summary plot
#' @param ylabs y axis labels
#' @param i plot position
#' @param DataType numeric indicating data type
#' @param maf mean annual flow series
#' @param Year1 start year of original time series
#' @param YearEnd end year of original time series
#' @param hyrstart numeric indicating month for start of the hydrologic year
#' @author Jennifer Dierauer


screen.summary.internal <- function(input, mparam, type, ylabs, i, DataType, maf, Year1, YearEnd,
                       hyrstart) {
    
    # set colors for data types
    cols <- c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6")
    
    # get plot number
    PlotNo <- i + 3
    
    MyY <- input
    MyX <- attr(MyY, "times")
    metricID <- mparam[[1]]
    
    NumYrs <- (YearEnd - Year1) + 1
    
    # convert x values to work with the slopes
    if (nchar(as.character(MyX[1])) > 5) {

        Start <- as.Date(paste(Year1, "-01-01", sep=""))
        MyX.mod <- c(1:length(MyX))
        
        for (j in 1:length(MyX)) {MyX.mod[j] <- (MyX[j]-Start)}
        MyXMax <- as.Date(paste(YearEnd, "-12-31", sep="")) - Start

    } else {
        
        MyX.mod <- c(1:length(MyX))
        for (j in 1:length(MyX)) {MyX.mod[j] <- (as.numeric(MyX[j]) - Year1) + 1}
        MyXMax <- NumYrs
        
    }

    ### set y axis limits based on y-axis label
    ifelse(identical(as.character(ylabs[i+3]), "Day of Year"), MyYlims <- c(0, 365),
                      ifelse(identical(as.character(ylabs[i+3]), "BFI"), MyYlims <- c(0, 1), 
                             MyYlims <- c(0, ceiling(1.2*max(MyY, na.rm=T)))))

    # set figure margins based on position in figure
    if (PlotNo %in% c(4,5,6,7,8)) {graphics::par(mar=c(0.5,4,0,0))}
    if (PlotNo == 9) {graphics::par(mar=c(0.5,4,0,0.5))}
    if (PlotNo %in% c(10,11,12)) {graphics::par(mar=c(3,4,0,0))}
    if (PlotNo == 13) {graphics::par(mar=c(3,4,0,0.5))}
    
    
    ## set y axis labels based on hydrologic year for cov plots
    if (metricID %in% c(2, 7, 8, 9, 27, 28, 29)) {
        
        graphics::plot(MyX.mod, MyY,
             ylab="", ylim=MyYlims, yaxt="n",
             xlim=c(0, MyXMax), xaxt="n", xlab="",
             type="p", col=cols[DataType[i]],
             cex=0.8, pch=19)
        
        if (hyrstart != 1) {
            mlabels <- c(month.abb[hyrstart:12], month.abb[1:(hyrstart-1)])
        } else {mlabels <- month.abb}
        
        myticks <- c(1, 153, 336)
        mlabels <- c(mlabels[1], mlabels[6], mlabels[12])
        graphics::axis(2, at=myticks, labels=mlabels)
        
    } else {
        
        graphics::plot(MyX.mod, MyY,
             ylab="", ylim=MyYlims, 
             xlim=c(0, MyXMax),
             xaxt="n", xlab="",
             type="p", col=cols[DataType[i]],
             cex=0.8, pch=19)
        
        graphics::title(ylab=ylabs[i+3], line=2)
    }
    
    
    # if plot is in bottom row, add axis labels
    if (PlotNo %in% c(11,12,13)) {
        
        ifelse(NumYrs > 100, mby <- 20, ifelse(NumYrs > 50, mby <- 15, mby <- 10))

        myticks <- seq(from = 1, to = max(MyX.mod), by = mby)
        mlabels <- seq(from = Year1, to = YearEnd, by = mby)
        graphics::axis(1, at=myticks, labels=mlabels)
    }
    
    if (PlotNo == 10) {
        if (type=="l") {
                
            series.length <- YearEnd - Year1
            ifelse(NumYrs > 100, mby <- 20, ifelse(series.length > 50, mby <- 15, mby <- 10))
            myticks <- seq(from=1, to=365.25*(series.length + 2), by = (365.25 * mby))
            mlabels <- seq(from = Year1, to = YearEnd, by = mby)
            graphics::axis(1, at=myticks, labels=mlabels)

        } else {

            ifelse(NumYrs > 100, mby <- 20, ifelse(NumYrs > 50, mby <- 15, mby <- 10))
            myticks <- seq(from = 1, to = max(MyX.mod), by = mby)
            mlabels <- seq(from = Year1, to = YearEnd, by = mby)
            graphics::axis(1, at=myticks, labels=mlabels)
            
        }
    }
    
    ### title placement based on value of first data point
    ypos <- ifelse(mean(MyY[1:3]) <= 0.5*MyYlims[2], 0.95*MyYlims[2], 0.05*MyYlims[2])
    graphics::text(1, ypos, mparam[[2]], adj=c(0,0))
    
    ## add trend lines and color according to increasing or decreasing
    slope <- mparam[[3]]
    
    if (!is.na(mparam[[6]]) && mparam[[6]] <= 0.1) {
        pval <- mparam[[6]]
        mcol <- ifelse(slope[2] < 0, "darkred", "darkblue")
        mlwd <- ifelse(pval <= 0.05, ifelse(pval<=0.01, 3, 2), 1)
        graphics::abline(coef=slope, col=mcol, lwd=mlwd)
        graphics::abline(coef=mparam[[4]], lty=3, col=mcol)
        graphics::abline(coef=mparam[[5]], lty=3, col=mcol)
        ypos2 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.1*MyYlims[2]), ypos - (0.1*MyYlims[2]))
        mypvalue <- ifelse(pval >= 0.01, round(as.numeric(pval), digits=2), "< 0.01")
        graphics::text(1, ypos2, paste("p = ", mypvalue),col=mcol, adj=c(0,0))
    }
  
    MyCpts <- mparam[[7]]
    MyMeans <- mparam[[8]]
    NumPoints <- length(MyX.mod)

    MyCpts <- MyCpts[MyCpts < NumPoints]
    
    # add changepoints and means
    if (!is.na(MyCpts) && length(MyCpts) > 0){
        
        for (j in 1:length(MyCpts)) {
            graphics::abline(v=MyX.mod[MyCpts[j]], lwd=2, lty=5)
            ifelse(j==1, xpts <- c(-10, MyX.mod[MyCpts[j]]),
                        xpts <- c(MyX.mod[MyCpts[j-1]], MyX.mod[MyCpts[j]]))
            ypts <- c(MyMeans[j], MyMeans[j])
            graphics::points(xpts, ypts, type="l", lwd=2)
        }
        
        xpts <- c(xpts[2], 1.1*max(MyX.mod))
        ypts <- c(MyMeans[length(MyMeans)], MyMeans[length(MyMeans)])
        graphics::points(xpts, ypts, type="l", lwd=2)
    }
    
    ## add r2 info to cov plots
    if (PlotNo == 11 && type %in% c("h", "b")) {
        
        covyrs <- stats::cor.test(MyX.mod, MyY, use="pair", method="pearson")
        mafyrs <- stats::cor.test(as.numeric(attr(maf, "times")), maf, use="pair", method="pearson")
        covR <- covyrs$estimate
        mafR <- mafyrs$estimate
        
        ratio <- round(covR/mafR, 2)
        covR <- round(covR, 2)
        mafR <- round(mafR, 2)
        
        # text placement based on plot titles and p-value text placement
        
        if (mparam[[6]] <= 0.1) {
            ypos3 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.2*MyYlims[2]), ypos - (0.2*MyYlims[2]))
            ypos4 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.3*MyYlims[2]), ypos - (0.3*MyYlims[2]))
            ypos5 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.4*MyYlims[2]), ypos - (0.4*MyYlims[2]))
        } else {
            ypos3 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.1*MyYlims[2]), ypos - (0.1*MyYlims[2]))
            ypos4 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.2*MyYlims[2]), ypos - (0.2*MyYlims[2]))
            ypos5 <- ifelse(ypos < 0.5*MyYlims[2], ypos + (0.3*MyYlims[2]), ypos - (0.3*MyYlims[2]))
        }

        
        graphics::text(1, ypos3, paste("Ratio", ratio), adj=c(0,0))
        ifelse(covyrs$p.value <= 0.05, pcol1 <- "red", pcol1 <- "black")
        ifelse(mafyrs$p.value <= 0.05, pcol2 <- "red", pcol2 <- "black")
        
        graphics::text(1, ypos4, paste("r cov-years:", covR), col=pcol1, adj=c(0,0))
        graphics::text(1, ypos5, paste("r maf-years:", mafR), col=pcol2, adj=c(0,0)) 
    }
}

