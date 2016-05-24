##================================================================
## Author : Y. Deville
##
## This function performs a Goodness-Of-Fit test for date of
## events. Homogeneous Poisson Process?
##
## It was heavily modified from versions > 0.5-2
##
##================================================================

gof.date <- function(date,
                     start = NULL,
                     end = NULL,
                     plot = TRUE,
                     main = NULL,
                     skip = NULL,
                     plot.type = "skip") {

    if (is(date, "POSIXct")) {
        TZ <- attr(date, "tz")
    } else {
        warning("'date' is not of class \"POSIXct\". Coerced ",
                " with time zone \"GMT\"")
        TZ <- "GMT"
        date <- as.POSIXct(date, tz = "TZ")
    }
    
    if (is.unsorted(date, strictly = TRUE))
    stop("'date' must be in strict increasing order") 
   
    restr <- FALSE
  
    if (is.null(start)) {
        start <- date[1]
    } else {
        start <- as.POSIXct(start, tz = TZ)
        restr <- TRUE
    }
    
    if (is.null(end)) {
        end <- date[length(date)]
    } else {
        end <- as.POSIXct(end, tz = TZ)
        restr <- TRUE
    }
    
    ## nevt is the number of INNER events and gives the size
    ## of the uniform sample to be used
    
    if (restr) date <- date[date > start & date < end]
    nevt <- length(date)
    
    if (nevt == 0) {
        warning("gof.date with zero evts returns NULL")
        return(NULL)
    } else if (nevt <= 5) {
        warning("gof.date called with a small number of evts (< 6)")
    }
    
    ## duration is in years and rate in events/years
    duration <- as.numeric(difftime(time1 = end, time2 = start,
                                    units = "days")) / 365.25 
    rate <- nevt / duration  
    
    ##=============================================================
    ## Perform KS test ignoring skips (if any)
    ##=============================================================
    
    ## theoretical distribution at the inner events = order
    ## statistics since F(x) = x
    F <- as.numeric(difftime(time1 = date, time2 = start,
                             units = "days")) / 365.25 / duration
    
    ## KS test from the stats package
    KS <- ks.test(x = F, y = punif)
    
    ## KS distance localisation
    F.emp <- (1:nevt) / nevt
    D <- cbind(F - c(0, F.emp[1:(nevt-1)]), F.emp - F)
    D1 <- apply(D, 1, max)
    i1 <- apply(D, 1, which.max)
    iD.max <- which.max(D1)
    D.max <- D1[iD.max]
    x.max <- date[iD.max]
    
    ## Compute F and F.emp at maximal distance for the plot.
    F.max <- F[iD.max]
    if (i1[iD.max] == 2) F.emp.max <- iD.max / nevt
    else F.emp.max <- (iD.max-1) / nevt
    
    ## cat(sprintf("x.max = %6.3f, F.emp.max = %6.3f F.max = %6.3f\n",
    ##            x.max, F.emp.max, F.max))
    
    ##=============================================================
    ## When skip is given the KS test will be computed in a special
    ## fashion 
    ##=============================================================
    
    if(!is.null(skip)) {
        
        noskip <- skip2noskip(skip = skip, start = start, end = end)
        nosn <- nrow(noskip)
        
        nosduration <-
            as.numeric(difftime(time1 = noskip$end,
                                time2 = noskip$start,
                                units = "days")) / 365.25
        effduration <- sum(nosduration) 
        
        ##============================================================
        ## Perform KS test and computations  on non-missing periods
        ## results will be output as noskip.new data.frame
        ##============================================================
        
        nosrate <- rep(NA, nosn)
        nosDn <- rep(NA, nosn)
        nosKS <- rep(NA, nosn)
        nosevt <- rep(NA, nosn)
        
        for (i in 1:nosn) {
            
            ## Fl is the leftside limit of the ecdf at xi
            indi <- (date >= noskip$start[i]) & (date <= noskip$end[i])
            ni <- sum(indi)        
            nosevt[i] <- ni
            nosrate[i] <- ni / nosduration[i]
            
            ## Here we have only INNER events
            if (ni >= 1) {
                tii <-
                    as.numeric(difftime(time1 = date[indi],
                                        time2 = noskip$start[i],
                                        units = "days")) / 365.25/ nosduration[i]
                KSi <- ks.test(x = tii, y = punif)
                nosDn[i] <- KSi$statistic
                nosKS[i] <- KSi$p.value
            }
        }
        
        ## number of inner events
        effnevt <- sum(nosevt)
        effrate <- effnevt / effduration
        
        noskip.new <- data.frame(start = I(noskip$start),
                                 end = I(noskip$end),
                                 duration = nosduration,
                                 nevt = nosevt,
                                 rate = nosrate,
                                 Dn = nosDn,
                                 KS = nosKS)
        
        ##==============================================================
        ## Now compute a modidied KS taking into account
        ## the gaps.
        ## This part is the new part in versions > 0.5-2
        ##==============================================================
        
        ie <- interevt(date, skip = skip)
        nie <- nrow(ie$interevt)
        nie1 <- nie - 1
        
        ## KS test...
        if (nie > 2) {
            
            x <- cumsum(ie$interevt$duration) / 365.25
            Fx <- x[1:nie1] / x[nie]
            
            ## KS test
            Fe <- (1:nie1) / nie1
            D <- cbind(Fx - c(0, Fe[1:(nie1-1)]), Fe - Fx)
            D1 <- apply(D, 1, max)
            i1 <- apply(D, 1, which.max)
            iD.max <- which.max(D1)
            D.max <- D1[iD.max]
            effx.max <- x[iD.max]
            
            ## Compute F and F.emp at maximal distance for the plot.
            effF.max <- Fx[iD.max]
            
            if (i1[iD.max] == 2) effF.emp.max <- iD.max / nie1
            else effF.emp.max <- (iD.max-1) / nie1
            
            ## cat(sprintf("iDn = %d x = %6.1f Dn = %6.4f\n", iDn, x[iDn], Dn))
            
            effKS <- ks.test(x = Fx, y = punif)
            
        } else {
            warning("not enough events KS is NULL")
            effKS <- NULL
        }
        
    } 
    
    ## Build the graph according to the choices
    if (plot) {
        
        if (is.null(main)) {
            main <- paste("from", format(start, "%Y-%m-%d"), "to",
                          format(end, "%Y-%m-%d"))
        }
        if ( (plot.type == "skip") || is.null(skip) ) {

            df <- rbind(data.frame(xg = start, yg = 0),
                        data.frame(xg = date, yg = F.emp),
                        data.frame(xg = end, yg = 1))
            ## xg <-  c(start, date, end)
            ## yg <-  c(0, F.emp, 1)
            
            x.left <- start
            x.right <- end
            
            plot(x = df$xg, y = df$yg, type = "n",
                 xlab = " ", ylab = " ",
                 main = main, ylim = c(0, 1))
            
            ## Rectangles
            if (!is.null(skip)){
                for (i in 1:nrow(skip)) {
                    polygon(x = c(skip$start[i], skip$end[i],
                                skip$end[i],  skip$start[i]),
                            y = c(0, 0, 1, 1),
                            col = "gray", bg = "lightgray")
                }
            }
            
            ## grid
            abline(h = pretty(c(0, 1)), col = "gray", lwd = 2, lty = "dotted")
            
            ## theroretical F
            lines(x = c(start, end), y = c(0, 1), col = "SteelBlue3", lwd = 2)
            
            ## steps
            lines(x = df$xg, y = df$yg, type = "s")
            
            ## materialize the KS distance and display results
            lines(x = c(x.max, x.max), y = c(F.max, F.emp.max),
                  col = "orangered", lwd = 2)
            
            text1 <- sprintf("p-value = %6.3f", KS$p.value)
            text2 <- sprintf("Dn = %6.3f", KS$statistic)
            
        } else {
            
            df <- data.frame(xg = c(0, cumsum(ie$interevt$duration)) / 365.25,
                             yg = c((0:nie1) / nie1, 1))
            
            x.left <- 0
            x.right <- df$xg[length(df$xg)]
            
            plot(x = df$xg, y = df$yg, type = "n",
                 xlab = "", ylab = "",
                 main = main,  ylim = c(0, 1),
                 xaxt = "n")
            
            ## axis must be done in custom style
            axis(side = 1, at = ie$axis$daysfrom/365.25,
                 labels = ie$axis$ticks)
            
            ## Flat rectangles (or segments)
            durs <- tapply(ie$interevt$duration, as.factor(ie$interevt$period), sum)
            cdurs <- c(0, cumsum(durs) / 365.25)
            
            segments(x0 = cdurs, x1 = cdurs,
                     y0 = rep(0, nie + 1), y1 =  rep(1, nie + 1),
                     col = "black")
            
            ## grid
            abline(h = pretty(c(0, 1)), col = "gray", lwd = 2, lty = "dotted")
            
            ## theoretical F      
            lines(x = c(0, sum(durs))/365.25, y = c(0, 1), col = "SteelBlue3", lwd = 2)
            
            durs <- tapply(ie$interevt$duration, as.factor(ie$interevt$period), sum)
            cdurs <- c(0, cumsum(durs) / 365.25)
            
            ## materialize the KS distance and display results
            lines(x = c(effx.max, effx.max), y = c(effF.max, effF.emp.max),
                  col = "orangered", lwd = 2)
            
            text1 <- sprintf("p-value = %6.3f", effKS$p.value)
            text2 <- sprintf("Dn = %6.3f", effKS$statistic)
            
        }
        
        ## Text annotation
        coords <- par()$usr
        
        wg <- strwidth(text1, units="user")
        hg <- strheight(text1, units="user")
        
        xgt  <- 0.25 * coords[1] + 0.75 * coords[2]
        ygt1 <- 0.8 * coords[3] + 0.2 * coords[4]
        ygt2 <- ygt1 - 1.8*hg
        
        wgt <- strwidth(text1)
        hgt <- strheight(text1)
        
        rect(xleft = xgt - wgt * 0.6, xright = xgt + wgt * 0.6,
             ytop = ygt1 + 1.5 * hgt,  ybottom = ygt2 - 1.5 * hgt,
             col = "white", border = "orange")
        
        text(x = xgt, y = ygt1, labels = text1)
        text(x = xgt, y = ygt2, labels = text2)
        
        ## These elements must apper aver all
        ## theroretical F
        lines(x = c(x.left, x.right), y = c(0, 1), col = "SteelBlue3", lwd = 2)
        
        ## steps
        lines(x = df$xg, y = df$yg, type = "s")
        
    } 
    
    if (is.null(skip)) 
        list(effKS.statistic = KS$statistic,
             effKS.pvalue = KS$p.value,
             KS.statistic = KS$statistic,
             KS.pvalue = KS$p.value,
             effnevt = nevt,
             nevt = nevt,
             rate = rate,
             effrate = rate,
             duration = duration,
             effduration = duration)
    else
        list(effKS.statistic = effKS$statistic,
             effKS.pvalue = effKS$p.value,
             KS.statistic = KS$statistic,
             KS.pvalue = KS$p.value,
             effnevt = effnevt,
             nevt = nevt,
             rate = rate,
             effrate = effrate,
             duration = duration,
             effduration = effduration,
             noskip = noskip.new)
}
