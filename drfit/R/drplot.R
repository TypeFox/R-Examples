if(getRversion() >= '2.15.1') utils::globalVariables(c("dose", "Substance", "mtype"))
drplot <- function(drresults, data, 
        dtype = "std", alpha = 0.95, ctype = "none",
        path = "./", fileprefix = "drplot", overlay = FALSE,
        xlim = c("auto","auto"), ylim = c("auto","auto"),
        xlab = paste("Decadic Logarithm of the dose in ", unit),
        ylab = "Normalized response",
        axes = TRUE, frame.plot = TRUE,
        postscript = FALSE, pdf = FALSE, png = FALSE, 
        bw = TRUE,
        pointsize = 12,
        colors = 1:8, ltys = 1:8, pchs = "auto",
        devoff=TRUE, lpos="topright")
{
    # Check if all data have the same unit
    unitlevels <- levels(as.factor(drresults$unit))
    if (length(unitlevels) == 1) {
        unit <- unitlevels
    } else {
        unit <- "different units"
    }

    # Determine the plot limits on the x-axis and y axis
    if(is.data.frame(data)) {
        # Get rid of pseudo substance names of controls
        nonzerodata <- subset(data,dose!=0)
        nonzerodata$substance <- factor(nonzerodata$substance) 
        zerodata <- subset(data,dose==0)
        nc <- length(zerodata$dose)     # Number of control points
        if (nc > 0) {
            sdc <- sd(zerodata$response)
            controlconf <- sdc * qt((1 + alpha)/2, nc - 1) / sqrt(nc)
            if (nc < 3) {
                ctype = "none"
            }
        } else {
            if (ctype != "none") {
                stop("There are no controls in the dataset, and therefore ",
                    "their scatter cannot be displayed\n")
            }
        }
        lld <- log10(min(nonzerodata$dose))
        lhd <- log10(max(nonzerodata$dose))
        hr <- max(nonzerodata$response)
        if (ctype == "std") hr <- max(hr,1 + sdc)
        if (ctype == "conf") hr <- max(hr,1 + controlconf)
        dsubstances <- levels(nonzerodata$substance)    
    } else {
        lld <- min(drresults[["logED50"]],na.rm=TRUE) - 2
        lhd <- max(drresults[["logED50"]],na.rm=TRUE) + 2
        if (length(subset(drresults,mtype=="linlogit")$Substance) != 0) {
            hr <- 1.8 
        } else {
            hr <- 1.0
        }
    }
    if (xlim[1] == "auto") xlim[1] <- lld - 0.5
    if (xlim[2] == "auto") xlim[2] <- lhd + 1
    if (ylim[1] == "auto") ylim[1] <- -0.1
    if (ylim[2] == "auto") ylim[2] <- hr + 0.2
    xlim <- as.numeric(xlim)
    ylim <- as.numeric(ylim)

    # Prepare overlay plot if requested
    if (overlay)
    {
        if (pchs == "auto") pchs = 1:length(dsubstances)
        if (postscript) {
            filename = paste(path,fileprefix,".eps",sep="")
            postscript(file=filename,
                    paper="special",width=7,height=7,horizontal=FALSE, pointsize=pointsize)
            message("Created File: ",filename,"\n")
        } 
        if (pdf) {
            filename = paste(path,fileprefix,".pdf",sep="")
            pdf(file=filename,
                    paper="special",width=7,height=7,horizontal=FALSE, pointsize=pointsize)
            message("Created File: ",filename,"\n")
        } 
        if (png) {
            filename = paste(path,fileprefix,".png",sep="")
            png(filename=filename,
                width=500, height=500, pointsize=pointsize)
            message("Created File: ",filename,"\n")
        }
            
        plot(0,type="n",
            xlim = xlim,
            ylim = ylim,
            xlab = xlab,
            ylab = ylab,
            axes = axes,
            frame.plot = frame.plot)
    } else {
        # If overlay plot is not requested, ask before showing multiple plots on the screen
        if (pchs == "auto") pchs = rep(1, length(dsubstances))
        if (!postscript && !png && !pdf && length(dsubstances) > 1) {
            op <- par(ask=TRUE)
            on.exit(par(op))
        } 
    }
    # nl is the overall number of fits to draw by different line types
    nl <- 0

    # Plot the data either as raw data or as error bars
    if(is.data.frame(data)) {
        splitted <- split(nonzerodata,nonzerodata$substance)
        # n is the index for the dose-response curves
        n <- 0
        if (bw) colors <- rep("black",length(dsubstances))
        # Loop over the substances in the data (index n)
        for (i in dsubstances) {
            n <- n + 1
            tmp <- splitted[[i]]
            if (length(tmp$response) != 0) {
                color <- colors[[n]]
                pch <- pchs[[n]]
                # Prepare the single graphs if an overlay is not requested
                if (!overlay)
                {
                    if (postscript) {
                        filename = paste(path,fileprefix,sub(" ","_",i),".eps",sep="")
                        postscript(file=filename,
                                paper="special",width=7,height=7,horizontal=FALSE,pointsize=pointsize)
                        message("Created File: ",filename,"\n")
                    } 
                    if (pdf) {
                        filename = paste(path,fileprefix,sub(" ","_",i),".pdf",sep="")
                        pdf(file=filename,
                                paper="special",width=7,height=7,horizontal=FALSE,pointsize=pointsize)
                        message("Created File: ",filename,"\n")
                    } 
                    if (png) {
                        filename = paste(path,fileprefix,sub(" ","_",i),".png",sep="")
                        png(filename=filename,
                            width=500, height=500, pointsize=pointsize)
                        message("Created File: ",filename,"\n")
                    }
                        
                    plot(0,type="n",
                        xlim = xlim,
                        ylim = ylim,
                        xlab = xlab,
                        ylab = ylab,
                        axes = axes,
                        frame.plot = frame.plot)
                }
                if (!overlay) legend(lpos, i, lty = 1, col = color, pch = pch, inset=0.05)
                tmp$dosefactor <- factor(tmp$dose)  # necessary because the old
                                                    # factor has all levels, not 
                                                    # only the ones tested with
                                                    # this substance

                # Plot the control lines, if requested
                if (ctype == "std") {
                    abline(h = 1 - sdc, lty = 2)
                    abline(h = 1 + sdc, lty = 2)
                }
                if (ctype == "conf") {
                    abline(h = 1 - controlconf, lty = 2) 
                    abline(h = 1 + controlconf, lty = 2)
                }

                # Plot the data, if requested
                if (dtype != "none") {
                    if (dtype == "raw") {
                        points(log10(tmp$dose), tmp$response, col = color, pch = pch)
                    } else {
                        splitresponses <- split(tmp$response,tmp$dosefactor)
                        means <- sapply(splitresponses,mean)
                        lengths <- sapply(splitresponses,length)
                        vars <- sapply(splitresponses,var)
                        standarddeviations <- sqrt(vars)
                    }
                    if (dtype == "std")
                    {
                        tops <- means + standarddeviations
                        bottoms <- means - standarddeviations
                    }
                    if (dtype == "conf")
                    {
                        confidencedeltas <- qt((1 + alpha)/2, lengths - 1) * sqrt(vars) 
                        tops <- means + confidencedeltas
                        bottoms <- means - confidencedeltas
                    }
                    if (dtype != "raw")
                    {
                        x <- log10(as.numeric(levels(tmp$dosefactor)))
                        segments(x,bottoms,x,tops,col=color)
                        points(x, means, col = color, pch = pch)
                        smidge <- 0.05
                        segments(x - smidge,bottoms,x + smidge,bottoms,col=color)
                        segments(x - smidge,tops,x + smidge,tops,col=color)
                    }
                }

                # Plot the fits for this substance, if there are any
                fits <- subset(drresults,Substance == i)
                fit.rows <- rownames(fits)
                nf <-  length(fit.rows)  # number of fits to plot for this substance
                if (nf > 0) {
                    for (j in 1:nf)
                    {
                        fit.row = fit.rows[j]
                        if (overlay) nl <- nl + 1 else nl = j
                        lty <- ltys[nl]
                        if (drresults[[fit.row, "mtype"]] %in% c("probit",
                                                                 "logit",
                                                                 "weibull",
                                                                 "linlogit")) 
                        {
                            m <- attr(drresults, "models")[[as.numeric(fit.row)]]
                            of <- function(x) {
                                predict(m, data.frame(dose = 10^x))
                            } 
                            plot(of, lld - 0.5, lhd + 2, 
                                 add = TRUE, col = color, lty = lty)
                        }
                    }
                }
                if (!overlay && (postscript || png || pdf)) dev.off()
            } else {
                message("No data for ",i,"\n")
            }
        }
    }
    if (overlay) legend(lpos, dsubstances, col = colors, pch = pchs, lty = ltys, inset=0.05)
    if (overlay && (postscript || png || pdf)) {
        if (devoff) {
            dev.off()
        }
    }
}
