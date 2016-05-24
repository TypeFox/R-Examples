siarproportionbygroupplot <-
function (siardata, siarversion = 0, probs = c(95, 75, 50), xlabels = NULL, 
    grp = NULL, type = "boxes", clr = gray((9:1)/10), scl = 1, 
    xspc = 0.5, prn = FALSE, leg = FALSE) 
{   

  #attach(siardata)
    if (siardata$SHOULDRUN == FALSE && siardata$GRAPHSONLY == 
        FALSE) {
        cat("You must load in some data first (via option 1) in order to use \n")
        cat("this feature of the program. \n")
        cat("Press <Enter> to continue")
        readline()
        invisible()
        return(NULL)
    }
    if (length(siardata$output) == 0) {
        cat("No output found - check that you have run the SIAR model. \n \n")
        return(NULL)
    }
    
    cat("Plot of proportions by group \n")
    cat("Producing plot..... \n \n")
    groupnames <- as.character(1:siardata$numgroups)
    
    if (is.null(grp)) {
        cat("Enter the group number you wish to plot \n")
        cat("The choices are:\n")
        title <- "The available options are:"
        choose2 <- menu(groupnames)
    }
    else {
        choose2 <- grp
    }

    groupseq <- seq(1, siardata$numsources, by = 1)
    shift <- siardata$numsources + siardata$numiso
    usepars <- siardata$output[, ((choose2-1)*(shift) + 1) : ((choose2-1)*(shift) + shift - siardata$numiso)]
    newgraphwindow()
    if (siardata$TITLE != "SIAR data") {
        plot(1, 1, xlab = "Source", ylab = "Proportion", main = paste(siardata$TITLE, 
            " by group: ", groupnames[choose2], sep = ""), 
            xlim = c(min(groupseq) - xspc, max(groupseq) + xspc), 
            ylim = c(0, 1), type = "n", xaxt = "n")
        if (is.null(xlabels)) {
            axis(side = 1, at = min(groupseq):max(groupseq), 
                labels = (as.character(siardata$sources[,1])))
        }
        else {
            axis(side = 1, at = min(groupseq):max(groupseq), 
                labels = (xlabels))
        }
    }
    else {
        plot(1, 1, xlab = "Source", ylab = "Proportion", main = paste("Proportions by group: ", 
            groupnames[choose2], sep = ""), xlim = c(min(groupseq) - 
            xspc, max(groupseq) + xspc), ylim = c(0, 1), type = "n", 
            xaxt = "n")
        if (is.null(xlabels)) {
            axis(side = 1, at = min(groupseq):max(groupseq), 
                labels = (as.character(siardata$sources[,1])))
        }
        else {
            axis(side = 1, at = min(groupseq):max(groupseq), 
                labels = (xlabels))
        }
    }
    if (siarversion > 0) 
        mtext(paste("siar v", siarversion), side = 1, line = 4, 
            adj = 1, cex = 0.6)
    clrs <- rep(clr, 5)
    for (j in 1:ncol(usepars)) {
        temp <- hdr(usepars[, j], probs, h = bw.nrd0(usepars[, 
            j]))$hdr
        line_widths <- seq(2, 20, by = 4) * scl
        bwd <- c(0.1, 0.15, 0.2, 0.25, 0.3) * scl
        if (prn == TRUE) {
            cat(paste("Probability values for Group", j, "\n"))
        }
        for (k in 1:length(probs)) {
            temp2 <- temp[k, ]
            if (type == "boxes") {
                polygon(c(groupseq[j] - bwd[k], groupseq[j] - 
                  bwd[k], groupseq[j] + bwd[k], groupseq[j] + 
                  bwd[k]), c(max(min(temp2[!is.na(temp2)]), 0), 
                  min(max(temp2[!is.na(temp2)]), 1), min(max(temp2[!is.na(temp2)]), 
                    1), max(min(temp2[!is.na(temp2)]), 0)), col = clrs[k])
            }
            if (type == "lines") {
                lines(c(groupseq[j], groupseq[j]), c(max(min(temp2[!is.na(temp2)]), 
                  0), min(max(temp2[!is.na(temp2)]), 1)), lwd = line_widths[k], 
                  lend = 2)
            }
            if (prn == TRUE) {
                cat(paste("\t", probs[k], "% lower =", format(max(min(temp2[!is.na(temp2)]), 
                  0), digits = 2, scientific = FALSE), "upper =", 
                  format(min(max(temp2[!is.na(temp2)]), 1), digits = 2, 
                    scientific = FALSE), "\n"))
            }
        }
    }
    if (leg == TRUE) {
        if (type == "lines") {
            legnames <- character(length = length(probs))
            for (i in 1:length(probs)) {
                legnames[i] <- paste(probs[i], "%", sep = "")
            }
            legend(mean(c(min(groupseq), max(groupseq))), 1.02, 
                legend = legnames, lwd = c(2, 6, 10), ncol = length(probs), 
                xjust = 0.5, text.width = strwidth(legnames)/2, 
                bty = "n")
        }
        if (type == "boxes") {
            print("Legends not yet supported for box style graph. Use type=lines with leg=TRUE instead.")
        }
    }
    cat("Please maximise this graph before saving or printing. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
}
