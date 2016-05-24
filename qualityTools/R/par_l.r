paretoChart = function(x, weight, showTable = TRUE, las = 0, main, col, border, xlab, ylab = "Frequency", percentVec, ...) {
    varName = deparse(substitute(x))[1]
    corp.col = "pink3"
    corp.border = "red3"
    if (!is.vector(x) & !is.data.frame(x) & !is.table(x)) 
        stop("x should be a vector, dataframe or a table")
    if (is.table(x)) {
        xtable = x
    }
    if (is.vector(x)) {
        if (!is.null(names(x))) 
            xtable = as.table(x)
        else xtable = table(x)
    }
    if (!missing(weight)) {
        if (!is.numeric(weight)) 
            stop("weight must be numeric!")
        if (is.null(names(weight))) 
            stop("weight is missing names for matching!")
        else {
            if (FALSE %in% (sort(names(weight)) == sort(names(xtable)))) 
                stop("names of weight and table do not match!")
            else {
                for (i in 1:length(xtable)) {
                  xtable[i] = weight[names(weight) == names(xtable)[i]] * xtable[i]
                }
            }
        }
    }
    else {
        weight = FALSE
    }
    if (missing(showTable)) 
        showTable = TRUE
    if (missing(xlab)) 
        xlab = ""
    if (missing(main)) 
        main = paste("Pareto Chart for ", varName)
    if (missing(col)) 
        col = corp.col
    if (missing(border)) 
        border = corp.border
    if (missing(percentVec)) 
        percentVec = seq(0, 1, by = 0.25)
    call <- match.call(expand.dots = TRUE)
    if (length(xtable) > 1) {
        ylim = c(min(xtable), max(xtable) * 1.025)
        xtable = c(sort(xtable, decreasing = TRUE, na.last = TRUE))
        cumFreq = cumsum(xtable)
        sumFreq = sum(xtable)
        percentage = xtable/sum(xtable) * 100
        cumPerc = cumFreq/sumFreq * 100
        frameOut = data.frame(Frequency = as.numeric(xtable), Cum.Freq = cumFreq, Percentage = percentage, Cum.Perc = cumPerc)
        names(frameOut) = c(ylab, paste("Cum.", ylab), "Percentage", "Cum. Percentage")
        row.names(frameOut) = names(xtable)
        frameInt = as.data.frame(t(frameOut))
        names(frameInt) = rep(" ", dim(frameInt)[2])
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        lineHeight = par("csi") * par("lheight") * par("cex")
        tablespace = (2 * 4) * (lineHeight/par("fin")[2])
        plot.new()
        if (las == 0 | las == 1) {
            mymai = par("mai")
            mymai[1] = max(strheight(names(xtable), units = "inches")) * 3
            mymai[4] = strheight("Cumulative Percentage", units = "inches") * 8
            mymai[2] = max(strwidth(names(frameOut), units = "inches")) * 1.2
            par(mai = mymai, new = TRUE)
        }
        if (las == 2 | las == 3) {
            mymai = par("mai")
            mymai[1] = max(strwidth(names(xtable), units = "inches")) * 1.4
            mymai[4] = strheight("Cumulative Percentage", units = "inches") * 8
            mymai[2] = max(strwidth(names(frameOut), units = "inches")) * 1.2
            par(mai = mymai, new = TRUE)
        }
        if (showTable) {
            par(fig = c(0, 1, tablespace, 1))
            xValue = barplot(xtable, axes = FALSE, las = las, width = 1, space = 0.2, xlim = c(0.2, 1.2 * length(xtable)), main = main, ylim = c(0, sum(xtable) + 
                0.01 * (sum(xtable))), ylab = ylab, xlab = xlab, col = col, border = border)
            axis(1, at = xValue, labels = names(xtable), las = las)
            axis(2)
            axis(4, at = percentVec * (sumFreq), labels = percentVec)
            mtext(4, text = "Cumulative Percentage", line = 3)
            lines(xValue, cumFreq, col = corp.border)
            points(xValue, cumFreq, col = corp.border, pch = 15)
            par(fig = c(0, 1, 0, tablespace), new = TRUE)
            mymai[1] = 0
            mymai[3] = 0
            par(mai = mymai)
            plot(xValue, rep(1, length(xValue)), xlim = c(0.2, 1.2 * length(xtable)), ylim = c(0, 5), axes = FALSE, ylab = "", type = "n")
            axis(2, pos = 0.2, at = 1:4, labels = rev(c(ylab, paste("Cum.", ylab), "Percentage", "Cum. Percentage")), tick = FALSE, las = 1)
            numCol = dim(frameInt)[2]
            numRow = dim(frameInt)[1]
            for (i in 1:numCol) {
                for (j in 1:numRow) {
                  text(xValue[i], numRow + 1 - j, round(frameInt[j, i]), adj = c(1, 0.5))
                }
            }
        }
        else {
            mymai[2] = mymai[4]
            par(mai = mymai)
            xValue = barplot(xtable, axes = FALSE, las = las, width = 1, space = 0.2, xlim = c(0.2, 1.2 * length(xtable)), main = main, ylim = c(0, sum(xtable) + 
                0.01 * (sum(xtable))), ylab = ylab, xlab = xlab, col = col, border = border)
            axis(1, at = xValue, labels = names(xtable), las = las)
            axis(2)
            axis(4, at = percentVec * (sumFreq), labels = percentVec)
            mtext(4, text = "Cumulative Percentage", line = 3)
            lines(xValue, cumFreq, col = corp.border)
            points(xValue, cumFreq, col = corp.border, pch = 15)
        }
    }
    else {
        warning("data should have at least two categories!")
    }
    frameOut = frameInt
    for (i in 3:nrow(frameInt)) {
        frameInt[i, ] = sprintf("%.1f%%", frameInt[i, ])
    }
#    cat(paste("Pareto Analysis for", varName, "\n"))
#    cat("---\n")
    print(format(frameInt, digits = 3))
    cat("\n")
    frameOut
} 
