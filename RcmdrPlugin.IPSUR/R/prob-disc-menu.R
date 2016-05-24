# Last modified Feb 16, 2008


`binomialDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Binomial Distribution"))
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    functionVar <- tclVar("Probability")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Probability")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        trials <- tclvalue(trialsVar)
        if (trials == "") {
            errorCondition(recall = binomialDistributionPlot.ipsur, 
                message = gettextRcmdr("Number of trials must be specified."))
            return()
        }
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = binomialDistributionPlot.ipsur, 
                message = gettextRcmdr("Probability of success must be specified."))
            return()
        }
        fun <- tclvalue(functionVar)
        command <- paste("xmin <- qbinom(.0005, size=", trials, " , prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qbinom(.9995, size=", trials, " , prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".x <- xmin:xmax", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        if (fun == "Probability") {
            doItAndPrint(paste("plot(.x, dbinom(.x, size=", trials, 
                ", prob=", prob, "), xlab=\"Number of Successes\", ylab=\"Probability Mass\", main=\"Binomial Distribution: Trials = ", 
                trials, ", Probability of success = ", prob, 
                "\", type=\"h\")", sep = ""))
            doItAndPrint(paste("points(.x, dbinom(.x, size=", 
                trials, ", prob=", prob, "), pch=16)", sep = ""))
        }
        else if (fun == "Cumulative Probability") {
            justDoIt(paste("plot( stepfun(.x, pbinom((xmin-1):xmax", 
                ", size=", trials, ", prob=", prob, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of Successes\", ylab=\"Cumulative Probability\", main=\"Binomial Distribution: Trials = ", 
                trials, ", Probability of success = ", prob, 
                "\")", sep = ""))
            logger(paste("plot( stepfun(.x, pbinom((xmin-1):xmax", 
                ", size=", trials, ", prob=", prob, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of Successes\", ylab=\"Cumulative Probability\", main=\"Binomial Distribution: Trials = ", 
                trials, ", Probability of success = ", prob, 
                "\")", sep = ""))
            doItAndPrint(paste("points( .x, pbinom(xmin:xmax,", 
                " size=", trials, ", prob=", prob, "), pch = 16, cex=1.2 )", 
                sep = ""))
            doItAndPrint(paste("points( .x, pbinom((xmin-1):(xmax-1),", 
                " size=", trials, ", prob=", prob, "), pch = 1, cex=1.2 )", 
                sep = ""))
            doItAndPrint("abline( h = 1, lty = 2, col = \"grey\" )")
        }
        else if (fun == "Quantile Function") {
            doItAndPrint(paste("q <- pbinom((xmin-1):xmax, size =", 
                trials, ",prob = ", prob, ")", sep = ""))
            main <- paste("\"Binomial QF\n size=", trials, ", prob=", 
                prob, "\"", sep = "")
            ylab <- "\"Number of Successes\""
            xlab <- "\"Cumulative Probability\""
            command <- paste("plot(q[1:(length(q)-1)], .x, xlim=c(0,1), main =", 
                main, ", xlab = ", xlab, ",ylab = ", ylab, ", cex = 1.2)", 
                sep = "")
            justDoIt(command)
            logger(command)
            doItAndPrint(paste("abline( v = 0, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("abline( v = 1, lty = 2, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("points(q[2:(length(q))], .x, pch=16, cex=1.2)", 
                sep = ""))
            doItAndPrint(paste("segments(  q[1:(length(q)-1)], .x, q[2:(length(q))], .x)",
                sep = ""))
        }
        doItAndPrint("abline( h = 0, lty = 2, col = \"grey\" )")
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("size (number of trials)")), 
        trialsEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot probability mass function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(trialsEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = trialsEntry)
}



`binomialMass.ipsur` <-
function () 
{
    checkTrials <- function(trials) {
        RcmdrTkmessageBox(message = sprintf(gettextRcmdr("Number of trials, %d, is large.\nCreate long output?"), 
            trials), icon = "warning", type = "yesno", default = "no")
    }
    initializeDialog(title = gettextRcmdr("Binomial Probabilities"))
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    onOK <- function() {
        closeDialog()
        trials <- tclvalue(trialsVar)
        prob <- tclvalue(probVar)
        if (trials == "") {
            errorCondition(recall = binomialMass.ipsur, message = gettextRcmdr("Number of trials not specified."))
            return()
        }
        if (prob == "") {
            errorCondition(recall = binomialMass.ipsur, message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        command <- paste("xmin <- qbinom(.00005, size=", trials, " , prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qbinom(.99995, size=", trials, " , prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".x <- xmin:xmax", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".Table <- data.frame(Pr=dbinom(xmin:xmax, size=", 
            trials, ", prob=", prob, "))", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        logger(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        justDoIt(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        doItAndPrint(".Table")
        logger("remove(.Table, xmin, xmax)")
        remove(.Table, xmin, xmax, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("size (number of trials)")), 
        trialsEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success)")), 
        probEntry, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(trialsEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    dialogSuffix(rows = 3, columns = 2, focus = trialsEntry)
}



`binomialProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Cumulative Binomial Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        trials <- tclvalue(trialsVar)
        prob <- tclvalue(probVar)
        if ("" == probabilities) {
            errorCondition(recall = binomialProbabilities.ipsur, 
                message = gettextRcmdr("Values not specified."))
            return()
        }
        if (trials == "") {
            errorCondition(recall = binomialProbabilities.ipsur, 
                message = gettextRcmdr("Number of trials not specified."))
            return()
        }
        if (prob == "") {
            errorCondition(recall = binomialProbabilities.ipsur, 
                message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("pbinom(c(", probabilities, "), size=", 
            trials, ", prob=", prob, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("size (number of trials)")), 
        trialsEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(trialsEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = probabilitiesEntry)
}



`binomialQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Binomial Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        trials <- tclvalue(trialsVar)
        prob <- tclvalue(probVar)
        if ("" == quantiles) {
            errorCondition(recall = binomialQuantiles.ipsur, 
                message = gettextRcmdr("Probabilities not specified."))
            return()
        }
        if (trials == "") {
            errorCondition(recall = binomialQuantiles.ipsur, 
                message = gettextRcmdr("Number of trials not specified."))
            return()
        }
        if (prob == "") {
            errorCondition(recall = binomialQuantiles.ipsur, 
                message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qbinom(c(", quantiles, "), size=", 
            trials, ", prob=", prob, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("size (number of trials)")), 
        trialsEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(trialsEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}



`geomDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Geometric Distribution"))
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    functionVar <- tclVar("Probability")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Probability")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = geomDistributionPlot.ipsur, 
                message = gettextRcmdr("Success probability was not specified."))
            return()
        }
        fun <- tclvalue(functionVar)
        command <- paste("xmin <- qgeom(.0005, prob=", prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qgeom(.9995, prob=", prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".x <- xmin:xmax", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        if (fun == "Probability") {
            doItAndPrint(paste("plot(.x, dgeom(.x, prob=", prob, 
                "), xlab=\"Number of Failures until Success\", ylab=\"Probability Mass\", main=\"Geometric Distribution: Prob of success = ", 
                prob, "\", type=\"h\")", sep = ""))
            doItAndPrint(paste("points(.x, dgeom(.x, prob=", 
                prob, "), pch=16)", sep = ""))
        }
        else if (fun == "Cumulative Probability") {
            justDoIt(paste("plot( stepfun(.x, pgeom((xmin-1):xmax", 
                ", prob=", prob, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of Failures until Success\", ylab=\"Cumulative Probability\", main=\"Geometric Distribution: Prob of success = ", 
                prob, "\")", sep = ""))
            logger(paste("plot( stepfun(.x, pgeom((xmin-1):xmax", 
                ", prob=", prob, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of Failures until Success\", ylab=\"Cumulative Probability\", main=\"Geometric Distribution: Prob of success = ", 
                prob, "\")", sep = ""))
            doItAndPrint(paste("points( .x, pgeom(xmin:xmax,", 
                ", prob=", prob, "), pch = 16, cex=1.2 )", sep = ""))
            doItAndPrint(paste("points( .x, pgeom((xmin-1):(xmax-1),", 
                ", prob=", prob, "), pch = 1, cex=1.2 )", sep = ""))
            doItAndPrint("abline( h = 1, lty = 2, col = \"grey\" )")
        }
        else if (fun == "Quantile Function") {
            doItAndPrint(paste("q <- pgeom((xmin-1):xmax, prob =", 
                prob, ")", sep = ""))
            main <- paste("\"Geometric QF: prob=", prob, "\"", 
                sep = "")
            ylab <- "\"Number of Failures until First Success\""
            xlab <- "\"Cumulative Probability\""
            command <- paste("plot(q[1:(length(q)-1)], .x, xlim=c(0,1), main =", 
                main, ", xlab = ", xlab, ",ylab = ", ylab, ", cex = 1.2)", 
                sep = "")
            justDoIt(command)
            logger(command)
            doItAndPrint(paste("abline( v = 0, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("abline( v = 1, lty = 2, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("points(q[2:(length(q))], .x, pch=16, cex=1.2)", 
                sep = ""))
            doItAndPrint(paste("segments(  q[1:(length(q)-1)], .x, q[2:(length(q))], .x)",
                sep = ""))
        }
        doItAndPrint("abline( h = 0, lty = 2, col = \"grey\" )")
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dgeom")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot probability mass function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = probEntry)
}



`geomMass.ipsur` <-
function () 
{
    checkRange <- function(range) {
        RcmdrTkmessageBox(message = sprintf(gettextRcmdr("Range of values over which to plot, %d, is large.\nCreate long output?"), 
            range), icon = "warning", type = "yesno", default = "no")
    }
    initializeDialog(title = gettextRcmdr("Geometric Probabilities"))
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    onOK <- function() {
        closeDialog()
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = geomMass.ipsur, message = gettextRcmdr("Success probability was not specified."))
            return()
        }
        command <- paste("xmin <- qgeom(.00005, prob=", prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qgeom(.99995, prob=", prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".Table <- data.frame(Pr=dgeom(xmin:xmax, prob=", 
            prob, "))", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        logger(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        justDoIt(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        doItAndPrint(".Table")
        logger("remove(.Table, xmin, xmax)")
        remove(.Table, xmin, xmax, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dgeom")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    dialogSuffix(rows = 2, columns = 2, focus = probEntry)
}



`geomProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Geometric Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = geomProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        prob <- tclvalue(probVar)
        tail <- tclvalue(tailVar)
        if (prob == "") {
            errorCondition(recall = geomProbabilities.ipsur, 
                message = gettextRcmdr("Success probability was not specified."))
            return()
        }
        doItAndPrint(paste("pgeom(c(", probabilities, "), prob=", 
            prob, ", lower.tail=", tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pgeom")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}



`geomQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Geometric Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = geomQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        prob <- tclvalue(probVar)
        tail <- tclvalue(tailVar)
        if (prob == "") {
            errorCondition(recall = geomQuantiles.ipsur, message = gettextRcmdr("Success probability not specified."))
            return()
        }
        doItAndPrint(paste("qgeom(c(", quantiles, "), prob=", 
            prob, ", lower.tail=", tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qgeom")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}



`hyperDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Hypergeometric Distribution"))
    mVar <- tclVar("1")
    mEntry <- tkentry(top, width = "6", textvariable = mVar)
    nVar <- tclVar("1")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    kVar <- tclVar("1")
    kEntry <- tkentry(top, width = "6", textvariable = kVar)
    functionVar <- tclVar("Probability")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Probability")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        m <- tclvalue(mVar)
        n <- tclvalue(nVar)
        k <- tclvalue(kVar)
        fun <- tclvalue(functionVar)
        if (m == "") {
            errorCondition(recall = hyperDistributionPlot.ipsur, 
                message = gettextRcmdr("The m parameter was not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = hyperDistributionPlot.ipsur, 
                message = gettextRcmdr("The n parameter was not specified."))
            return()
        }
        if (k == "") {
            errorCondition(recall = hyperDistributionPlot.ipsur, 
                message = gettextRcmdr("The k parameter was not specified."))
            return()
        }
        command <- paste("xmin <- qhyper(.0005, m=", m, ", n=", n, ", k=", 
            k, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qhyper(.9995, m=", m, ", n=", n, ", k=", 
            k, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".x <- xmin:xmax", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        if (fun == "Probability") {
            doItAndPrint(paste("plot(.x, dhyper(.x, m=", m, ", n=", 
                n, ", k=", k, "), xlab=\"Number of White Balls in Sample\", ylab=\"Probability Mass\", main=\"Hypergeometric Distribution: m=", 
                m, ", n=", n, ", k=", k, "\", type=\"h\")", sep = ""))
            doItAndPrint(paste("points(.x, dhyper(.x, m=", m, 
                ", n=", n, ", k=", k, "), pch=16)", sep = ""))
        }
        else if (fun == "Cumulative Probability") {
            justDoIt(paste("plot( stepfun(.x, phyper((xmin-1):xmax", 
                ", m=", m, ", n=", n, ", k=", k, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of White Balls in Sample\", ylab=\"Cumulative Probability\", main=\"Hypergeometric Distribution: m=", 
                m, ", n=", n, ", k=", k, "\")", sep = ""))
            logger(paste("plot( stepfun(.x, phyper((xmin-1):xmax", 
                ", m=", m, ", n=", n, ", k=", k, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of White Balls in Sample\", ylab=\"Cumulative Probability\", main=\"Hypergeometric Distribution: m=", 
                m, ", n=", n, ", k=", k, "\")", sep = ""))
            doItAndPrint(paste("points( .x, phyper(xmin:xmax,", 
                ", m=", m, ", n=", n, ", k=", k, "), pch = 16, cex=1.2 )", 
                sep = ""))
            doItAndPrint(paste("points( .x, phyper((xmin-1):(xmax-1),", 
                ", m=", m, ", n=", n, ", k=", k, "), pch = 1, cex=1.2 )", 
                sep = ""))
            doItAndPrint("abline( h = 1, lty = 2, col = \"grey\" )")
        }
        else if (fun == "Quantile Function") {
            doItAndPrint(paste("q <- phyper((xmin-1):xmax, m=", 
                m, ", n=", n, ", k=", k, ")", sep = ""))
            main <- paste("\"Hypergeometric QF\n m=", m, ", n=", 
                n, ", k=", k, "\"", sep = "")
            ylab <- "\"Number of White Balls in Sample\""
            xlab <- "\"Cumulative Probability\""
            command <- paste("plot(q[1:(length(q)-1)], .x, xlim=c(0,1), main =", 
                main, ", xlab = ", xlab, ",ylab = ", ylab, ", cex = 1.2)", 
                sep = "")
            justDoIt(command)
            logger(command)
            doItAndPrint(paste("abline( v = 0, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("abline( v = 1, lty = 2, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("points(q[2:(length(q))], .x, pch=16, cex=1.2)", 
                sep = ""))
            doItAndPrint(paste("segments(  q[1:(length(q)-1)], .x, q[2:(length(q))], .x)",
                sep = ""))
        }
        doItAndPrint("abline( h = 0, lty = 2, col = \"grey\" )")
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dhyper")
    tkgrid(tklabel(top, text = gettextRcmdr("m (num of white balls in the urn)")), 
        mEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("n (num of black balls in the urn)")), 
        nEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("k (num of balls drawn from the urn)")), 
        kEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot probability mass function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(mEntry, sticky = "w")
    tkgrid.configure(nEntry, sticky = "w")
    tkgrid.configure(kEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = mEntry)
}



`hyperMass.ipsur` <-
function () 
{
    checkRange <- function(range) {
        RcmdrTkmessageBox(message = sprintf(gettextRcmdr("Range of values over which to plot, %d, is large.\nCreate long output?"), 
            range), icon = "warning", type = "yesno", default = "no")
    }
    initializeDialog(title = gettextRcmdr("Hypergeometric  Probabilities"))
    mVar <- tclVar("1")
    mEntry <- tkentry(top, width = "6", textvariable = mVar)
    nVar <- tclVar("1")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    kVar <- tclVar("1")
    kEntry <- tkentry(top, width = "6", textvariable = kVar)
    onOK <- function() {
        closeDialog()
        m <- tclvalue(mVar)
        n <- tclvalue(nVar)
        k <- tclvalue(kVar)
        if (m == "") {
            errorCondition(recall = hyperMass.ipsur, message = gettextRcmdr("The m parameter was not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = hyperMass.ipsur, message = gettextRcmdr("The n parameter was not specified."))
            return()
        }
        if (k == "") {
            errorCondition(recall = hyperMass.ipsur, message = gettextRcmdr("The k parameter was not specified."))
            return()
        }
        command <- paste("xmin <- qhyper(.00005, m=", m, ", n=", n, ", k=", 
            k, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qhyper(.99995, m=", m, ", n=", n, ", k=", 
            k, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".Table <- data.frame(Pr=dhyper(xmin:xmax, m=", 
            m, ", n=", n, ", k=", k, "))", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        logger(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        justDoIt(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        doItAndPrint(".Table")
        logger("remove(.Table, xmin, xmax)")
        remove(.Table, xmin, xmax, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dhyper")
    tkgrid(tklabel(top, text = gettextRcmdr("m (num of white balls in the urn)")), 
        mEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("n (num of black balls in the urn)")), 
        nEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("k (num of balls drawn from the urn)")), 
        kEntry, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(mEntry, sticky = "w")
    tkgrid.configure(nEntry, sticky = "w")
    tkgrid.configure(kEntry, sticky = "w")
    dialogSuffix(rows = 2, columns = 2, focus = mEntry)
}



`hyperProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Hypergeometric Probabilities"))
    ProbabilitiesVar <- tclVar("")
    ProbabilitiesEntry <- tkentry(top, width = "30", textvariable = ProbabilitiesVar)
    mVar <- tclVar("1")
    mEntry <- tkentry(top, width = "6", textvariable = mVar)
    nVar <- tclVar("1")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    kVar <- tclVar("1")
    kEntry <- tkentry(top, width = "6", textvariable = kVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(ProbabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = hyperProbabilities.ipsur, 
                message = gettextRcmdr("No probabilities specified."))
            return()
        }
        m <- tclvalue(mVar)
        n <- tclvalue(nVar)
        k <- tclvalue(kVar)
        if (m == "") {
            errorCondition(recall = hyperProbabilities.ipsur, 
                message = gettextRcmdr("The m parameter was not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = hyperProbabilities.ipsur, 
                message = gettextRcmdr("The n parameter was not specified."))
            return()
        }
        if (k == "") {
            errorCondition(recall = hyperProbabilities.ipsur, 
                message = gettextRcmdr("The k parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("phyper(c(", probabilities, "), m=", 
            m, ", n=", n, ", k=", k, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "phyper")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        ProbabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("m (num of white balls in the urn)")), 
        mEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("n (num of black balls in the urn)")), 
        nEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("k (num of balls drawn from the urn)")), 
        kEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(ProbabilitiesEntry, sticky = "w")
    tkgrid.configure(mEntry, sticky = "w")
    tkgrid.configure(nEntry, sticky = "w")
    tkgrid.configure(kEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = ProbabilitiesEntry)
}



`hyperQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Hypergeometric Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    mVar <- tclVar("1")
    mEntry <- tkentry(top, width = "6", textvariable = mVar)
    nVar <- tclVar("1")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    kVar <- tclVar("1")
    kEntry <- tkentry(top, width = "6", textvariable = kVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = hyperQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        m <- tclvalue(mVar)
        n <- tclvalue(nVar)
        k <- tclvalue(kVar)
        if (m == "") {
            errorCondition(recall = hyperQuantiles.ipsur, message = gettextRcmdr("The m parameter was not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = hyperQuantiles.ipsur, message = gettextRcmdr("The n parameter was not specified."))
            return()
        }
        if (k == "") {
            errorCondition(recall = hyperQuantiles.ipsur, message = gettextRcmdr("The k parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qhyper(c(", quantiles, "), m=", m, 
            ", n=", n, ", k=", k, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qhyper")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("m (num of white balls in the urn)")), 
        mEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("n (num of black balls in the urn) ")), 
        nEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("k (num of balls drawn from the urn)")), 
        kEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(mEntry, sticky = "w")
    tkgrid.configure(nEntry, sticky = "w")
    tkgrid.configure(kEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}



`negbinomialDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Negative Binomial Distribution"))
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    functionVar <- tclVar("Probability")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Probability")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        trials <- tclvalue(trialsVar)
        if (trials == "") {
            errorCondition(recall = negbinomialDistributionPlot.ipsur, 
                message = gettextRcmdr("Target number of successes was not specified."))
            return()
        }
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = negbinomialDistributionPlot.ipsur, 
                message = gettextRcmdr("Probability of success was not specified."))
            return()
        }
        fun <- tclvalue(functionVar)
        command <- paste("xmin <- qnbinom(.0005, size=", trials, ", prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qnbinom(.9995, size=", trials, ", prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".x <- xmin:xmax", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        if (fun == "Probability") {
            doItAndPrint(paste("plot(.x, dnbinom(.x, size=", 
                trials, ", prob=", prob, "), xlab=\"Number of Failures until Target Successes\", ylab=\"Probability Mass\", main=\"Negative Binomial Distribution: Target = ", 
                trials, ", Prob of success = ", prob, "\", type=\"h\")", 
                sep = ""))
            doItAndPrint(paste("points(.x, dnbinom(.x, size=", 
                trials, ", prob=", prob, "), pch=16)", sep = ""))
        }
        else if (fun == "Cumulative Probability") {
            justDoIt(paste("plot( stepfun(xmin:xmax, pnbinom((xmin-1):xmax", 
                ", size=", trials, ", prob=", prob, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of Failures until Target Successes\", ylab=\"Cumulative Probability\", main=\"Negative Binomial Distribution: Target = ", 
                trials, ", Probability of success = ", prob, 
                "\")", sep = ""))
            logger(paste("plot( stepfun(xmin:xmax, pnbinom((xmin-1):xmax", 
                ", size=", trials, ", prob=", prob, ")), verticals=F, do.p=F,", 
                " xlab=\"Number of Failures until Target Successes\", ylab=\"Cumulative Probability\", main=\"Negative Binomial Distribution: Target = ", 
                trials, ", Probability of success = ", prob, 
                "\")", sep = ""))
            doItAndPrint("abline( h = 1, lty = 2, col = \"grey\" )")
            doItAndPrint("abline( h = 0, lty = 2, col = \"grey\" )")
            doItAndPrint(paste("points( .x, pnbinom(.x,", " size=", 
                trials, ", prob=", prob, "), pch = 16, cex=1.2 )", 
                sep = ""))
            doItAndPrint(paste("points( .x, pnbinom((xmin-1):(xmax-1),", 
                " size=", trials, ", prob=", prob, "), pch = 1, cex=1.2 )", 
                sep = ""))
        }
        else if (fun == "Quantile Function") {
            doItAndPrint(paste("q <- pnbinom((xmin-1):xmax, size =", 
                trials, ",prob = ", prob, ")", sep = ""))
            main <- paste("\"Negative Binomial QF\n size=", trials, 
                ", prob=", prob, "\"", sep = "")
            ylab <- "\"Number of Failures until Target Successes\""
            xlab <- "\"Cumulative Probability\""
            command <- paste("plot(q[1:(length(q)-1)], .x, xlim=c(0,1), main =", 
                main, ", xlab = ", xlab, ",ylab = ", ylab, ", cex = 1.2)", 
                sep = "")
            justDoIt(command)
            logger(command)
            doItAndPrint(paste("abline( v = 0, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("abline( v = 1, lty = 2, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("points(q[2:(length(q))], .x, pch=16, cex=1.2)", 
                sep = ""))
            doItAndPrint(paste("segments(  q[1:(length(q)-1)], .x, q[2:(length(q))], .x)",
                sep = ""))
        }
        doItAndPrint("abline( h = 0, lty = 2, col = \"grey\" )")
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dnbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("size (target number of successes)")), 
        trialsEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot probability mass function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(trialsEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = trialsEntry)
}


`negbinomialMass.ipsur` <-
function () 
{
    checkTrials <- function(trials) {
        RcmdrTkmessageBox(message = sprintf(gettextRcmdr("Number of trials, %d, is large.\nCreate long output?"), 
            trials), icon = "warning", type = "yesno", default = "no")
    }
    initializeDialog(title = gettextRcmdr("Negative Binomial Probabilities"))
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    onOK <- function() {
        closeDialog()
        trials <- tclvalue(trialsVar)
        if (trials == "") {
            errorCondition(recall = negbinomialMass.ipsur, message = gettextRcmdr("Number of trials not specified."))
            return()
        }
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = negbinomialMass.ipsur, message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        command <- paste("xmin <- qnbinom(.00005, size=", trials, ", prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qnbinom(.99995, size=", trials, ", prob=", 
            prob, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".Table <- data.frame(Pr=dnbinom(xmin:xmax, size=", 
            trials, ", prob=", prob, "))", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        logger(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        justDoIt(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        doItAndPrint(".Table")
        logger("remove(.Table, xmin, xmax)")
        remove(.Table, xmin, xmax, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dnbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("size (target number of successes)")), 
        trialsEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(trialsEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    dialogSuffix(rows = 3, columns = 2, focus = trialsEntry)
}



`negbinomialProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Negative Binomial Probabilities "))
    ProbabilitiesVar <- tclVar("")
    ProbabilitiesEntry <- tkentry(top, width = "30", textvariable = ProbabilitiesVar)
    sizeVar <- tclVar("1")
    sizeEntry <- tkentry(top, width = "6", textvariable = sizeVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(ProbabilitiesVar))
        if ("" == quantiles) {
            errorCondition(recall = negbinomialProbabilities.ipsur, 
                message = gettextRcmdr("No probabilities specified."))
            return()
        }
        size <- tclvalue(sizeVar)
        prob <- tclvalue(probVar)
        if (size == "") {
            errorCondition(recall = negbinomialProbabilities.ipsur, 
                message = gettextRcmdr("The size parameter was not specified."))
            return()
        }
        if (prob == "") {
            errorCondition(recall = negbinomialProbabilities.ipsur, 
                message = gettextRcmdr("The probability parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("pnbinom(c(", quantiles, "), size=", 
            size, ", prob=", prob, ",  lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pnbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        ProbabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("size (target number of successes)")), 
        sizeEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(ProbabilitiesEntry, sticky = "w")
    tkgrid.configure(sizeEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = ProbabilitiesEntry)
}


`negbinomialQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Negative Binomial Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    sizeVar <- tclVar("1")
    sizeEntry <- tkentry(top, width = "6", textvariable = sizeVar)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = negbinomialQuantiles.ipsur, 
                message = gettextRcmdr("No probabilities specified."))
            return()
        }
        size <- tclvalue(sizeVar)
        prob <- tclvalue(probVar)
        if (size == "") {
            errorCondition(recall = negbinomialQuantiles.ipsur, 
                message = gettextRcmdr("The size parameter was not specified."))
            return()
        }
        if (prob == "") {
            errorCondition(recall = negbinomialQuantiles.ipsur, 
                message = gettextRcmdr("The probability parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qnbinom(c(", quantiles, "), size=", 
            size, ", prob=", prob, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qnbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("size (target number of successes)")), 
        sizeEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
        probEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(sizeEntry, sticky = "w")
    tkgrid.configure(probEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`PoissonDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Poisson Distribution"))
    meanVar <- tclVar("1")
    meanEntry <- tkentry(top, width = "6", textvariable = meanVar)
    functionVar <- tclVar("Probability")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Probability")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        mean <- tclvalue(meanVar)
        if (mean == "") {
            errorCondition(recall = PoissonDistributionPlot.ipsur, 
                message = gettextRcmdr("Mean not specified."))
            return()
        }
        fun <- tclvalue(functionVar)
        command <- paste("xmin <- qpois(.0005, lambda=", mean, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qpois(.9995, lambda=", mean, ")", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".x <- xmin:xmax", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        if (fun == "Probability") {
            doItAndPrint(paste("plot(.x, dpois(.x, lambda=", 
                mean, "), xlab=\"x\", ylab=\"Probability Mass\", main=\"Poisson Distribution: lambda = ", 
                mean, "\", type=\"h\")", sep = ""))
            doItAndPrint(paste("points(.x, dpois(.x, lambda=", 
                mean, "), pch=16)", sep = ""))
        }
        else if (fun == "Cumulative Probability") {
            justDoIt(paste("plot( stepfun(xmin:xmax, ppois((xmin-1):xmax", 
                ", lambda=", mean, ")), verticals=F, do.p=F,", 
                " xlab=\"x\", ylab=\"Cumulative Probability\", main=\"Poisson Distribution: lambda = ", 
                mean, "\")", sep = ""))
            logger(paste("plot( stepfun(.x, ppois((xmin-1):xmax", 
                ", lambda=", mean, ")), verticals=F, do.p=F,", 
                " xlab=\"x\", ylab=\"Cumulative Probability\", main=\"Poisson Distribution: lambda = ", 
                mean, "\")", sep = ""))
            doItAndPrint(paste("points( .x, ppois(xmin:xmax,", 
                " lambda=", mean, "), pch = 16, cex=1.2 )", sep = ""))
            doItAndPrint(paste("points( .x, ppois((xmin-1):(xmax-1),", 
                " lambda=", mean, "), pch = 1, cex=1.2 )", sep = ""))
            doItAndPrint("abline( h = 1, lty = 2, col = \"grey\" )")
        }
        else if (fun == "Quantile Function") {
            doItAndPrint(paste("q <- ppois((xmin-1):xmax, lambda =", 
                mean, ")", sep = ""))
            main <- paste("\"Poisson QF: lambda=", mean, "\"", 
                sep = "")
            ylab <- "\"Poisson Quantiles\""
            xlab <- "\"Cumulative Probability\""
            command <- paste("plot(q[1:(length(q)-1)], .x, xlim=c(0,1), main =", 
                main, ", xlab = ", xlab, ",ylab = ", ylab, ", cex = 1.2)", 
                sep = "")
            justDoIt(command)
            logger(command)
            doItAndPrint(paste("abline( v = 0, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("abline( v = 1, lty = 2, col = \"grey\" )", 
                sep = ""))
            doItAndPrint(paste("points(q[2:(length(q))], .x, pch=16, cex=1.2)", 
                sep = ""))
            doItAndPrint(paste("segments(  q[1:(length(q)-1)], .x, q[2:(length(q))], .x)",
                sep = ""))
        }
        doItAndPrint("abline( h = 0, lty = 2, col = \"grey\" )")
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dpois")
    tkgrid(tklabel(top, text = gettextRcmdr("lambda (mean)")), 
        meanEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot probability mass function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(meanEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 4, columns = 2, focus = meanEntry)
}


`PoissonMass.ipsur` <-
function () 
{
    checkRange <- function(range) {
        RcmdrTkmessageBox(message = sprintf(gettextRcmdr("Range of values over which to plot, %d, is large.\nCreate long output?"), 
            range), icon = "warning", type = "yesno", default = "no")
    }
    initializeDialog(title = gettextRcmdr("Poisson Probabilities"))
    meanVar <- tclVar("1")
    meanEntry <- tkentry(top, width = "6", textvariable = meanVar)
    onOK <- function() {
        closeDialog()
        mean <- tclvalue(meanVar)
        if (mean == "") {
            errorCondition(recall = PoissonMass.ipsur, message = gettextRcmdr("The mean parameter was not specified."))
            return()
        }
        command <- paste("xmin <- qpois(.00005, lambda=", mean, ")", 
            sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste("xmax <- qpois(.99995, lambda=", mean, ")", 
            sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        command <- paste(".Table <- data.frame(Pr=dpois(xmin:xmax, lambda=", 
            mean, "))", sep = "")
        logger(paste(command, sep = ""))
        justDoIt(command)
        logger(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        justDoIt(paste("rownames(.Table) <- xmin:xmax", sep = ""))
        doItAndPrint(".Table")
        logger("remove(.Table, xmin, xmax)")
        remove(.Table, xmin, xmax, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dpois")
    tkgrid(tklabel(top, text = gettextRcmdr("lambda (mean)")), 
        meanEntry, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(meanEntry, sticky = "w")
    dialogSuffix(rows = 2, columns = 2, focus = meanEntry)
}


`poissonProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Poisson Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    lambdaVar <- tclVar("1")
    lambdaEntry <- tkentry(top, width = "6", textvariable = lambdaVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = poissonProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        lambda <- tclvalue(lambdaVar)
        tail <- tclvalue(tailVar)
        if (lambda == "") {
            errorCondition(recall = poissonProbabilities.ipsur, 
                message = gettextRcmdr("The mean parameter was not specified."))
            return()
        }
        doItAndPrint(paste("ppois(c(", probabilities, "), lambda=", 
            lambda, ", lower.tail=", tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "ppois")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("lambda (mean)")), 
        lambdaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(lambdaEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`poissonQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Poisson Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    lambdaVar <- tclVar("1")
    lambdaEntry <- tkentry(top, width = "6", textvariable = lambdaVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = poissonQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        lambda <- tclvalue(lambdaVar)
        tail <- tclvalue(tailVar)
        if (lambda == "") {
            errorCondition(recall = poissonQuantiles.ipsur, message = gettextRcmdr("The mean parameter was not specified."))
            return()
        }
        doItAndPrint(paste("qpois(c(", quantiles, "), lambda=", 
            lambda, ", lower.tail=", tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qpois")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("lambda (mean)")), 
        lambdaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(lambdaEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}
