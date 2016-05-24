
# Last modified Feb 14, 2008

`betaDistributionPlot.ipsur` <-
function()
{
    initializeDialog(title = gettextRcmdr("Beta Distribution"))
    shape1Var <- tclVar("1")
    shape1Entry <- tkentry(top, width = "6", textvariable = shape1Var)
    shape2Var <- tclVar("1")
    shape2Entry <- tkentry(top, width = "6", textvariable = shape2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        shape1 <- tclvalue(shape1Var)
        shape2 <- tclvalue(shape2Var)
        ncp <- tclvalue(ncpVar)
        fun <- tclvalue(functionVar)
        if (shape1 == "") {
            errorCondition(recall = betaDistributionPlot.ipsur, 
                message = gettextRcmdr("Shape1 parameter not specified."))
            return()
        }
        if (shape2 == "") {
            errorCondition(recall = betaDistributionPlot.ipsur, 
                message = gettextRcmdr("Shape2 parameter not specified."))
            return()
        }
        if (ncp == "") {
            errorCondition(recall = betaDistributionPlot.ipsur, 
                message = gettextRcmdr("Noncentrality parameter not specified."))
            return()
        }
        fn <- if (fun == "Density") 
            "dbeta"
        else "pbeta"
        if (fun == "Density") {
            command <- paste("xmin <- qbeta(.00005, shape1=", shape1, 
                ", shape2=", shape2, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qbeta(.99995, shape1=", shape1, 
                ", shape2=", shape2, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dbeta(.x, shape1=", 
                shape1, ", shape2=", shape2, ", ncp=", ncp, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Beta Distribution: \", alpha, \" = ", 
                shape1, ", \", beta, \" = ", shape2, ", \", delta, \" = ", 
                ncp, "\")), type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qbeta(.00005, shape1=", shape1, 
                ", shape2=", shape2, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qbeta(.99995, shape1=", shape1, 
                ", shape2=", shape2, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pbeta(.x, shape1=", 
                shape1, ", shape2=", shape2, ", ncp=", ncp, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Beta Distribution: \", alpha, \" = ", 
                shape1, ", \", beta, \" = ", shape2, ", \", delta, \" = ", 
                ncp, "\")), type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qbeta(.00005, shape1=", shape1, 
                ", shape2=", shape2, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qbeta(.99995, shape1=", shape1, 
                ", shape2=", shape2, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qbeta(.x, shape1=", 
                shape1, ", shape2=", shape2, ", ncp=", ncp, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Beta QF: \", alpha, \" = ", 
                shape1, ", \", beta, \" = ", shape2, ", \", delta, \" = ", 
                ncp, "\")), type=\"l\")", sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dbeta")
    tkgrid(tklabel(top, text = gettextRcmdr("shape1 ")), shape1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape2 ")), shape2Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (non-centrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(shape1Entry, sticky = "w")
    tkgrid.configure(shape2Entry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = shape1Entry)
}

`betaProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Beta Probabilities "))
    ProbabilitiesVar <- tclVar("")
    ProbabilitiesEntry <- tkentry(top, width = "30", textvariable = ProbabilitiesVar)
    shape1Var <- tclVar("1")
    shape1Entry <- tkentry(top, width = "6", textvariable = shape1Var)
    shape2Var <- tclVar("1")
    shape2Entry <- tkentry(top, width = "6", textvariable = shape2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(ProbabilitiesVar))
        if ("" == quantiles) {
            errorCondition(recall = betaProbabilities.ipsur, 
                message = gettextRcmdr("No probabilities specified."))
            return()
        }
        shape1 <- tclvalue(shape1Var)
        shape2 <- tclvalue(shape2Var)
        ncp <- tclvalue(ncpVar)
        if (shape1 == "") {
            errorCondition(recall = betaProbabilities.ipsur, 
                message = gettextRcmdr("Shape1 parameter not specified."))
            return()
        }
        if (shape2 == "") {
            errorCondition(recall = betaProbabilities.ipsur, 
                message = gettextRcmdr("Shape2 parameter not specified"))
            return()
        }
        if (ncp == "") {
            errorCondition(recall = betaProbabilities.ipsur, 
                message = gettextRcmdr("Noncentrality parameter not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("pbeta(c(", quantiles, "), shape1=", 
            shape1, ", shape2=", shape2, ", ncp=", ncp, ", lower.tail=", 
            tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pbeta")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        ProbabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape1 ")), shape1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape2 ")), shape2Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (non-centrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(ProbabilitiesEntry, sticky = "w")
    tkgrid.configure(shape1Entry, sticky = "w")
    tkgrid.configure(shape2Entry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = ProbabilitiesEntry)
}

`betaQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Beta Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    shape1Var <- tclVar("1")
    shape1Entry <- tkentry(top, width = "6", textvariable = shape1Var)
    shape2Var <- tclVar("1")
    shape2Entry <- tkentry(top, width = "6", textvariable = shape2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = betaQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        shape1 <- tclvalue(shape1Var)
        shape2 <- tclvalue(shape2Var)
        ncp <- tclvalue(ncpVar)
        if (shape1 == "") {
            errorCondition(recall = betaQuantiles.ipsur, message = gettextRcmdr("Shape1 parameter not specified."))
            return()
        }
        if (shape2 == "") {
            errorCondition(recall = betaQuantiles.ipsur, message = gettextRcmdr("Shape2 parameter not specified."))
            return()
        }
        if (ncp == "") {
            errorCondition(recall = betaQuantiles.ipsur, message = gettextRcmdr("Noncentrality parameter not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qbeta(c(", quantiles, "), shape1=", 
            shape1, ", shape2=", shape2, ", ncp=", ncp, ", lower.tail=", 
            tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qbeta")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape1 ")), shape1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape2 ")), shape2Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (non-centrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(shape1Entry, sticky = "w")
    tkgrid.configure(shape2Entry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}

`birthdayProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Birthday Problem"))
    classesVar <- tclVar("365")
    classesEntry <- tkentry(top, width = "12", textvariable = classesVar)
    coincidentVar <- tclVar("2")
    coincidentEntry <- tkentry(top, width = "12", textvariable = coincidentVar)
    optionVar <- tclVar("probability")
    distFrame <- tkframe(top)
    probButton <- tkradiobutton(distFrame, variable = optionVar, 
        value = "probability")
    quantileButton <- tkradiobutton(distFrame, variable = optionVar, 
        value = "quantile")
    probFrame <- tkframe(distFrame)
    nVar <- tclVar("")
    probField <- tkentry(probFrame, width = "6", textvariable = nVar)
    quantileFrame <- tkframe(distFrame)
    probVar <- tclVar("")
    quantileField <- tkentry(quantileFrame, width = "6", textvariable = probVar)
    onOK <- function() {
        closeDialog()
        classes <- as.numeric(tclvalue(classesVar))
        coincident <- as.numeric(tclvalue(coincidentVar))
        n <- as.numeric(tclvalue(nVar))
        prob <- as.numeric(tclvalue(probVar))
        option <- tclvalue(optionVar)
        if (is.na(classes)) {
            errorCondition(recall = birthdayProbabilities.ipsur, 
                message = gettextRcmdr("Class parameter not specified."))
            return()
        }
        if (is.na(coincident)) {
            errorCondition(recall = birthdayProbabilities.ipsur, 
                message = gettextRcmdr("Coincident parameter not specified."))
            return()
        }
        if (option == "probability") {
            if (is.na(n)) {
                errorCondition(recall = birthdayProbabilities.ipsur, 
                  message = gettextRcmdr("Number of people not specified."))
                return()
            }
            doItAndPrint(paste("pbirthday.ipsur(n=", n, ", classes=", 
                classes, ", coincident=", coincident, ")", sep = ""))
            tkfocus(CommanderWindow())
        }
        if (option == "quantile") {
            if (is.na(prob)) {
                errorCondition(recall = birthdayProbabilities.ipsur, 
                  message = gettextRcmdr("Probability not specified."))
                return()
            }
            doItAndPrint(paste("qbirthday.ipsur(prob=", prob, 
                ", classes=", classes, ", coincident=", coincident, 
                ")", sep = ""))
            tkfocus(CommanderWindow())
        }
    }
    OKCancelHelp(helpSubject = "qbirthday")
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters"), fg = "red"), 
        columnspan = 6, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("classes (distinct categories)")), 
        classesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("coincident (num people falling in same category) ")), 
        coincidentEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr(""), fg = "red"), 
        columnspan = 6, sticky = "w")
    tkgrid(tklabel(distFrame, text = gettextRcmdr("Select Option"), 
        fg = "red"), columnspan = 6, sticky = "w")
    tkgrid.configure(classesEntry, sticky = "w")
    tkgrid.configure(coincidentEntry, sticky = "w")
    tkgrid(tklabel(probFrame, text = gettextRcmdr("n (num of people) = ")), 
        probField, sticky = "w")
    tkgrid(tklabel(distFrame, text = gettextRcmdr("Probabilities"), 
        fg = "blue"), probButton, probFrame, sticky = "w")
    tkgrid(tklabel(quantileFrame, text = gettextRcmdr("prob (of coincidence) = ")), 
        quantileField, sticky = "w")
    tkgrid(tklabel(distFrame, text = gettextRcmdr("Quantiles"), 
        fg = "blue"), quantileButton, quantileFrame, sticky = "w")
    tkgrid(distFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = classesEntry)
}

`cauchyDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Cauchy Distribution 03"))
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        location <- tclvalue(locationVar)
        scale1 <- tclvalue(scale1Var)
        fun <- tclvalue(functionVar)
        if (location == "") {
            errorCondition(recall = cauchyDistributionPlot.ipsur, 
                message = gettextRcmdr("The location parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = cauchyDistributionPlot.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qcauchy(.00005, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qcauchy(.99995, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dcauchy(.x, location=", 
                location, ", scale=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Cauchy Distribution: \", location, \" = ", 
                location, ", \", scale, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qcauchy(.00005, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qcauchy(.99995, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pcauchy(.x, location=", 
                location, ", scale=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Cauchy Distribution: \", location, \" = ", 
                location, ", \", scale, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qcauchy(.00005, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qcauchy(.99995, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qcauchy(.x, location=", 
                location, ", scale=", scale1, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", \n\t\t\t\tmain=expression(paste(\"Cauchy Distribution: \", location, \" = ", 
                location, ", \", scale, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dcauchy")
    tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(locationEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = locationEntry)
}


`cauchyProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Cauchy Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = cauchyProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        location <- tclvalue(locationVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (location == "") {
            errorCondition(recall = cauchyProbabilities.ipsur, 
                message = gettextRcmdr("The location parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = cauchyProbabilities.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        doItAndPrint(paste("pcauchy(c(", probabilities, "), location=", 
            location, ", scale=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pcauchy")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(locationEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`cauchyQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Cauchy Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = cauchyQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        location <- tclvalue(locationVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (location == "") {
            errorCondition(recall = cauchyQuantiles.ipsur, message = gettextRcmdr("The location parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = cauchyQuantiles.ipsur, message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        doItAndPrint(paste("qcauchy(c(", quantiles, "), location=", 
            location, ", scale=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qcauchy")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(locationEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`chisqProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Chi-Squared Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = chisqProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        df <- tclvalue(dfVar)
        if (df == "") {
            errorCondition(recall = chisqProbabilities.ipsur, 
                message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = chisqProbabilities.ipsur, 
                message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("pchisq(c(", probabilities, "), df=", 
            df, ", ncp=", ncp, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pchisq")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(OKbutton, cancelButton, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(dfEntry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = probabilitiesEntry)
}


`chisqQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Chi-Squared Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = chisqQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        df <- tclvalue(dfVar)
        if (df == "") {
            errorCondition(recall = chisqQuantiles.ipsur, message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = chisqQuantiles.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qchisq(c(", quantiles, "), df=", 
            df, ", ncp=", ncp, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qchisq")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(dfEntry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = quantilesEntry)
}


`chisquareDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Chi-squared Distribution"))
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        df <- tclvalue(dfVar)
        if (df == "") {
            errorCondition(recall = chisquareDistributionPlot.ipsur, 
                message = gettextRcmdr("The degrees of freedom were not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = chisquareDistributionPlot.ipsur, 
                message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        fun <- tclvalue(functionVar)
        if (fun == "Density") {
            command <- paste("xmin <- qchisq(.00005, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qchisq(.99995, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dchisq(.x, df=", df, 
                ", ncp=", ncp, "), xlab=expression(chi^2), ylab=\"", 
                fun, "\", main=\"Chi-Squared Distribution: df = ", 
                df, ", ncp=", ncp, "\", type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qchisq(.00005, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qchisq(.99995, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pchisq(.x, df=", df, 
                ", ncp=", ncp, "), xlab=expression(chi^2), ylab=\"", 
                fun, "\", main=\"Chi-Squared Distribution: df = ", 
                df, ", ncp=", ncp, "\", type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qchisq(.00005, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qchisq(.99995, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qchisq(.x, df=", df, 
                ", ncp=", ncp, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=\"Chi-Squared Distribution: df = ", 
                df, ", ncp=", ncp, "\", type=\"l\")", sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dchisq")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(dfEntry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 4, columns = 2, focus = dfEntry)
}


`expDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Exponential Distribution"))
    rateVar <- tclVar("1")
    rateEntry <- tkentry(top, width = "6", textvariable = rateVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        rate <- tclvalue(rateVar)
        fun <- tclvalue(functionVar)
        if (rate == "") {
            errorCondition(recall = expDistributionPlot.ipsur, 
                message = gettextRcmdr("The rate parameter was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qexp(.00005, rate=", rate, ")", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qexp(.99995, rate=", rate, ")", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dexp(.x, rate=", rate, 
                " ), xlab=\"x\", ylab=\"", fun, "\", main=expression(paste(\"Exponential Distribution: \", lambda, \" =", 
                rate, "\")), type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qexp(.00005, rate=", rate, ")", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qexp(.99995, rate=", rate, ")", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pexp(.x, rate=", rate, 
                " ), xlab=\"x\", ylab=\"", fun, "\", main=expression(paste(\"Exponential Distribution: \", lambda, \" =", 
                rate, "\")), type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qexp(.00005, rate=", rate, ")", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qexp(.99995, rate=", rate, ")", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qexp(.x, rate=", rate, 
                " ), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Exponential Distribution: \", lambda, \" =", 
                rate, "\")), type=\"l\")", sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dexp")
    tkgrid(tklabel(top, text = gettextRcmdr("rate (of arrivals in unit time)")), 
        rateEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(rateEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = rateEntry)
}


`expProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Exponential Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    rateVar <- tclVar("1")
    rateEntry <- tkentry(top, width = "6", textvariable = rateVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = expProbabilities.ipsur, message = gettextRcmdr("No values specified."))
            return()
        }
        rate <- tclvalue(rateVar)
        tail <- tclvalue(tailVar)
        if (rate == "") {
            errorCondition(recall = expProbabilities.ipsur, message = gettextRcmdr("The rate parameter was not specified."))
            return()
        }
        doItAndPrint(paste("pexp(c(", probabilities, "), rate=", 
            rate, ", lower.tail=", tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pexp")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("rate (of arrivals in unit time)")), 
        rateEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(rateEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`expQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Exponential Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    rateVar <- tclVar("1")
    rateEntry <- tkentry(top, width = "6", textvariable = rateVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = expQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        rate <- tclvalue(rateVar)
        tail <- tclvalue(tailVar)
        if (rate == "") {
            errorCondition(recall = expQuantiles.ipsur, message = gettextRcmdr("The rate parameter was not specified."))
            return()
        }
        doItAndPrint(paste("qexp(c(", quantiles, "), rate=", 
            rate, ", lower.tail=", tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qexp")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("rate (of arrivals in unit time)")), 
        rateEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(rateEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`FDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("F Distribution"))
    df1Var <- tclVar("1")
    df2Var <- tclVar("1")
    df1Entry <- tkentry(top, width = "6", textvariable = df1Var)
    df2Entry <- tkentry(top, width = "6", textvariable = df2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        df1 <- tclvalue(df1Var)
        df2 <- tclvalue(df2Var)
        if (df1 == "" || df2 == "") {
            errorCondition(recall = FDistributionPlot.ipsur, 
                message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = FDistributionPlot.ipsur, 
                message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        fun <- tclvalue(functionVar)
        if (fun == "Density") {
            command <- paste("xmin <- qf(.00005, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qf(.99995, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, df(.x, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, "), xlab=\"f\", ylab=\"", 
                fun, "\", main=\"F Distribution: Numerator df = ", 
                df1, ", Denominator df = ", df2, "ncp = ", ncp, 
                "\", type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qf(.00005, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qf(.99995, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pf(.x, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, "), xlab=\"f\", ylab=\"", 
                fun, "\", main=\"F Distribution: Numerator df = ", 
                df1, ", Denominator df = ", df2, "ncp = ", ncp, 
                "\", type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qf(.00005, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qf(.99995, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qf(.x, df1=", df1, ", df2=", 
                df2, ", ncp=", ncp, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=\"F Distribution: Numerator df = ", 
                df1, ", Denominator df = ", df2, "ncp = ", ncp, 
                "\", type=\"l\")", sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "df")
    tkgrid(tklabel(top, text = gettextRcmdr("df1 (numerator degrees of freedom)")), 
        df1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df2 (denominator degrees of freedom)")), 
        df2Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Pot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(df1Entry, sticky = "w")
    tkgrid.configure(df2Entry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = df1Entry)
}


`FProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("F Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    df1Var <- tclVar("1")
    df1Entry <- tkentry(top, width = "6", textvariable = df1Var)
    df2Var <- tclVar("1")
    df2Entry <- tkentry(top, width = "6", textvariable = df2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = FProbabilities.ipsur, message = gettextRcmdr("Values not specified."))
            return()
        }
        df1 <- tclvalue(df1Var)
        df2 <- tclvalue(df2Var)
        if (df1 == "" || df2 == "") {
            errorCondition(recall = FProbabilities.ipsur, message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = FProbabilities.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("pf(c(", probabilities, "), df1=", 
            df1, ", df2=", df2, ", ncp=", ncp, ", lower.tail=", 
            tail == "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pf")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df1 (numerator degrees of freedom)")), 
        df1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df2 (denominator degrees of freedom)")), 
        df2Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(df1Entry, sticky = "w")
    tkgrid.configure(df2Entry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = probabilitiesEntry)
}


`FQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("F Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    df1Var <- tclVar("1")
    df1Entry <- tkentry(top, width = "6", textvariable = df1Var)
    df2Var <- tclVar("1")
    df2Entry <- tkentry(top, width = "6", textvariable = df2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = FQuantiles.ipsur, message = gettextRcmdr("Probabilities not specified"))
            return()
        }
        df1 <- tclvalue(df1Var)
        df2 <- tclvalue(df2Var)
        if (df1 == "" || df2 == "") {
            errorCondition(recall = FQuantiles.ipsur, message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = FQuantiles.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qf(c(", quantiles, "), df1=", df1, 
            ", df2=", df2, ", ncp=", ncp, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qf")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df1 (numerator degrees of freedom)")), 
        df1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df2 (denominator degrees of freedom)")), 
        df2Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(df1Entry, sticky = "w")
    tkgrid.configure(df2Entry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`gammaDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Gamma Distribution"))
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        shape <- tclvalue(shapeVar)
        scale1 <- tclvalue(scale1Var)
        fun <- tclvalue(functionVar)
        if (shape == "") {
            errorCondition(recall = gammaDistributionPlot.ipsur, 
                message = gettextRcmdr("The shape parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = gammaDistributionPlot.ipsur, 
                message = gettextRcmdr("The rate parameter was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qgamma(.00005, shape=", shape, 
                ", rate=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qgamma(.99995, shape=", shape, 
                ", rate=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dgamma(.x, shape=", 
                shape, ", rate=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Gamma Distribution: \", alpha, \" = ", 
                shape, ", \", lambda, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qgamma(.00005, shape=", shape, 
                ", rate=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qgamma(.99995, shape=", shape, 
                ", rate=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pgamma(.x, shape=", 
                shape, ", rate=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Gamma Distribution: \", alpha, \" = ", 
                shape, ", \", lambda, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qgamma(.00005, shape=", shape, 
                ", rate=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qgamma(.99995, shape=", shape, 
                ", rate=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qgamma(.x, shape=", 
                shape, ", rate=", scale1, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Gamma Distribution: \", alpha, \" = ", 
                shape, ", \", lambda, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dgamma")
    tkgrid(tklabel(top, text = gettextRcmdr("shape ")), shapeEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale (=1/scale)")), 
        scale1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(shapeEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = shapeEntry)
}


`gammaProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Gamma Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = gammaProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        shape <- tclvalue(shapeVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (shape == "") {
            errorCondition(recall = gammaProbabilities.ipsur, 
                message = gettextRcmdr("The shape parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = gammaProbabilities.ipsur, 
                message = gettextRcmdr("The rate parameter was not specified."))
            return()
        }
        doItAndPrint(paste("pgamma(c(", probabilities, "), shape=", 
            shape, ", rate=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pgamma")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape ")), shapeEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale (=1/scale)")), 
        scale1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(shapeEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`gammaQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Gamma Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = gammaQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        shape <- tclvalue(shapeVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (shape == "") {
            errorCondition(recall = gammaQuantiles.ipsur, message = gettextRcmdr("The shape parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = gammaQuantiles.ipsur, message = gettextRcmdr("The rate parameter was not specified."))
            return()
        }
        doItAndPrint(paste("qgamma(c(", quantiles, "), shape=", 
            shape, ", rate=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qgamma")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape ")), shapeEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("rate (=1/scale)")), 
        scale1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(shapeEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`logisticDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Logistic Distribution"))
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        location <- tclvalue(locationVar)
        scale1 <- tclvalue(scale1Var)
        fun <- tclvalue(functionVar)
        if (location == "") {
            errorCondition(recall = logisticDistributionPlot.ipsur, 
                message = gettextRcmdr("The location parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = logisticDistributionPlot.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qlogis(.00005, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qlogis(.99995, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dlogis(.x, location=", 
                location, ", scale=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Logistic Distribution: \", mu, \" = ", 
                location, ", \", beta, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qlogis(.00005, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qlogis(.99995, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, plogis(.x, location=", 
                location, ", scale=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Logistic Distribution: \", mu, \" = ", 
                location, ", \", beta, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qlogis(.00005, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qlogis(.99995, location=", location, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qlogis(.x, location=", 
                location, ", scale=", scale1, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Logistic Distribution: \", mu, \" = ", 
                location, ", \", beta, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dlogis")
    tkgrid(tklabel(top, text = gettextRcmdr("location ")), locationEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale ")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantil function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(locationEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = locationEntry)
}


`logisticProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Logistic Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = logisticProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        location <- tclvalue(locationVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (location == "") {
            errorCondition(recall = logisticProbabilities.ipsur, 
                message = gettextRcmdr("The location parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = logisticProbabilities.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        doItAndPrint(paste("plogis(c(", probabilities, "), location=", 
            location, ", scale=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "plogis")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("location ")), locationEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale ")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(locationEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`logisticQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Logistic Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = logisticQuantiles.ipsur, 
                message = gettextRcmdr("No probabilities specified."))
            return()
        }
        location <- tclvalue(locationVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (location == "") {
            errorCondition(recall = logisticQuantiles.ipsur, 
                message = gettextRcmdr("The location parameter must be a real number."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = logisticQuantiles.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        doItAndPrint(paste("qlogis(c(", quantiles, "), location=", 
            location, ", scale=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qlogis")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("location ")), locationEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale ")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(locationEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`lognormalDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Log Normal Distribution"))
    mulogVar <- tclVar("0")
    mulogEntry <- tkentry(top, width = "6", textvariable = mulogVar)
    sigmalogVar <- tclVar("1")
    sigmalogEntry <- tkentry(top, width = "6", textvariable = sigmalogVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        mulog <- tclvalue(mulogVar)
        sigmalog <- tclvalue(sigmalogVar)
        fun <- tclvalue(functionVar)
        if (mulog == "") {
            errorCondition(recall = lognormalDistributionPlot.ipsur, 
                message = gettextRcmdr("The mean was not specified."))
            return()
        }
        if (sigmalog == "") {
            errorCondition(recall = lognormalDistributionPlot.ipsur, 
                message = gettextRcmdr("The standard deviation was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qlnorm(.00005, meanlog=", mulog, 
                ", sdlog=", sigmalog, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qlnorm(.99995, meanlog=", mulog, 
                ", sdlog=", sigmalog, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dlnorm(.x, meanlog=", 
                mulog, ", sdlog=", sigmalog, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Log Normal Distribution: \", mulog, \" = ", 
                mulog, ", \", sigmalog, \" = ", sigmalog, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qlnorm(.00005, meanlog=", mulog, 
                ", sdlog=", sigmalog, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qlnorm(.99995, meanlog=", mulog, 
                ", sdlog=", sigmalog, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, plnorm(.x, meanlog=", 
                mulog, ", sdlog=", sigmalog, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Log Normal Distribution: \", mulog, \" = ", 
                mulog, ", \", sigmalog, \" = ", sigmalog, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qlnorm(.00005, meanlog=", mulog, 
                ", sdlog=", sigmalog, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qlnorm(.99995, meanlog=", mulog, 
                ", sdlog=", sigmalog, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qlnorm(.x, meanlog=", 
                mulog, ", sdlog=", sigmalog, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Log Normal Distribution: \", mulog, \" = ", 
                mulog, ", \", sigmalog, \" = ", sigmalog, "\")), type=\"l\")", 
                sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dlnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("meanlog (mean of dist'n on log scale)")), 
        mulogEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("sdlog (std dev of dist'n on log scale)")), 
        sigmalogEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(mulogEntry, sticky = "w")
    tkgrid.configure(sigmalogEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = mulogEntry)
}


`lognormalProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Log Normal Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    logmuVar <- tclVar("0")
    logmuEntry <- tkentry(top, width = "6", textvariable = logmuVar)
    logsigmaVar <- tclVar("1")
    logsigmaEntry <- tkentry(top, width = "6", textvariable = logsigmaVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = lognormalProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        logmu <- tclvalue(logmuVar)
        logsigma <- tclvalue(logsigmaVar)
        tail <- tclvalue(tailVar)
        if (logmu == "") {
            errorCondition(recall = lognormalProbabilities.ipsur, 
                message = gettextRcmdr("The mean was not specified."))
            return()
        }
        if (logsigma == "") {
            errorCondition(recall = lognormalProbabilities.ipsur, 
                message = gettextRcmdr("The standard deviation was not specified."))
            return()
        }
        doItAndPrint(paste("plnorm(c(", probabilities, "), meanlog=", 
            logmu, ", sdlog=", logsigma, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "plnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("meanlog (mean of dist'n on log scale )")), 
        logmuEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("sdlog (std dev of dist'n on log scale)")), 
        logsigmaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(logmuEntry, sticky = "w")
    tkgrid.configure(logsigmaEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`lognormalQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Log Normal Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    logmuVar <- tclVar("0")
    logmuEntry <- tkentry(top, width = "6", textvariable = logmuVar)
    logsigmaVar <- tclVar("1")
    logsigmaEntry <- tkentry(top, width = "6", textvariable = logsigmaVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = lognormalQuantiles.ipsur, 
                message = gettextRcmdr("No probabilities specified."))
            return()
        }
        logmu <- tclvalue(logmuVar)
        logsigma <- tclvalue(logsigmaVar)
        tail <- tclvalue(tailVar)
        if (logmu == "") {
            errorCondition(recall = lognormalQuantiles.ipsur, 
                message = gettextRcmdr("The mean was not specified."))
            return()
        }
        if (logsigma == "") {
            errorCondition(recall = lognormalQuantiles.ipsur, 
                message = gettextRcmdr("The standard deviation was not specified."))
            return()
        }
        doItAndPrint(paste("qlnorm(c(", quantiles, "), meanlog=", 
            logmu, ", sdlog=", logsigma, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qlnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("meanlog (mean of dist'n on log scale)")), 
        logmuEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("sdlog (std dev of dist'n on log scale)")), 
        logsigmaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(logmuEntry, sticky = "w")
    tkgrid.configure(logsigmaEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`normalDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Normal Distribution"))
    muVar <- tclVar("0")
    muEntry <- tkentry(top, width = "6", textvariable = muVar)
    sigmaVar <- tclVar("1")
    sigmaEntry <- tkentry(top, width = "6", textvariable = sigmaVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        mu <- tclvalue(muVar)
        sigma <- tclvalue(sigmaVar)
        fun <- tclvalue(functionVar)
        if (mu == "") {
            errorCondition(recall = normalDistributionPlot.ipsur, 
                message = gettextRcmdr("The mean was not specified."))
            return()
        }
        if (sigma == "") {
            errorCondition(recall = normalDistributionPlot.ipsur, 
                message = gettextRcmdr("The standard deviation was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qnorm(.00005, mean=", mu, ", sd=", 
                sigma, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qnorm(.99995, mean=", mu, ", sd=", 
                sigma, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dnorm(.x, mean=", mu, 
                ", sd=", sigma, "), xlab=\"x\", ylab=\"", fun, 
                "\", main=expression(paste(\"Normal Distribution: \", mu, \" = ", 
                mu, ", \", sigma, \" = ", sigma, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qnorm(.00005, mean=", mu, ", sd=", 
                sigma, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qnorm(.99995, mean=", mu, ", sd=", 
                sigma, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pnorm(.x, mean=", mu, 
                ", sd=", sigma, "), xlab=\"x\", ylab=\"", fun, 
                "\", main=expression(paste(\"Normal Distribution: \", mu, \" = ", 
                mu, ", \", sigma, \" = ", sigma, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qnorm(.00005, mean=", mu, ", sd=", 
                sigma, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qnorm(.99995, mean=", mu, ", sd=", 
                sigma, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qnorm(.x, mean=", mu, 
                ", sd=", sigma, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Normal Distribution: \", mu, \" = ", 
                mu, ", \", sigma, \" = ", sigma, "\")), type=\"l\")", 
                sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("mu (mean)")), muEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("sigma (standard deviation)")), 
        sigmaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(muEntry, sticky = "w")
    tkgrid.configure(sigmaEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = muEntry)
}


`normalProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Normal Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    muVar <- tclVar("0")
    muEntry <- tkentry(top, width = "6", textvariable = muVar)
    sigmaVar <- tclVar("1")
    sigmaEntry <- tkentry(top, width = "6", textvariable = sigmaVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        mu <- tclvalue(muVar)
        sigma <- tclvalue(sigmaVar)
        tail <- tclvalue(tailVar)
        if ("" == probabilities) {
            errorCondition(recall = normalProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        if (mu == "") {
            errorCondition(recall = normalProbabilities.ipsur, 
                message = gettextRcmdr("The mean was not specified."))
            return()
        }
        if (sigma == "") {
            errorCondition(recall = normalProbabilities.ipsur, 
                message = gettextRcmdr("The standard deviation was not specified."))
            return()
        }
        doItAndPrint(paste("pnorm(c(", probabilities, "), mean=", 
            mu, ", sd=", sigma, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("mu (mean)")), muEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("sigma (standard deviation)")), 
        sigmaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(muEntry, sticky = "w")
    tkgrid.configure(sigmaEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`normalQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Normal Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    muVar <- tclVar("0")
    muEntry <- tkentry(top, width = "6", textvariable = muVar)
    sigmaVar <- tclVar("1")
    sigmaEntry <- tkentry(top, width = "6", textvariable = sigmaVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        mu <- tclvalue(muVar)
        sigma <- tclvalue(sigmaVar)
        tail <- tclvalue(tailVar)
        if ("" == quantiles) {
            errorCondition(recall = normalQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        if (mu == "") {
            errorCondition(recall = normalQuantiles.ipsur, message = gettextRcmdr("The mean was not specified."))
            return()
        }
        if (sigma == "") {
            errorCondition(recall = normalQuantiles.ipsur, message = gettextRcmdr("The standard deviation was not specified."))
            return()
        }
        doItAndPrint(paste("qnorm(c(", quantiles, "), mean=", 
            mu, ", sd=", sigma, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("mu (mean)")), muEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("sigma (standard deviation)")), 
        sigmaEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(muEntry, sticky = "w")
    tkgrid.configure(sigmaEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`tDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("t Distribution"))
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        df <- tclvalue(dfVar)
        fun <- tclvalue(functionVar)
        if (df == "") {
            errorCondition(recall = tDistributionPlot.ipsur, 
                message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = tDistributionPlot.ipsur, 
                message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qt(.00005, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qt(.99995, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dt(.x, df=", df, ", ncp=", 
                ncp, "), xlab=\"t\", ylab=\"", fun, "\", main=\"t Distribution: df = ", 
                df, ", ncp=", ncp, "\", type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qt(.00005, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qt(.99995, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pt(.x, df=", df, ", ncp=", 
                ncp, "), xlab=\"t\", ylab=\"", fun, "\", main=\"t Distribution: df = ", 
                df, ", ncp=", ncp, "\", type=\"l\")", sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qt(.00005, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qt(.99995, df=", df, ", ncp=", 
                ncp, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qt(.x, df=", df, ", ncp=", 
                ncp, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=\"t Distribution: df = ", 
                df, ", ncp=", ncp, "\", type=\"l\")", sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dt")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(dfEntry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 4, columns = 2, focus = dfEntry)
}


`tProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("t Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        df <- as.numeric(tclvalue(dfVar))
        if ("" == probabilities) {
            errorCondition(recall = tProbabilities.ipsur, message = gettextRcmdr("No values specified."))
            return()
        }
        df <- tclvalue(dfVar)
        if (df == "") {
            errorCondition(recall = tProbabilities.ipsur, message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = tProbabilities.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("pt(c(", probabilities, "), df=", 
            df, ", ncp=", ncp, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pt")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(dfEntry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = probabilitiesEntry)
}


`tQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("t Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = tQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        df <- tclvalue(dfVar)
        if (df == "") {
            errorCondition(recall = tQuantiles.ipsur, message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        ncp <- tclvalue(ncpVar)
        if (ncp == "") {
            errorCondition(recall = tQuantiles.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
            return()
        }
        tail <- tclvalue(tailVar)
        doItAndPrint(paste("qt(c(", quantiles, "), df=", df, 
            ", ncp=", ncp, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qt")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(dfEntry, sticky = "w")
    tkgrid.configure(ncpEntry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = quantilesEntry)
}


`unifDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Uniform Distribution"))
    min1Var <- tclVar("0")
    min1Entry <- tkentry(top, width = "6", textvariable = min1Var)
    max1Var <- tclVar("1")
    max1Entry <- tkentry(top, width = "6", textvariable = max1Var)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        min1 <- tclvalue(min1Var)
        max1 <- tclvalue(max1Var)
        fun <- tclvalue(functionVar)
        if (min1 == "") {
            errorCondition(recall = unifDistributionPlot.ipsur, 
                message = gettextRcmdr("The lower limit(min) was not specified."))
            return()
        }
        if (max1 == "") {
            errorCondition(recall = unifDistributionPlot.ipsur, 
                message = gettextRcmdr("The upper limit(max) was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qunif(.00005, min=", min1, ", max=", 
                max1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qunif(.99995, min=", min1, ", max=", 
                max1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dunif(.x, min=", min1, 
                ", max=", max1, "), xlab=\"x\", ylab=\"", fun, 
                "\", main=expression(paste(\"Uniform Distribution: \", min, \" = ", 
                min1, ", \", max, \" = ", max1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qunif(.00005, min=", min1, ", max=", 
                max1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qunif(.99995, min=", min1, ", max=", 
                max1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, punif(.x, min=", min1, 
                ", max=", max1, "), xlab=\"x\", ylab=\"", fun, 
                "\", main=expression(paste(\"Uniform Distribution: \", min, \" = ", 
                min1, ", \", max, \" = ", max1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qunif(.00005, min=", min1, ", max=", 
                max1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qunif(.99995, min=", min1, ", max=", 
                max1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qunif(.x, min=", min1, 
                ", max=", max1, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Uniform Distribution: \", min, \" = ", 
                min1, ", \", max, \" = ", max1, "\")), type=\"l\")", 
                sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dunif")
    tkgrid(tklabel(top, text = gettextRcmdr("min (lower limit of the distribution)")), 
        min1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("max (upper limit of the distribution)")), 
        max1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(min1Entry, sticky = "w")
    tkgrid.configure(max1Entry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = min1Entry)
}


`uniformProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Uniform Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    min1Var <- tclVar("0")
    min1Entry <- tkentry(top, width = "6", textvariable = min1Var)
    max1Var <- tclVar("1")
    max1Entry <- tkentry(top, width = "6", textvariable = max1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = uniformProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        min1 <- tclvalue(min1Var)
        max1 <- tclvalue(max1Var)
        tail <- tclvalue(tailVar)
        if (min1 == "") {
            errorCondition(recall = uniformProbabilities.ipsur, 
                message = gettextRcmdr("The lower limit(min) was not specified."))
            return()
        }
        if (max1 == "") {
            errorCondition(recall = uniformProbabilities.ipsur, 
                message = gettextRcmdr("The upper limit(max) was not specified."))
            return()
        }
        doItAndPrint(paste("punif(c(", probabilities, "), min=", 
            min1, ", max=", max1, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "punif")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("min (lower limit of the distribution)")), 
        min1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("max (upper limit of the distribution)")), 
        max1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(min1Entry, sticky = "w")
    tkgrid.configure(max1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`uniformQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Uniform Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    min1Var <- tclVar("0")
    min1Entry <- tkentry(top, width = "6", textvariable = min1Var)
    max1Var <- tclVar("1")
    max1Entry <- tkentry(top, width = "6", textvariable = max1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = uniformQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        min1 <- tclvalue(min1Var)
        max1 <- tclvalue(max1Var)
        tail <- tclvalue(tailVar)
        if (min1 == "") {
            errorCondition(recall = uniformQuantiles.ipsur, message = gettextRcmdr("The lower limit(min) was not specified."))
            return()
        }
        if (max1 == "") {
            errorCondition(recall = uniformQuantiles.ipsur, message = gettextRcmdr("The upper limit(max) was not specified."))
            return()
        }
        doItAndPrint(paste("qunif(c(", quantiles, "), min=", 
            min1, ", max=", max1, ", lower.tail=", tail == "lower", 
            ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qunif")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("min (lower limit of the distribution)")), 
        min1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("max (upper limit of the distribution)")), 
        max1Entry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(min1Entry, sticky = "w")
    tkgrid.configure(max1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}


`weibullDistributionPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Weibull Distribution"))
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    functionVar <- tclVar("Density")
    densityButton <- tkradiobutton(top, variable = functionVar, 
        value = "Density")
    distributionButton <- tkradiobutton(top, variable = functionVar, 
        value = "Cumulative Probability")
    quantileButton <- tkradiobutton(top, variable = functionVar, 
        value = "Quantile Function")
    onOK <- function() {
        closeDialog()
        shape <- tclvalue(shapeVar)
        scale1 <- tclvalue(scale1Var)
        fun <- tclvalue(functionVar)
        if (shape == "") {
            errorCondition(recall = weibullDistributionPlot.ipsur, 
                message = gettextRcmdr("The shape parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = weibullDistributionPlot.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        if (fun == "Density") {
            command <- paste("xmin <- qweibull(.00005, shape=", shape, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qweibull(.99995, shape=", shape, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, dweibull(.x, shape=", 
                shape, ", scale=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Weibull Distribution: \", gamma, \" = ", 
                shape, ", \", beta, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Cumulative Probability") {
            command <- paste("xmin <- qweibull(.00005, shape=", shape, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qweibull(.99995, shape=", shape, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(xmin, xmax, length=100)", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, pweibull(.x, shape=", 
                shape, ", scale=", scale1, "), xlab=\"x\", ylab=\"", 
                fun, "\", main=expression(paste(\"Weibull Distribution: \", gamma, \" = ", 
                shape, ", \", beta, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            doItAndPrint("abline(h=0, col=\"gray\")")
        }
        else if (fun == "Quantile Function") {
            command <- paste("xmin <- qweibull(.00005, shape=", shape, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste("xmax <- qweibull(.99995, shape=", shape, 
                ", scale=", scale1, ")", sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            command <- paste(".x <- seq(.00005, .99995, length=100)", 
                sep = "")
            logger(paste(command, sep = ""))
            justDoIt(command)
            doItAndPrint(paste("plot(.x, qweibull(.x, shape=", 
                shape, ", scale=", scale1, "), xlab=\"Cumulative Probability\", ylab=\"Quantile\", main=expression(paste(\"Weibull Distribution: \", gamma, \" = ", 
                shape, ", \", beta, \" = ", scale1, "\")), type=\"l\")", 
                sep = ""))
            justDoIt(paste("abline(h=0, col = \"grey\")"))
            justDoIt(paste("abline(v=0, col = \"grey\")"))
            justDoIt(paste("abline(v=1, lty = 2, col = \"grey\")"))
        }
        remove(.x, xmin, xmax, envir = .GlobalEnv)
        logger("remove(.x, xmin, xmax)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "dweibull")
    tkgrid(tklabel(top, text = gettextRcmdr("shape ")), shapeEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale ")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot density function")), 
        densityButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot distribution function")), 
        distributionButton, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Plot quantile function")), 
        quantileButton, sticky = "e")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    tkgrid.configure(shapeEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(densityButton, sticky = "w")
    tkgrid.configure(distributionButton, sticky = "w")
    tkgrid.configure(quantileButton, sticky = "w")
    dialogSuffix(rows = 5, columns = 2, focus = shapeEntry)
}


`weibullProbabilities.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Weibull Probabilities"))
    probabilitiesVar <- tclVar("")
    probabilitiesEntry <- tkentry(top, width = "30", textvariable = probabilitiesVar)
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        probabilities <- gsub(" ", ",", tclvalue(probabilitiesVar))
        if ("" == probabilities) {
            errorCondition(recall = weibullProbabilities.ipsur, 
                message = gettextRcmdr("No values specified."))
            return()
        }
        shape <- tclvalue(shapeVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (shape == "") {
            errorCondition(recall = weibullProbabilities.ipsur, 
                message = gettextRcmdr("The shape parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = weibullProbabilities.ipsur, 
                message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        doItAndPrint(paste("pweibull(c(", probabilities, "), shape=", 
            shape, ", scale=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pweibull")
    tkgrid(tklabel(top, text = gettextRcmdr("Variable value(s)")), 
        probabilitiesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape ")), shapeEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale ")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(probabilitiesEntry, sticky = "w")
    tkgrid.configure(shapeEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 1, focus = probabilitiesEntry)
}


`weibullQuantiles.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Weibull Quantiles"))
    quantilesVar <- tclVar("")
    quantilesEntry <- tkentry(top, width = "30", textvariable = quantilesVar)
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    scale1Var <- tclVar("1")
    scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
    tailVar <- tclVar("lower")
    lowerTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "lower")
    upperTailButton <- tkradiobutton(top, variable = tailVar, 
        value = "upper")
    onOK <- function() {
        closeDialog()
        quantiles <- gsub(" ", ",", tclvalue(quantilesVar))
        if ("" == quantiles) {
            errorCondition(recall = weibullQuantiles.ipsur, message = gettextRcmdr("No probabilities specified."))
            return()
        }
        shape <- tclvalue(shapeVar)
        scale1 <- tclvalue(scale1Var)
        tail <- tclvalue(tailVar)
        if (shape == "") {
            errorCondition(recall = weibullQuantiles.ipsur, message = gettextRcmdr("The shape parameter was not specified."))
            return()
        }
        if (scale1 == "") {
            errorCondition(recall = weibullQuantiles.ipsur, message = gettextRcmdr("The scale parameter was not specified."))
            return()
        }
        doItAndPrint(paste("qweibull(c(", quantiles, "), shape=", 
            shape, ", scale=", scale1, ", lower.tail=", tail == 
                "lower", ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qweibull")
    tkgrid(tklabel(top, text = gettextRcmdr("Probabilities")), 
        quantilesEntry, sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("shape ")), shapeEntry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("scale ")), scale1Entry, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Lower tail")), lowerTailButton, 
        sticky = "e")
    tkgrid(tklabel(top, text = gettextRcmdr("Upper tail")), upperTailButton, 
        sticky = "e")
    tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
    tkgrid.configure(quantilesEntry, sticky = "w")
    tkgrid.configure(shapeEntry, sticky = "w")
    tkgrid.configure(scale1Entry, sticky = "w")
    tkgrid.configure(lowerTailButton, sticky = "w")
    tkgrid.configure(upperTailButton, sticky = "w")
    dialogSuffix(rows = 6, columns = 2, focus = quantilesEntry)
}
