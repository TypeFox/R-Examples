# Last modified Feb 14, 2008
# simulations optimized by Tyler Drombosky 2007

`betaDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from beta population"))
    dsname <- tclVar(gettextRcmdr("BetaSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    shape1Var <- tclVar("1")
    shape1Entry <- tkentry(top, width = "6", textvariable = shape1Var)
    shape2Var <- tclVar("1")
    shape2Entry <- tkentry(top, width = "6", textvariable = shape2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = betaDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = betaDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                betaDistributionSamples.ipsur()
                return()
            }
        }
        shape1 <- tclvalue(shape1Var)
        shape2 <- tclvalue(shape2Var)
        ncp <- tclvalue(ncpVar)
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = betaDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = betaDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rbeta(", 
            samples, "*", n, ", shape1=", shape1, ", shape2=", 
            shape2, ", ncp=", ncp, "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rbeta")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("shape1 ")), shape1Entry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("shape2 ")), shape2Entry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = shape1Entry)
}


`binomialDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from binomial population"))
    dsname <- tclVar(gettextRcmdr("BinomialSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = binomialDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = binomialDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                binomialDistributionSamples.ipsur()
                return()
            }
        }
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = binomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        trials <- tclvalue(trialsVar)
        if (trials == "") {
            errorCondition(recall = binomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of trials not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = binomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = binomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rbinom(", 
            samples, "*", n, ", size=", trials, ", prob=", prob, 
            "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("size (number of trials)")), 
        trialsEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (probability of success)")), 
        probEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = trialsEntry)
}


`CauchyDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from Cauchy population"))
    dsname <- tclVar(gettextRcmdr("CauchySamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    sVar <- tclVar("1")
    sEntry <- tkentry(top, width = "6", textvariable = sVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = CauchyDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = CauchyDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                CauchyDistributionSamples.ipsur()
                return()
            }
        }
        location <- tclvalue(locationVar)
        s <- tclvalue(sVar)
        if (s == "") {
            errorCondition(recall = CauchyDistributionSamples.ipsur, 
                message = gettextRcmdr("Scale not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = CauchyDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = CauchyDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rcauchy(", 
            samples, "*", n, ", location=", location, ", scale=", 
            s, "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rcauchy")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("scale")), sEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = locationEntry)
}


`chisquareDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from chi-squared population"))
    dsname <- tclVar(gettextRcmdr("ChisquareSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = chisquareDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = chisquareDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                chisquareDistributionSamples.ipsur()
                return()
            }
        }
        df <- tclvalue(dfVar)
        ncp <- tclvalue(ncpVar)
        if (df == "") {
            errorCondition(recall = chisquareDistributionSamples.ipsur, 
                message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        if (ncp == "") {
            errorCondition(recall = chisquareDistributionSamples.ipsur, 
                message = gettextRcmdr("Noncentrality parameter not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = chisquareDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = chisquareDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rchisq(", 
            samples, "*", n, ", df=", df, ", ncp=", ncp, "), ncol=", 
            n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rchisq")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 10, columns = 2, focus = dfEntry)
}


`disunifDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from discrete uniform population"))
    dsname <- tclVar(gettextRcmdr("DisUnifSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    from1Var <- tclVar("1")
    from1Entry <- tkentry(top, width = "6", textvariable = from1Var)
    to1Var <- tclVar("10")
    to1Entry <- tkentry(top, width = "6", textvariable = to1Var)
    by1Var <- tclVar("1")
    by1Entry <- tkentry(top, width = "6", textvariable = by1Var)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                disunifDistributionSamples.ipsur()
                return()
            }
        }
        from1 <- tclvalue(from1Var)
        to1 <- tclvalue(to1Var)
        by1 <- tclvalue(by1Var)
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (from1 == "") {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = gettextRcmdr("Parameter 'from' not specified."))
            return()
        }
        if (to1 == "") {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = gettextRcmdr("Parameter 'to' not specified."))
            return()
        }
        if (by1 == "") {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = gettextRcmdr("Parameter 'by' not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = disunifDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste("support <- seq(", from1, ", ", to1, ", by=", by1, 
            ")", sep = "")
        justDoIt(command)
        command <- paste(dsnameValue, " <- as.data.frame(matrix(sample(support, size=", 
            samples, "*", n, ", replace=TRUE), ncol=", n, "))", 
            sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "sample")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 2, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("from (lower limit)")), 
        from1Entry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("to (upper limit)")), 
        to1Entry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("by (step size)")), 
        by1Entry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = from1Entry)
}


`exponentialDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from exponential population"))
    dsname <- tclVar(gettextRcmdr("ExponentialSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    rateVar <- tclVar("1")
    rateEntry <- tkentry(top, width = "6", textvariable = rateVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = exponentialDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = exponentialDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                exponentialDistributionSamples.ipsur()
                return()
            }
        }
        rate <- tclvalue(rateVar)
        if (rate == "") {
            errorCondition(recall = exponentialDistributionSamples.ipsur, 
                message = gettextRcmdr("Rate not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = exponentialDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = exponentialDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rexp(", 
            samples, "*", n, ", rate=", rate, "), ncol=", n, 
            "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rexp")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameter:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("rate (of arrivals in unit time)")), 
        rateEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 10, columns = 2, focus = rateEntry)
}


`FDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from F population"))
    dsname <- tclVar(gettextRcmdr("FSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    df1Var <- tclVar("1")
    df2Var <- tclVar("1")
    df1Entry <- tkentry(top, width = "6", textvariable = df1Var)
    df2Entry <- tkentry(top, width = "6", textvariable = df2Var)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                FDistributionSamples.ipsur()
                return()
            }
        }
        df1 <- tclvalue(df1Var)
        df2 <- tclvalue(df2Var)
        ncp <- tclvalue(ncpVar)
        if (df1 == "") {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = gettextRcmdr("Numerator degrees of freedom not specified."))
            return()
        }
        if (df2 == "") {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = gettextRcmdr("Denominator degrees of freedom not specified."))
            return()
        }
        if (ncp == "") {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = gettextRcmdr("Noncentrality parameter not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = FDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rf(", 
            samples, "*", n, ", df1=", df1, ", df2=", df2, ", ncp=", 
            ncp, "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rf")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("df1 (numerator degrees of freedom)")), 
        df1Entry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("df2 (denominator degrees of freedom)")), 
        df2Entry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = df1Entry)
}


`gammaDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from gamma population"))
    dsname <- tclVar(gettextRcmdr("GammaSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    sVar <- tclVar("1")
    sEntry <- tkentry(top, width = "6", textvariable = sVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = gammaDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = gammaDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                gammaDistributionSamples.ipsur()
                return()
            }
        }
        shape <- tclvalue(shapeVar)
        if (shape == "") {
            errorCondition(recall = gammaDistributionSamples.ipsur, 
                message = gettextRcmdr("Shape not specified."))
            return()
        }
        s <- tclvalue(sVar)
        if (s == "") {
            errorCondition(recall = gammaDistributionSamples.ipsur, 
                message = gettextRcmdr("Rate not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = gammaDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = gammaDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rgamma(", 
            samples, "*", n, ", shape=", shape, ", rate=", s, 
            "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rgamma")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("shape")), shapeEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("rate (= 1/scale)")), 
        sEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = shapeEntry)
}


`geomDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from geometric population"))
    dsname <- tclVar(gettextRcmdr("GeometricSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = geomDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = geomDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                geomDistributionSamples.ipsur()
                return()
            }
        }
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = geomDistributionSamples.ipsur, 
                message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = geomDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = geomDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rgeom(", 
            samples, "*", n, ", prob=", prob, "), ncol=", n, 
            "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rgeom")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (probability of success)")), 
        probEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 10, columns = 2, focus = probEntry)
}


`hyperDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from hypergeometric population"))
    dsname <- tclVar(gettextRcmdr("HypergeometricSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    mVar <- tclVar("1")
    mEntry <- tkentry(top, width = "6", textvariable = mVar)
    nVar <- tclVar("1")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    kVar <- tclVar("1")
    kEntry <- tkentry(top, width = "6", textvariable = kVar)
    sampleSizeVar <- tclVar("10")
    sampleSizeEntry <- tkentry(top, width = "6", textvariable = sampleSizeVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                hyperDistributionSamples.ipsur()
                return()
            }
        }
        m <- tclvalue(mVar)
        n <- tclvalue(nVar)
        k <- tclvalue(kVar)
        if (m == "") {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = gettextRcmdr("The m parameter was not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = gettextRcmdr("The n parameter was not specified."))
            return()
        }
        if (k == "") {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = gettextRcmdr("The k parameter was not specified."))
            return()
        }
        sampleSize <- tclvalue(sampleSizeVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (sampleSize == "") {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = hyperDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rhyper(", 
            samples, "*", sampleSize, ", n=", n, ", m=", m, ", k=", 
            k, "), ncol=", sampleSize, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (sampleSize == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            sampleSize, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", sampleSize, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", sampleSize, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", sampleSize, "]), 1, sd)", 
                sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rhyper")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        sampleSizeEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("m (number of white balls in the urn)")), 
        mEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("n (number of black balls in the urn)")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("k (number of balls drawn from the urn)")), 
        kEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 12, columns = 2, focus = mEntry)
}


`logisticDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from logistic population"))
    dsname <- tclVar(gettextRcmdr("LogisticSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    locationVar <- tclVar("0")
    locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
    sVar <- tclVar("1")
    sEntry <- tkentry(top, width = "6", textvariable = sVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = logisticDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = logisticDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                logisticDistributionSamples.ipsur()
                return()
            }
        }
        location <- tclvalue(locationVar)
        s <- tclvalue(sVar)
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = logisticDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = logisticDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rlogis(", 
            samples, "*", n, ", location=", location, ", scale=", 
            s, "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rlogis")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("scale")), sEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = locationEntry)
}



`lognormalDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from log-normal population"))
    dsname <- tclVar(gettextRcmdr("LogNormalSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    meanlogVar <- tclVar("0")
    meanlogEntry <- tkentry(top, width = "6", textvariable = meanlogVar)
    sdlogVar <- tclVar("1")
    sdlogEntry <- tkentry(top, width = "6", textvariable = sdlogVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = lognormalDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = lognormalDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                lognormalDistributionSamples.ipsur()
                return()
            }
        }
        meanlog <- tclvalue(meanlogVar)
        sdlog <- tclvalue(sdlogVar)
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = lognormalDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = lognormalDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rlnorm(", 
            samples, "*", n, ", meanlog=", meanlog, ", sdlog=", 
            sdlog, "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rlnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("meanlog (mean of dist'n on log scale)")), 
        meanlogEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("sdlog (std dev of dist'n on log scale)")), 
        sdlogEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = meanlogEntry)
}



`negbinomialDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from negative binomial population"))
    dsname <- tclVar(gettextRcmdr("NegativeBinomialSamples"))
    entryDsname <- tkentry(top, width = "22", textvariable = dsname)
    probVar <- tclVar("0.5")
    probEntry <- tkentry(top, width = "6", textvariable = probVar)
    trialsVar <- tclVar("1")
    trialsEntry <- tkentry(top, width = "6", textvariable = trialsVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = negbinomialDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = negbinomialDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                negbinomialDistributionSamples.ipsur()
                return()
            }
        }
        prob <- tclvalue(probVar)
        if (prob == "") {
            errorCondition(recall = negbinomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Probability of success not specified."))
            return()
        }
        trials <- tclvalue(trialsVar)
        if (trials == "") {
            errorCondition(recall = negbinomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Target number of successes not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = negbinomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = negbinomialDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rnbinom(", 
            samples, "*", n, ", size=", trials, ", prob=", prob, 
            "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rnbinom")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("size (target number of successes)")), 
        trialsEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("prob (probability of success)")), 
        probEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = trialsEntry)
}


`normalDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from Normal population"))
    dsname <- tclVar(gettextRcmdr("NormalSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    muVar <- tclVar("0")
    muEntry <- tkentry(top, width = "6", textvariable = muVar)
    sigmaVar <- tclVar("1")
    sigmaEntry <- tkentry(top, width = "6", textvariable = sigmaVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = normalDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = normalDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                normalDistributionSamples.ipsur()
                return()
            }
        }
        mu <- tclvalue(muVar)
        sigma <- tclvalue(sigmaVar)
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (sigma == "") {
            errorCondition(recall = normalDistributionSamples.ipsur, 
                message = gettextRcmdr("Standard deviation not specified."))
            return()
        }
        if (mu == "") {
            errorCondition(recall = normalDistributionSamples.ipsur, 
                message = gettextRcmdr("Mean not specified."))
            return()
        }
        if (n == "") {
            errorCondition(recall = normalDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = normalDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rnorm(", 
            samples, "*", n, ", mean=", mu, ", sd=", sigma, "), ncol=", 
            n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rnorm")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 2, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("mu (mean)")), muEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("sigma (standard deviation)")), 
        sigmaEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = muEntry)
}


`PoissonDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from Poisson population"))
    dsname <- tclVar(gettextRcmdr("PoissonSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    meanVar <- tclVar("1")
    meanEntry <- tkentry(top, width = "6", textvariable = meanVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = PoissonDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = PoissonDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                PoissonDistributionSamples.ipsur()
                return()
            }
        }
        mean <- tclvalue(meanVar)
        if (mean == "") {
            errorCondition(recall = PoissonDistributionPlot, 
                message = gettextRcmdr("Mean not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = PoissonDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = PoissonDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rpois(", 
            samples, "*", n, ", lambda=", mean, "), ncol=", n, 
            "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rpois")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameter:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("lambda (mean)")), 
        meanEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 10, columns = 2, focus = meanEntry)
}


`tDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from Student's t population"))
    dsname <- tclVar(gettextRcmdr("tSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    dfVar <- tclVar("1")
    dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
    ncpVar <- tclVar("0")
    ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = tDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = tDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                tDistributionSamples.ipsur()
                return()
            }
        }
        df <- tclvalue(dfVar)
        ncp <- tclvalue(ncpVar)
        if (df == "") {
            errorCondition(recall = tDistributionSamples.ipsur, 
                message = gettextRcmdr("Degrees of freedom not specified."))
            return()
        }
        if (ncp == "") {
            errorCondition(recall = tDistributionSamples.ipsur, 
                message = gettextRcmdr("Noncentrality parameter not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = tDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = tDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rt(", 
            samples, "*", n, ", df=", df, ", ncp=", ncp, "), ncol=", 
            n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rt")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
        dfEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
        ncpEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 10, columns = 2, focus = dfEntry)
}


`uniformDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from uniform population"))
    dsname <- tclVar(gettextRcmdr("UniformSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    minVar <- tclVar("0")
    maxVar <- tclVar("1")
    minEntry <- tkentry(top, width = "6", textvariable = minVar)
    maxEntry <- tkentry(top, width = "6", textvariable = maxVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = uniformDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = uniformDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                uniformDistributionSamples.ipsur()
                return()
            }
        }
        minValue <- tclvalue(minVar)
        maxValue <- tclvalue(maxVar)
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = uniformDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = uniformDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(runif(", 
            samples, "*", n, ", min=", minValue, ", max=", maxValue, 
            "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "runif")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("min (lower limit of the distribution)")), 
        minEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("max (upper limit of the distribution)")), 
        maxEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = minEntry)
}


`WeibullDistributionSamples.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Sample from weibull population"))
    dsname <- tclVar(gettextRcmdr("WeibullSamples"))
    entryDsname <- tkentry(top, width = "20", textvariable = dsname)
    shapeVar <- tclVar("1")
    shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
    sVar <- tclVar("1")
    sEntry <- tkentry(top, width = "6", textvariable = sVar)
    nVar <- tclVar("10")
    nEntry <- tkentry(top, width = "6", textvariable = nVar)
    samplesVar <- tclVar("100")
    samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sum", 
        "sd"), initialValues = c("1", "0", "0"), labels = gettextRcmdr(c("Sample means", 
        "Sample sums", "Sample standard deviations")))
    otherstatFrame <- tkframe(top)
    otherstatVariable <- tclVar("0")
    otherstatCheckBox <- tkcheckbutton(otherstatFrame, variable = otherstatVariable)
    otherstat <- tclVar("median")
    otherstat2 <- tclVar("IQR")
    otherstat3 <- tclVar("< function >")
    otherstatEntry <- tkentry(otherstatFrame, width = "12", textvariable = otherstat)
    otherstat2Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat2)
    otherstat3Entry <- tkentry(otherstatFrame, width = "12", 
        textvariable = otherstat3)
    savesimVariable <- tclVar("1")
    checkBoxes(frame = "discardcheckBoxFrame", boxes = c("discard"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Discard the original observations")))
    onOK <- function() {
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall = WeibullDistributionSamples.ipsur, 
                message = gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall = WeibullDistributionSamples.ipsur, 
                message = paste("\"", dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
                  sep = ""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))) {
                WeibullDistributionSamples.ipsur()
                return()
            }
        }
        shape <- tclvalue(shapeVar)
        if (shape == "") {
            errorCondition(recall = WeibullDistributionSamples.ipsur, 
                message = gettextRcmdr("Shape not specified."))
            return()
        }
        s <- tclvalue(sVar)
        if (s == "") {
            errorCondition(recall = WeibullDistributionSamples.ipsur, 
                message = gettextRcmdr("Scale not specified."))
            return()
        }
        n <- tclvalue(nVar)
        samples <- tclvalue(samplesVar)
        otherst <- tclvalue(otherstat)
        otherst2 <- tclvalue(otherstat2)
        otherst3 <- tclvalue(otherstat3)
        if (n == "") {
            errorCondition(recall = WeibullDistributionSamples.ipsur, 
                message = gettextRcmdr("Sample size not specified."))
            return()
        }
        if (samples == "") {
            errorCondition(recall = WeibullDistributionSamples.ipsur, 
                message = gettextRcmdr("Number of samples not specified."))
            return()
        }
        command <- paste(dsnameValue, " <- as.data.frame(matrix(rweibull(", 
            samples, "*", n, ", shape=", shape, ", scale=", s, 
            "), ncol=", n, "))", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (samples == 1) 
            paste("rownames(", dsnameValue, ") <- \"sample\"", 
                sep = "")
        else paste("rownames(", dsnameValue, ") <- paste(\"sample\", 1:", 
            samples, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        command <- if (n == 1) 
            paste("colnames(", dsnameValue, ") <- \"obs\"", sep = "")
        else paste("colnames(", dsnameValue, ") <- paste(\"obs\", 1:", 
            n, ", sep=\"\")", sep = "")
        justDoIt(command)
        logger(command)
        if (tclvalue(meanVariable) == "1") {
            command <- paste(dsnameValue, "$mean <- rowMeans(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sumVariable) == "1") {
            command <- paste(dsnameValue, "$sum <- rowSums(as.matrix(", 
                dsnameValue, "[,1:", n, "]))", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(sdVariable) == "1") {
            command <- paste(dsnameValue, "$sd <- apply(as.matrix(", 
                dsnameValue, "[,1:", n, "]), 1, sd)", sep = "")
            justDoIt(command)
            logger(command)
        }
        if (tclvalue(otherstatVariable) == "1") {
            if (otherst != "") {
                command <- paste(dsnameValue, "$", otherst, " <- apply(as.matrix(", 
                  dsnameValue, "[,1:", n, "]), 1, ", otherst, 
                  ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst2 != "") {
                command <- paste(dsnameValue, "$", otherst2, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst2, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
            if (otherst3 != "< function >" & otherst3 != "") {
                command <- paste(dsnameValue, "$", otherst3, 
                  " <- apply(as.matrix(", dsnameValue, "[,1:", 
                  n, "]), 1, ", otherst3, ")", sep = "")
                justDoIt(command)
                logger(command)
            }
        }
        if (tclvalue(discardVariable) == "1") {
            if (n == 1) {
                command <- paste(dsnameValue, "$obs <- NULL", 
                  sep = "")
                justDoIt(command)
            }
            else {
                variables = paste("obs", 1:as.numeric(n), sep = "")
                for (k in 1:as.numeric(n)) {
                  eval(parse(text = paste(dsnameValue, "$obs", 
                    k, "<- NULL", sep = "")), envir = .GlobalEnv)
                }
            }
            logger(paste("The original observations were discarded.", 
                sep = ""))
        }
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rweibull")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:"), 
        fg = "blue"), entryDsname, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("")), columnspan = 4, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Dimensions:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (rows) ")), 
        samplesEntry, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Sample size (columns) ")), 
        nEntry, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
        columnspan = 4, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("shape")), shapeEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("scale")), sEntry, 
        sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(tklabel(top, text = gettextRcmdr("Add to Data Set:"), 
        fg = "blue"), sticky = "w")
    tkgrid(checkBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("Other sample statistic(s)      ")), 
        otherstatCheckBox, tklabel(otherstatFrame, text = gettextRcmdr(" Specify:"), 
            fg = "blue"), otherstatEntry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat2Entry, sticky = "w")
    tkgrid(tklabel(otherstatFrame, text = gettextRcmdr("")), 
        tklabel(otherstatFrame, text = gettextRcmdr("")), tklabel(otherstatFrame, 
            text = gettextRcmdr("")), otherstat3Entry, sticky = "w")
    tkgrid(otherstatFrame, columnspan = 2, sticky = "w")
    tkgrid(discardcheckBoxFrame, columnspan = 2, sticky = "w")
    tkgrid(tklabel(top, text = ""))
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 11, columns = 2, focus = shapeEntry)
}
