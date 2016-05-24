"Scatter3DDialog.HH" <-
function(){
    initializeDialog(title=gettextRcmdr("3D Scatterplot"))
    variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Explanatory variables (pick two)"), selectmode="multiple",
        initialSelection=NULL, listHeight=7)
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Response variable (pick one)"), listHeight=7)
    surfacesFrame <- tkframe(top)
    identifyPoints <- tclVar("0")
    identifyPointsCheckBox <- tkcheckbutton(surfacesFrame, variable=identifyPoints)
    gridLines <- tclVar("1")
    gridLinesCheckBox <- tkcheckbutton(surfacesFrame, variable=gridLines)
    rglOpen <- tclVar("0")
    rglOpenCheckBox <- tkcheckbutton(surfacesFrame, variable=rglOpen)
    linearLSSurface <- tclVar("1")
    linearLSCheckBox <- tkcheckbutton(surfacesFrame, variable=linearLSSurface)
    betaMultiplierVariable <- tclVar(gettextRcmdr("1.0"))
    betaMultiplierField <- tkentry(surfacesFrame, width="6", textvariable=betaMultiplierVariable)
    quadLSSurface <- tclVar("0")
    quadLSCheckBox <- tkcheckbutton(surfacesFrame, variable=quadLSSurface)
    nonparSurface <- tclVar("0")
    nonparCheckBox <- tkcheckbutton(surfacesFrame, variable=nonparSurface)
    dfNonparVariable <- tclVar(gettextRcmdr("<auto>"))
    dfNonparField <- tkentry(surfacesFrame, width="6", textvariable=dfNonparVariable)
    additiveSurface <- tclVar("0")
    additiveCheckBox <- tkcheckbutton(surfacesFrame, variable=additiveSurface)
    dfAddVariable <- tclVar(gettextRcmdr("<auto>"))
    dfAddField <- tkentry(surfacesFrame, width="6", textvariable=dfAddVariable)
    bgFrame <- tkframe(top)
    bgVariable <- tclVar("white")
    whiteButton <- tkradiobutton(bgFrame, variable=bgVariable, value="white")
    blackButton <- tkradiobutton(bgFrame, variable=bgVariable, value="black")
    residFrame <- tkframe(top)
    residVariable <- tclVar("line")
    lineButton <- tkradiobutton(residFrame, variable=residVariable, value="line")
    squareButton <- tkradiobutton(residFrame, variable=residVariable, value="square")
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (length(y) == 0) {
            errorCondition(recall=Scatter3DDialog.HH, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (2 != length(x)) {
            errorCondition(recall=Scatter3DDialog.HH, message=gettextRcmdr("You must select 2 explanatory variables."))
            return()
            }
        if (is.element(y, x)) {
            errorCondition(recall=Scatter3DDialog.HH, message=gettextRcmdr("Response and explanatory variables must be different."))
            return()
            }
        grid <- if (tclvalue(gridLines) == 1) "TRUE" else "FALSE"
        lin <- if(tclvalue(linearLSSurface) == 1) '"linear"'
        quad <- if(tclvalue(quadLSSurface) == 1) '"quadratic"'
        nonpar <- if (tclvalue(nonparSurface) == 1) '"smooth"'
        additive <- if (tclvalue(additiveSurface) == 1) '"additive"'
        surfaces <- c(lin, quad, nonpar, additive)
        nsurfaces <- length(surfaces)
        opts <- options(warn=-1)
        betaMultiplier <- tclvalue(betaMultiplierVariable)
        betaMultiplier <- if (betaMultiplier == gettextRcmdr("1.0")) "" else paste(", coef.ratio=", as.numeric(betaMultiplier), sep="")
        dfNonpar <- tclvalue(dfNonparVariable)
        dfNonpar <- if (dfNonpar == gettextRcmdr("<auto>")) "" else paste(", df.smooth=", as.numeric(dfNonpar), sep="")
        dfAdd <- tclvalue(dfAddVariable)
        dfAdd <- if (dfAdd == gettextRcmdr("<auto>")) "" else paste(", df.additive=", as.numeric(dfAdd), sep="")
        options(opts)
        fit <- if (nsurfaces == 0) ", surface=FALSE"
            else if (nsurfaces == 1) paste(", fit=", surfaces, sep="")
            else paste(", fit=c(", paste(surfaces, collapse=","), ")", sep="")
        bg <- tclvalue(bgVariable)
        resid <- tclvalue(residVariable)
        .activeDataSet <- ActiveDataSet()
        if (.groups != FALSE){
            groups <- paste(", groups=", .activeDataSet, "$", .groups, sep="")
            parallel <- paste(", parallel=", .linesByGroup, sep="")
            }
        else groups <- parallel <- ""
        if (tclvalue(rglOpen) == 1) doItAndPrint("rgl.open()")
        command <- paste("scatter3dHH(", .activeDataSet, "$", x[1], ", ",
            .activeDataSet, "$", y, ", ", .activeDataSet, "$", x[2], fit, dfNonpar,
            dfAdd, groups, parallel, ', bg="', bg, '", grid=', grid,
                         ', squares=', resid=="square", betaMultiplier,
            ', xlab="', x[1], '", ylab="', y, '", zlab="', x[2], '")', sep="")
        doItAndPrint(command)
        putRcmdr("rgl", TRUE)
        command <- paste("Identify3d(", .activeDataSet, "$", x[1], ", ",
            .activeDataSet, "$", y, ", ", .activeDataSet, "$", x[2], groups,
            ", labels=row.names(", .activeDataSet, "))", sep="")
        putRcmdr("Identify3d", command)
        .Tcl("update")
        if (tclvalue(identifyPoints) == 1) doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
        rgl.bringtotop()
        }
    groupsBox(Scatter3DDialog.HH, plotLinesByGroup=TRUE, plotLinesByGroupsText=gettextRcmdr("Parallel regression surfaces"))
    OKCancelHelp(helpSubject="Scatter3DDialog")
    tkgrid(getFrame(yBox), tklabel(variablesFrame, text="  "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="nw")
    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Identify observations\nwith mouse")), identifyPointsCheckBox, sticky="w")
    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Show surface grid lines")), gridLinesCheckBox, sticky="w")
    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Open new 3D graphics window")), rglOpenCheckBox, sticky="w")
    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Surfaces to Fit"), fg="blue"), sticky="w")

    betaMultiplierLabel <- tklabel(surfacesFrame, text=gettextRcmdr("coef multiple = "))
    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Linear least-squares")), linearLSCheckBox, betaMultiplierLabel, betaMultiplierField, sticky="w")
    tkgrid.configure(betaMultiplierLabel, sticky="e")

    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Quadratic least-squares")), quadLSCheckBox, sticky="w")

    dfLabel <- tklabel(surfacesFrame, text=gettextRcmdr("df = "))
    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Smooth regression")), nonparCheckBox,
        dfLabel, dfNonparField, sticky="w")
    tkgrid.configure(dfLabel, sticky="e")

    tkgrid(tklabel(surfacesFrame, text=gettextRcmdr("Additive regression")), additiveCheckBox,
        tklabel(surfacesFrame, text=gettextRcmdr("df(each term) = ")), dfAddField, sticky="w")
    tkgrid(surfacesFrame, sticky="w")
    tkgrid(tklabel(bgFrame, text=gettextRcmdr("Background Color"), fg="blue"), sticky="w", columnspan=2)
    tkgrid(tklabel(bgFrame, text=gettextRcmdr("Black")), blackButton, sticky="w")
    tkgrid(tklabel(bgFrame, text=gettextRcmdr("White")), whiteButton, sticky="w")
    tkgrid(bgFrame, sticky="w")

    tkgrid(tklabel(residFrame, text=gettextRcmdr("Residuals Display"), fg="blue"), sticky="w", columnspan=2)
    tkgrid(tklabel(residFrame, text=gettextRcmdr("Line")), lineButton, sticky="w")
    tkgrid(tklabel(residFrame, text=gettextRcmdr("Square")), squareButton, sticky="w")
    tkgrid(residFrame, sticky="w")

    tkgrid(groupsFrame, sticky="w")
    tkgrid(buttonsFrame, stick="w")
    dialogSuffix(rows=5, columns=1)
    }

