Regr1Plot <- function() {
    initializeDialog(title=gettextRcmdr("Squared Residuals"))
    .numeric <- Numeric()
    variablesFrame <- tkframe(top)
    xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("x-variable (pick one)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("y-variable (pick one)"))
    gBox <- variableListBox(variablesFrame, Factors(),
                            title=gettextRcmdr("group-variable (pick zero or one)"),
                            selectmode="single", initialSelection=-1)
    workingFrame <- tkframe(top)
    woFrame <- tkframe(workingFrame)
    checkBoxes(window=woFrame,
               frame="optionsFrame",
               boxes=c("jitterX", "points.yhat"),
               initialValues=c(0, 1),
               labels=gettextRcmdr(c("Jitter x-variable", "Display Y.hat points")))
    radioButtons(window=workingFrame,
                 name="residuals",
                 buttons=c("square",
                   "line",
                   "none"),
                 values=c("'square'", "'line'", "FALSE"),
                 initialValue="'square'",
                 labels=gettextRcmdr(c
                   ("residuals---squares",
                    "residuals---straight line",
                    "point estimate only")),
                 title=gettextRcmdr("Residual Display"))
    radioButtons(window=workingFrame,
                 name="model",
                 buttons=c("linear", "active"),
                 values=c("0", "1"),
                 initialValue="0",
                 labels=gettextRcmdr(c
                   ("default linear model",
                    "Active model")),
                 title=gettextRcmdr("Model"))
    subsetBox()    
    labelsFrame <- tkframe(top)
    xlabVar <- tclVar(gettextRcmdr("<auto>"))
    ylabVar <- tclVar(gettextRcmdr("<auto>"))
    xlabFrame <- tkframe(labelsFrame)
    xlabEntry <- tkentry(xlabFrame, width="25", textvariable=xlabVar)
    xlabScroll <- tkscrollbar(xlabFrame, orient="horizontal",
        repeatinterval=5, command=function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand=function(...) tkset(xlabScroll, ...))
    tkgrid(tklabel(xlabFrame, text=gettextRcmdr("x-axis label"), fg="blue"), sticky="w")
    tkgrid(xlabEntry, sticky="w")
    tkgrid(xlabScroll, sticky="ew")
    ylabFrame <- tkframe(labelsFrame)
    ylabEntry <- tkentry(ylabFrame, width="25", textvariable=ylabVar)
    ylabScroll <- tkscrollbar(ylabFrame, orient="horizontal",
        repeatinterval=5, command=function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand=function(...) tkset(ylabScroll, ...))
    tkgrid(tklabel(ylabFrame, text=gettextRcmdr("y-axis label"), fg="blue"), sticky="w")
    tkgrid(ylabEntry, sticky="w")
    tkgrid(ylabScroll, sticky="ew")
    tkgrid(xlabFrame, tklabel(labelsFrame, text="     "), ylabFrame, sticky="w")    
    parFrame <- tkframe(top) 


    limitsFrame <- tkframe(top)
    xminVar <- tclVar("")
    xminEntry <- tkentry(limitsFrame, width="6", textvariable=xminVar)
    xmaxVar <- tclVar("")
    xmaxEntry <- tkentry(limitsFrame, width="6", textvariable=xmaxVar)
    yminVar <- tclVar("")
    yminEntry <- tkentry(limitsFrame, width="6", textvariable=yminVar)
    ymaxVar <- tclVar("")
    ymaxEntry <- tkentry(limitsFrame, width="6", textvariable=ymaxVar)
   
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        g <- getSelection(gBox)

        xmin <- as.numeric(tclvalue(xminVar))
        xmax <- as.numeric(tclvalue(xmaxVar))
        ymin <- as.numeric(tclvalue(yminVar))
        ymax <- as.numeric(tclvalue(ymaxVar))

        number.xna <- is.na(xmin) + is.na(xmax)
        if (number.xna==1) {
            errorCondition(recall=Regr1Plot,
                           message=gettextRcmdr("Both or neither xlim values must be numbers."))
            return()
          }
        xlim.command <- ifelse(number.xna==0, deparse(c(xmin, xmax)), "")

        number.yna <- is.na(ymin) + is.na(ymax)
        if (number.yna==1) {
            errorCondition(recall=Regr1Plot,
                           message=gettextRcmdr("Both or neither ylim values must be numbers."))
            return()
          }
         ylim.command <- ifelse(number.yna==0, deparse(c(ymin, ymax)), "")



        closeDialog()
        if (length(x) == 0 || length(y) == 0){
            errorCondition(recall=Regr1Plot, message=gettextRcmdr("You must select two variables"))
            return()
            }
        if (x == y) {
            errorCondition(recall=Regr1Plot, message=gettextRcmdr("x and y variables must be different"))
            return()
            }
        .activeDataSet <- ActiveDataSet()
        Dx <- paste(.activeDataSet, "$", x, sep="")
        Dy <- paste(.activeDataSet, "$", y, sep="")
        Dg <- paste(.activeDataSet, "$", g, sep="")
        if ("1" == tclvalue(jitterXVariable)) Dx <- paste("jitter(", Dx, ")", sep="")
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" 
            else paste(", subset=", subset, sep="")

        points.yhat <- ("1" == tclvalue(points.yhatVariable))
        points.yhat.command <- ifelse(points.yhat, "", ", points.yhat=FALSE")

        xlab <- trim.blanks(tclvalue(xlabVar))
        if(xlab == gettextRcmdr("<auto>")) xlab <- x
        xlab <- paste(', xlab="', xlab, '"', sep="")
        ylab <- trim.blanks(tclvalue(ylabVar))
        if(ylab == gettextRcmdr("<auto>")) ylab <- y
        ylab <- paste(', ylab="', ylab, '"', sep="")

        ADS.x <- get(.activeDataSet)[[x]]
        ADS.y <- get(.activeDataSet)[[y]]
        ADS.lm <- lm(ADS.y ~ ADS.x)
        ADS.ylim <- range(ADS.y, predict(ADS.lm))

        aspect.x.y <- diff(range(ADS.x)) / diff(range(ADS.y))
        ADS.xlim <- range(ADS.x, ADS.x + resid(ADS.lm)*aspect.x.y)
        
        if(xlim.command == "") xlim.command <- deparse(ADS.xlim)
        xlim <- paste(', xlim=', xlim.command, sep="")
        if(ylim.command == "") ylim.command <- deparse(ADS.ylim)
        ylim <- paste(', ylim=', ylim.command, sep="")

        model <- if (0==tclvalue(modelVariable)) ""
        else {
          .activeModel <- ActiveModel()
          if (is.null(.activeModel)) {
            errorCondition(recall=Regr1Plot,
                           message=gettextRcmdr("No Active model"))
            return()
          }
          paste(', model=', .activeModel, sep="")
        }
        main <- if (model == "")
          paste('Residuals from model: ', y, ' ~ ', x)
        else
          paste('Residuals from model: ', deparse(as.formula(get(ActiveModel()))))
        
        if (length(g) == 0) {
          doItAndPrint(paste("regr1.plot(", "x=", Dx, ", y=", Dy,
                             ", resid.plot=", tclvalue(residualsVariable),
                             xlab, ylab, ", cex=1.3",
                             xlim, ylim, model,
                             subset,
                             points.yhat.command,
                             ", main='", main, "')", sep=""))
        }
        else {
          doItAndPrint(paste("regr1.plot(", "x=", Dx, ", y=", Dy,
                             ", resid.plot=", tclvalue(residualsVariable),
                             ", col=", "match(", Dg, ", levels(", Dg, "))+1",
                             xlab, ylab, ", cex=1.3",
                             xlim, ylim, model,
                             subset,
                             points.yhat.command,
                             ", main='", main, "')", sep=""))
        }
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="regr1.plot")
    tkgrid(getFrame(xBox), getFrame(yBox),
           getFrame(gBox), columnspan=1, sticky="nw") 
    tkgrid(variablesFrame, sticky="w")   
    tkgrid(subsetFrame, sticky="w")

    tkgrid(tklabel(woFrame, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(modelFrame, residualsFrame,
           woFrame, columnspan=1, sticky="nw") 
    tkgrid(workingFrame, sticky="w")   
    
    limitsnamesFrame <- tkframe(top)
    tkgrid(tklabel(limitsnamesFrame,
                   text=gettextRcmdr("xlim and ylim (try defaults first)"),
                   fg="blue"), sticky="w")
    tkgrid(tklabel(limitsnamesFrame,
                   text="xmin       xmax       ymin       ymax"), sticky="w")
    tkgrid(limitsnamesFrame, sticky="w")
    tkgrid(xminEntry, xmaxEntry, yminEntry, ymaxEntry)
    tkgrid(limitsFrame, sticky="w")
    tkgrid(labelsFrame, sticky="w")
    tkgrid(tklabel(top, text=" "))    
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=8, columns=2)
    }

## source("~/HH-R.package/RcmdrPlugin.HH/R/Regr1Plot.R")
