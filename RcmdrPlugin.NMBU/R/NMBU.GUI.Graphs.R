# $Id: NMBU.GUI.Graphs.R 35 2014-01-10 21:17:26Z khliland $

##
## GUI functions for the Graphs menu
##



## GUI tips
#
# Usual code structure:
#    1. Intitialise window and prepare graphical elements
#    2. onOK function contianing actions to perform
#       2.1 Collect values from GUI
#       2.2 Test if combination of values is usable
#       2.3 Perform main calculations, print, update models/data etc.
#    3. Set up GUI.
#       - tkgrid() adds elements. Explicit placement and width/heigth by colum, row, columnspan and rowspan
#       - Frames with graphical elements are safer than direct placement of elements due to version conflicts.
#       - dialogSuffix() defines the final size of the grid used for elements.

#####################################
# Dotplot (Minitab) using dotPlot()
dotplotGUI <- function(){
  defaults <- list(initial.x = NULL, initialGroup=NULL, initial.stacked = 0, initial.commonscale = 1) 
  dialog.values <- getDialog("dotplot", defaults)
  initializeDialog(title=gettextRcmdr("Dotplot"))
  variablesFrame <- tkframe(top)
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Variable to be ploted (pick one)"),
                          initialSelection = varPosn (dialog.values$initial.x, "numeric"))
  initial.group <- dialog.values$initial.group
  .groups <- if (is.null(initial.group)) FALSE else initial.group
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    stacked <- tclvalue(stackedVariable)
    commonscale <- tclvalue(commonscaleVariable)
    if (length(x) != 1){
      errorCondition(recall=dotplotGUI, message=gettextRcmdr("You must select one variable."))
      return()
    }
    putDialog ("dotplot", list(initial.x = x, 
                               initial.group=if (.groups == FALSE) NULL else .groups, 
                               initial.stacked = "1"==stacked,
                               initial.commonscale = "1"==commonscale))
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    if (is.null(.groups) || .groups == FALSE) {
      doItAndPrint(paste("dotPlot(", .activeDataSet, "$", x, ", xlab='", x, "',cex=1)", sep=""))
    } else {
      ngr <- eval(parse(text=paste("length(levels(",.activeDataSet,"$",.groups,"))")))
      if(stacked=="1"){
        pch.tmp <- c(21,20,17,15,5,4,3)
        doItAndPrint(paste("pch <- c(", paste(pch.tmp[1:ngr],sep="",collapse=","),")",sep=""))
        doItAndPrint(paste("ord <- order(as.numeric(", .activeDataSet, "$", .groups, "))", sep=""))
        doItAndPrint(paste("dotPlot(", .activeDataSet, "$", x, "[ord], xlab='", x, "', pch=pch[",.activeDataSet,"$",.groups,"[ord]], cex=1)", sep=""))
        doItAndPrint(paste("legend('topright', legend=c('",paste(eval(parse(text=paste("levels(factor(",.activeDataSet,"$",.groups,"[ord]))",sep=""))),sep="",collapse="','"),"'), pch=pch)",sep=""))
        doItAndPrint("rm(list=c('pch','ord'))")
      } else {
        doItAndPrint(paste("par(mfrow=c(",ngr,",1), mar=c(4,4,1,1), cex=1)",sep=""))
        if(commonscale=="1"){
          for(i in 1:ngr){
            doItAndPrint(paste("dotPlot(", .activeDataSet, "$", x, "[", .activeDataSet,"$",.groups,"==levels(", .activeDataSet,"$",.groups,")[",i,"]],xlim = range(",.activeDataSet,"$",x,", na.rm = TRUE), xlab='", x, " (",eval(parse(text=paste("levels(", .activeDataSet,"$",.groups,")[",i,"]",sep=""))),")')", sep=""))
          }
        }else{
          for(i in 1:ngr){
            doItAndPrint(paste("dotPlot(", .activeDataSet, "$", x, "[", .activeDataSet,"$",.groups,"==levels(", .activeDataSet,"$",.groups,")[",i,"]], xlab='", x, " (",eval(parse(text=paste("levels(", .activeDataSet,"$",.groups,")[",i,"]",sep=""))),")')", sep=""))
          }
        }
        doItAndPrint(paste("par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1, cex=1)",sep=""))
      }
    }
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  groupsBox(dotPlot, initialGroup=initial.group, 
            initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
  OKCancelHelp(helpSubject="dotPlot")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw")
  tkgrid(groupsFrame, sticky = "w")
  stackedFrame <- tkframe(top)
  commonscaleFrame <- tkframe(top)
  checkBoxes(frame="stackedFrame", boxes=c("stacked"), initialValues=dialog.values$initial.stacked, labels=gettextRcmdr(c("Stack groups")))
  checkBoxes(frame="commonscaleFrame", boxes=c("commonscale"), initialValues=dialog.values$initial.commonscale, labels=gettextRcmdr(c("Use common x-scale")))
  tkgrid(stackedFrame, sticky="w")
  tkgrid(commonscaleFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=3, columns=1)
}

#####################################
# Fitted line (regression) plot
fittedLinePlot <- function(){
  initializeDialog(title=gettextRcmdr("Fitted regression plot"))
  variablesXFrame <- tkframe(top)
  variablesYFrame <- tkframe(top)
  .numeric <- Numeric()
  .activeDataSet <- ActiveDataSet()
  xBox <- variableListBox(variablesXFrame, .numeric, title=gettextRcmdr("X: regressor variable (pick one)"))
  yBox <- variableListBox(variablesYFrame, .numeric, title=gettextRcmdr("Y: response variable (pick one)"))
  comboXFrame <- tkframe(top)
  # comboYFrame <- tkframe(top)
  comboXVar <- tclVar()
  # comboYVar <- tclVar()
  # FIXME
  valuesX <- c('none', 'no intercept', 'x^2, x', 'x^3, x^2, x', 'log(x)', '1/x', 'exp(x)')
  # valuesY <- c('none', 'x^2, x', 'x^3, x^2, x', 'log(x)', '1/x', 'exp(x)')
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    if (length(x) == 0 | length(y) == 0){
      errorCondition(recall=fittedLinePlot, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y){
      errorCondition(recall=fittedLinePlot, message=gettextRcmdr("Variables must be different."))
      return()
    }
    level <- tclvalue(confidenceLevel)
    range.X <- justDoIt(paste("range(", .activeDataSet, "$", x, ")", sep=""))
    new.x   <- seq(range.X[1], range.X[2], length.out=200)
    if(tclvalue(comboXVar)==""){
      selectedTrans <- 1
    } else {
      selectedTrans <- which(valuesX==tclvalue(comboXVar))}
    if(selectedTrans==1){
      x.val <- paste(x, sep="")
    }
    if(selectedTrans==2){
      x.val <- paste(x, "-1", sep="")
    }
    if(selectedTrans==3){
      x.val <- paste("I(", x, "^2) + ", x, sep="")
    }
    if(selectedTrans==4){
      x.val <- paste("I(", x, "^3) + I(", x, "^2) + ", x, sep="")
    }
    if(selectedTrans==5){
      x.val <- paste("log(", x, ")", sep="")
    }
    if(selectedTrans==6){
      x.val <- paste("1/", x, sep="")
    }
    if(selectedTrans==7){
      x.val <- paste("exp(", x, ")", sep="")
    }
    closeDialog()
    justDoIt(paste("my.lm <- lm(", y, "~", x.val, ", data=list(", y, "=", .activeDataSet, "$", y, ", ", x, "=", .activeDataSet, "$", x, "))", sep=""))
    
    y.conf <- justDoIt(paste("predict(my.lm, data.frame(", x, "=c(", paste(new.x, sep="", collapse=","), ")), interval='conf', level=", level, ")", sep=""))
    y.pred <- justDoIt(paste("predict(my.lm, data.frame(", x, "=c(", paste(new.x, sep="", collapse=","), ")), interval='prediction', level=", level, ")", sep=""))
    
    lm.coef <- coef(my.lm)
    fit.text <- paste(y, " ~ ", format(lm.coef[1],digits=2), sep="")
    for(i in 2:length(lm.coef))
      fit.text <- paste(fit.text, " + ", format(lm.coef[i],digits=2), "*", names(lm.coef)[i], sep="")
    
    sub.text <- paste(100 * as.numeric(level), "% CI and PI, R^2=", format(summary(my.lm)[]$r.squared, digits=2), sep = "")
    justDoIt(paste("plot(", .activeDataSet, "$", x, ", ", .activeDataSet, "$", y, ", xlab='", x, "', ylab='", y, "', ylim=c(",min(y.pred),",",max(y.pred),"), main=c('",fit.text,"', '", sub.text,"'))", sep=""))
    justDoIt(paste("lines(c(", paste(new.x, sep="", collapse=","), "), c(", paste(y.conf[,1], sep="", collapse=","), "))", sep=""))
    justDoIt(paste("lines(c(", paste(new.x, sep="", collapse=","), "), c(", paste(y.conf[,2], sep="", collapse=","), "), lty=2, col='blue')", sep=""))
    justDoIt(paste("lines(c(", paste(new.x, sep="", collapse=","), "), c(", paste(y.conf[,3], sep="", collapse=","), "), lty=2, col='blue')", sep=""))
    justDoIt(paste("lines(c(", paste(new.x, sep="", collapse=","), "), c(", paste(y.pred[,2], sep="", collapse=","), "), lty=2, col='red')", sep=""))
    justDoIt(paste("lines(c(", paste(new.x, sep="", collapse=","), "), c(", paste(y.pred[,3], sep="", collapse=","), "), lty=2, col='red')", sep=""))
    ##
    
    
    # doItAndPrint(paste("CIplot(my.lm, conf.level=", level,")"))
    # doItAndPrint("rm(my.lm)")
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="lm")
  optionsFrame <- tkframe(top)
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(variablesXFrame, column=1, row=1, columnspan=1, sticky="nw")
  tkgrid(getFrame(yBox), sticky="nw")
  tkgrid(variablesYFrame, column=2, row=1, columnspan=1, sticky="nw")
  comboX <- ttkcombobox(comboXFrame, values=valuesX, textvariable=comboXVar)
  tkgrid(labelRcmdr(comboXFrame, text=gettextRcmdr("X transf.:"), fg="blue"), comboX, sticky="w")
  tkgrid(comboXFrame, sticky="w", column=1, row=2, columnspan=1)
  # comboY <- ttkcombobox(comboYFrame, values=valuesY, textvariable=comboYVar)
  # tkgrid(labelRcmdr(comboYFrame, text=gettextRcmdr("Y transf.:"), fg="blue"), comboY, sticky="w")
  # tkgrid(comboYFrame, sticky="w", column=2, row=2, columnspan=1)
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(confidenceFrame, labelRcmdr(optionsFrame, text="    "),sticky="nw")
  tkgrid(optionsFrame, column=1, row=3, columnspan=2, sticky="nw")
  tkgrid(buttonsFrame, column=1, row=4, columnspan=2, sticky="w")
  dialogSuffix(rows=4, columns=2)
}

######################################
## Customized line plot
linePlotNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Line and point plot"))
  variablesFrame <- tkframe(top)
  .numeric <- Numeric()
  .variable <- Variables()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("x variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("y variables (pick one or more)"),
                          selectmode="multiple", initialSelection=NULL)
  zBox <- variableListBox(variablesFrame, .variable, title=gettextRcmdr("Groups (optional)"))
  axisLabelVariable <- tclVar(gettextRcmdr("<use y-variable names>"))
  axisLabelFrame <- tkframe(top)
  axisLabelEntry <- ttkentry(axisLabelFrame, width="40", textvariable=axisLabelVariable)
  axisLabelScroll <- ttkscrollbar(axisLabelFrame, orient="horizontal",
                                  command=function(...) tkxview(axisLabelEntry, ...))
  tkconfigure(axisLabelEntry, xscrollcommand=function(...) tkset(axisLabelScroll, ...))
  legendFrame <- tkframe(top)
  legendVariable <- tclVar("0")
  legendCheckBox <- tkcheckbutton(legendFrame, variable=legendVariable)
  onOK <- function(){
    z <- getSelection(zBox)
    y <- getSelection(yBox)
    x <- getSelection(xBox)
    lineType <- as.character(tclvalue(lineTypeVariable))
    closeDialog()
    if (0 == length(x)) {
      errorCondition(recall=linePlotNMBU, message=gettextRcmdr("No x variable selected."))
      return()
    }
    if (0 == length(y)) {
      errorCondition(recall=linePlotNMBU, message=gettextRcmdr("No y variables selected."))
      return()
    }
    if (1 < length(y) && length(z) == 1) {
      errorCondition(recall=linePlotNMBU, message=gettextRcmdr("Only one y variable can be plotted with groups."))
      return()
    }
    if (is.element(x, y)) {
      errorCondition(recall=linePlotNMBU, message=gettextRcmdr("x and y variables must be different."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    .x <- na.omit(eval(parse(text=paste(.activeDataSet, "$", x, sep="")), envir=.GlobalEnv))
    if (!identical(order(.x), seq(along.with=.x)) && length(z)==0){
      response <- tclvalue(RcmdrTkmessageBox(message=gettextRcmdr("x-values are not in order.\nContinue?"),
                                             icon="warning", type="okcancel", default="cancel"))
      if (response == "cancel") {
        onCancel()
        return()
      }
    }
    axisLabel <- tclvalue(axisLabelVariable)
    legend <- tclvalue(legendVariable) == "1"
    if (axisLabel == gettextRcmdr("<use y-variable names>")){
      axisLabel <- if (legend && length(z)==0) ""
      else if(length(y) == 1) y
      else paste(paste("(", 1:length(y), ") ", y, sep=""), collapse=", ")
    }
    pch <- if (length(y) == 1) ", pch=1" else ""
    if (legend && length(y) > 1){
      mar <- par("mar")
      top <- 3.5 + length(y)
      command <- paste(".mar <- par(mar=c(", mar[1], ",", mar[2], ",", top, ",", mar[4], "))", sep="")
      logger(command)
      justDoIt(command)
    }
    if(length(z)==0){
      command <- paste("matplot(", .activeDataSet, "$", x, ", ", .activeDataSet, "[, ",
                       paste("c(", paste(paste('"', y, '"', sep=""), collapse=","), ")", sep=""),
                       '], type="',lineType,'", lty=1:',length(y),', ylab="', axisLabel, '"', pch, ")", sep="")
      logger(command)
      justDoIt(command)
    } else {
      doItAndPrint(paste("plotByGroups(x.name='",x,"', y.name='",y,"', z.name='",z,"', lineType='",lineType,"', axisLabel='",axisLabel,"', legend=", legend, ")", sep="")) # ORDNE en legend!!!!!!!!!!!!!!
    }
    if (legend && length(y) > 1){
      n <- length(y)
      cols <- rep(1:6, 1 + n %/% 6)[1:n]
      logger(".xpd <- par(xpd=TRUE)")
      justDoIt(".xpd <- par(xpd=TRUE)")
      usr <- par("usr")
      command <- paste("legend(", usr[1], ", ", usr[4] + 1.2*top*strheight("x"), ", legend=",
                       paste("c(", paste(paste('"', y, '"', sep=""), collapse=","), ")", sep=""),
                       ", col=c(", paste(cols, collapse=","), "), lty=1:",length(y),", pch=c(",
                       paste(paste('"', as.character(1:n), '"', sep=""), collapse=","), "))", sep="")
      logger(command)
      justDoIt(command)
      logger("par(mar=.mar)")
      justDoIt("par(mar=.mar)")
      logger("par(xpd=.xpd)")
      justDoIt("par(xpd=.xpd)")
    }
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="matplot")
  radioButtons(name="lineType",
               buttons=c("b", "l", "p"),
               values=c("b", "l", "p"), initialValue="b",
               labels=gettextRcmdr(c("Lines and points", "Lines", "Points")), title=gettextRcmdr("Plot type"))
  tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw", row=1, column=1, columnspan=2)
  tkgrid(labelRcmdr(axisLabelFrame, text=gettextRcmdr("Label for y-axis"), fg="blue"), sticky="w")
  tkgrid(axisLabelEntry, sticky="w")
  tkgrid(axisLabelScroll, sticky="ew")
  tkgrid(axisLabelFrame, sticky="w", row=2, column=1, columnspan=2)
  tkgrid(labelRcmdr(legendFrame, text=gettextRcmdr("Plot legend")),
         legendCheckBox, sticky="w")
  tkgrid(lineTypeFrame, sticky="w", row=3, column=1, columnspan=1)
  tkgrid(legendFrame, sticky="w", row=3, column=2, columnspan=1)
  tkgrid(buttonsFrame, stick="w", row=4, column=1, columnspan=2)
  dialogSuffix(rows=4, columns=2)
}

#####################################
# Plot mixture points
simplexGUI <- function(){
  initializeDialog(title=gettextRcmdr("Plot points in three component mixture designs"))
  .numeric <- Numeric()
  comboVar <- tclVar()
  comboFrame <- tkframe(top)
  comboVar2 <- tclVar()
  comboFrame2 <- tkframe(top)
  formatFrame <- tkframe(top)
  zoomFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  x1Box <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Left variable (pick one)"))
  x2Box <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Top variable (pick one)"))
  x3Box <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Right variable (pick one)"))
  
  onOK <- function(){ # Actions to perform
    x1 <- getSelection(x1Box)
    x2 <- getSelection(x2Box)
    x3 <- getSelection(x3Box)  	
    mix.format <- as.character(tclvalue(label.formatVariable))
    zoomed <- tclvalue(zoomedVariable)
    if(zoomed == gettextRcmdr("1")){
      zoomed <- "TRUE"
    } else {
      zoomed <- "FALSE"
    }
    n.ticks <- which(as.character(2:10)==tclvalue(comboVar))+1
    if(length(n.ticks)==0) n.ticks <- "6"
    n.grade <- which(as.character(seq(5,25,5))==tclvalue(comboVar2))*5
    if(length(n.grade)==0) n.grade <- "15"
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    formula1 <- justDoIt(paste("formula(~ ",x1," + ",x2," + ",x3, ")", sep=""))
    doItAndPrint(paste("mixture.contour(", .activeDataSet, ", ", paste(formula1[1],formula1[2],sep=" "), ", n.tick=", n.ticks, ", n.grade=", n.grade,", mix.format='",mix.format,"', show.points=TRUE, show.contour=FALSE, zoomed=",zoomed,", pch=21, cex=1.25, col.points='black', fill.points='white')", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="plot")
  tkgrid(getFrame(x1Box), labelRcmdr(variablesFrame, text="    "), getFrame(x2Box), labelRcmdr(variablesFrame, text="    "), getFrame(x3Box), sticky="nw")
  tkgrid(variablesFrame, sticky="nw", row=1, column=1, columnspan=2)
  combo <- ttkcombobox(comboFrame, values=as.character(2:10), textvariable=comboVar, width=3)
  combo2 <- ttkcombobox(comboFrame2, values=as.character(seq(5,25,5)), textvariable=comboVar2, width=3)
  tkgrid(labelRcmdr(comboFrame, text=gettextRcmdr("Plot ticks (default=6):")), combo, sticky="w")
  tkgrid(labelRcmdr(comboFrame2, text=gettextRcmdr("Plot gradings (approximate, default=15):")), combo2, sticky="w")
  tkgrid(comboFrame, sticky="w", column=1, row=2, columnspan=1)
  tkgrid(comboFrame2, sticky="w", column=1, row=3, columnspan=1)
  radioButtonsNMBU(formatFrame,name="label.format", buttons=c("frac", "dec"), values=c("frac", "dec"), initialValue = "frac",
                  labels=gettextRcmdr(c("Fraction", "Decimal")))
  tkgrid(formatFrame, row=4, column=1, columnspan=1, rowspan=1, sticky="w")
  checkBoxes(frame="zoomFrame", boxes=c("zoomed"), initialValues=c("0"), labels=gettextRcmdr(c("Zoom on samples")))
  tkgrid(zoomFrame, row=5, column=1, columnspan=1, rowspan=1, sticky="w")
  tkgrid(buttonsFrame, sticky="w", row=6, column=1, columnspan=1)
  dialogSuffix(rows=6, columns=1)
}


#####################################
# Histogram (discrete frequency)
histogram_discrete <- function () {
	defaults <- list(initial.x = NULL, initial.scale = "frequency", 
			initial.bins = gettextRcmdr ("<auto>"), initial.discrete = "0") 
	dialog.values <- getDialog("histogram_discrete", defaults)
	initializeDialog(title = gettextRcmdr("Histogram"))
	xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	discreteFrame <- tkframe(top)
    comboVar <- tclVar()
	values <- c("normal", "exponential", "gamma", "geometric", "log-normal", "lognormal", "logistic", "negative binomial", "Poisson", "t", "weibull")
	comboFrame <- tkframe(top)
	onOK <- function() {
		x <- getSelection(xBox)
		discrete <- tclvalue(discreteVariable)
		closeDialog()
		if (length(x) == 0) {
			errorCondition(recall = histogram_discrete, message = gettextRcmdr("You must select a variable"))
			return()
		}
		fitDens <- tclvalue(comboVar)
		bins <- tclvalue(binsVariable)
		opts <- options(warn = -1)
		binstext <- if (bins == gettextRcmdr("<auto>")) 
					"\"Sturges\""
				else as.numeric(bins)
		options(opts)
		scale <- tclvalue(scaleVariable)
		putDialog ("histogram_discrete", list (initial.x = x, initial.bins = bins, initial.scale = scale, initial.discrete=discrete))
		if(discrete == gettextRcmdr("0")){
			command <- paste("Hist(", ActiveDataSet(), "$", x, ", scale=\"", 
				scale, "\", breaks=", binstext, ", col=\"darkgray\", xlab='", x, "')", 
				sep = "")
			doItAndPrint(command)
			if(fitDens %in% values){
				if(scale=="density"){
					command <- paste("plotFitDens(", ActiveDataSet(), "$", x, ", '", fitDens, "')", sep="")
				} else {
					eval(parse(text=paste("tmp <- hist(", ActiveDataSet(), "$", x, ", breaks=", binstext, ", plot=FALSE)", sep="")))
					scaling <- max(tmp$counts)/max(tmp$density)
					command <- paste("plotFitDens(", ActiveDataSet(), "$", x, ", '", fitDens, "', scaling=", scaling, ")", sep="")
				}
				doItAndPrint(command)
			}
		} else {
			if(fitDens %in% values)
				errorCondition(recall = histogram_discrete, message = gettextRcmdr("'Density fit' not available for 'Discrete variable'"))
			command <- paste("barplot(.tmp, col=\"darkgray\", space=0, xlab='", x, "',", 
				sep = "")
			doItAndPrint(paste(".tmp <- table(factor(", ActiveDataSet(), "$", x, "))",sep=""))
			if(scale=="percent"){
				doItAndPrint(paste(".tmp <- .tmp/sum(.tmp)*100",sep=""))
				command <- paste(command, "ylab='percent')")
			} else {
				if(scale=="density"){
					doItAndPrint(paste(".tmp <- .tmp/sum(.tmp)",sep=""))
					command <- paste(command, "ylab='density')")
				} else {
					command <- paste(command, "ylab='frequency')")
				}
			}
			
			doItAndPrint(command)
			doItAndPrint("rm(.tmp)")
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "Hist", reset = "Histogram")
	radioButtons(name = "scale", buttons = c("frequency", "percent", 
					"density"), labels = gettextRcmdr(c("Frequency counts", 
							"Percentages", "Densities")), title = gettextRcmdr("Axis Scaling"), 
			initialValue = dialog.values$initial.scale)
	binsFrame <- tkframe(top)
	binsVariable <- tclVar(dialog.values$initial.bins)
	binsField <- ttkentry(binsFrame, width = "8", textvariable = binsVariable)
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(binsFrame, text = gettextRcmdr("Number of bins: ")), 
			binsField, sticky = "w")
	tkgrid(binsFrame, sticky = "w")
	tkgrid(scaleFrame, sticky = "w")
	checkBoxes(frame="discreteFrame", boxes=c("discrete"), initialValues=dialog.values$initial.discrete, labels=gettextRcmdr(c("Discrete variable")))
	tkgrid(discreteFrame, sticky="w")
	combo <- ttkcombobox(comboFrame, values=values, textvariable=comboVar)
    tkgrid(labelRcmdr(comboFrame, text=gettextRcmdr("Density fit:")), combo, sticky="w")
    tkgrid(comboFrame, sticky="w")
	tkgrid(buttonsFrame, sticky = "w")
	tkgrid.configure(binsField, sticky = "e")
	dialogSuffix(rows = 5, columns = 1)
}


#####################################
# Boxplot
boxPlotNMBU <- function () {
	defaults <- list(initial.x = NULL, initial.identifyPoints = 0, initialGroup=NULL, initial.means = 0) 
	dialog.values <- getDialog("boxPlot", defaults)
	initializeDialog(title = gettextRcmdr("Boxplot"))
	xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	identifyVariable <- tclVar(dialog.values$initial.identifyPoints)
	identifyFrame <- tkframe(top)
	identifyCheckBox <- tkcheckbutton(identifyFrame, variable = identifyVariable)
	initial.group <- dialog.values$initial.group
	.groups <- if (is.null(initial.group)) FALSE else initial.group
	onOK <- function() {
		x <- getSelection(xBox)
		identifyPoints <- "1" == tclvalue(identifyVariable)
		means <- tclvalue(meansVariable)
		putDialog ("boxPlot", list(initial.x = x, initial.identifyPoints = identifyPoints, 
						initial.group=if (.groups == FALSE) NULL else .groups, initial.means = "1"==means))
		closeDialog()
		if (length(x) == 0) {
			errorCondition(recall = boxPlotNMBU, message = gettextRcmdr("You must select a variable"))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		var <- paste(.activeDataSet, "$", x, sep = "")
		if (is.null(.groups) || .groups == FALSE) {
			command <- (paste("boxplot(", var, ", ylab=\"", x, 
								"\")", sep = ""))
			logger(command)
			justDoIt(command)
			if(means=="1")
				doItAndPrint(paste("points(1, mean(", var, ", na.rm=TRUE), pch=4)", sep=""))
			if (identifyPoints) {
				RcmdrTkmessageBox(title = "Identify Points", 
						message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
								gettextRcmdr(if (MacOSXP()) 
													"esc key to exit."
												else "right button to exit."), sep = ""), 
						icon = "info", type = "ok")
				doItAndPrint(paste("identify(rep(1, length(", 
								var, ")), ", var, ", rownames(", .activeDataSet, 
								"))", sep = ""))
			}
		}
		else {
			command <- (paste("boxplot(", x, "~", .groups, ", ylab=\"", 
								x, "\", xlab=\"", .groups, "\"", ", data=", .activeDataSet, 
								")", sep = ""))
			logger(command)
			justDoIt(command)
			if(means=="1")
				doItAndPrint(paste("points(1:length(levels(",ActiveDataSet(),"$",.groups,")), tapply(", var, ", ", ActiveDataSet(), "$", .groups, ", mean, na.rm=TRUE), pch=4)", sep=""))
			if (identifyPoints) {
				RcmdrTkmessageBox(title = "Identify Points", 
						message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
								gettextRcmdr(if (MacOSXP()) 
													"esc key to exit."
												else "right button to exit."), sep = ""), 
						icon = "info", type = "ok")
				doItAndPrint(paste("identify(", .activeDataSet, 
								"$", .groups, ", ", var, ", rownames(", .activeDataSet, 
								"))", sep = ""))
			}
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(boxPlot, initialGroup=initial.group, 
			initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
	OKCancelHelp(helpSubject = "boxplot", reset = "boxPlotNMBU")
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(identifyFrame, text = gettextRcmdr("Identify outliers with mouse"), 
					justify = "left"), identifyCheckBox, sticky = "w")
	tkgrid(identifyFrame, stick = "w")
	meansFrame <- tkframe(top)
	checkBoxes(frame="meansFrame", boxes=c("means"), initialValues=dialog.values$initial.means, labels=gettextRcmdr(c("Plot mean value(s)")))
	tkgrid(meansFrame, sticky="w")
	tkgrid(groupsFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 4, columns = 1)
}