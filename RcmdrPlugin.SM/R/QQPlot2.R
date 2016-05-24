
QQPlot2 <- function () {
# this function modified by Martin Maechler
	require("car")
	defaults <- list(initial.x = NULL, initial.identify = 0, initial.dist = "norm", initial.df = "",
			initial.chisqdf = "", initial.fdf1 = "", initial.fdf2 = "", initial.othername = "", 
			initial.otherparam = "")
	dialog.values <- getDialog("QQPlot2", defaults)
	initializeDialog(title = gettextRcmdr("Quantile-Comparison (QQ) Plot"))
	xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	onOK <- function() {
		x <- getSelection(xBox)
		initial.dist <-dist <- tclvalue(distVariable)
		identify <- tclvalue(identifyVariable)
		tdf <- tclvalue(tDfVariable)
		chisqdf <- tclvalue(chisqDfVariable)
		fdf1 <- tclvalue(FDf1Variable)
		fdf2 <- tclvalue(FDf2Variable)
		othername <- tclvalue(otherNameVariable)
		otherparam <- tclvalue(otherParamsVariable)
		putDialog ("QQPlot2", list (initial.x = x, initial.dist = initial.dist,
						initial.identify = identify, initial.df = tdf, initial.chisqdf = chisqdf,
						initial.fdf1 = fdf1, initial.fdf2 = fdf2, initial.othername = othername, 
						initial.otherparam = otherparam))
		closeDialog()
		if (0 == length(x)) {
			errorCondition(recall = QQPlot2, message = gettextRcmdr("You must select a variable."))
			return()
		}
		save <- options(warn = -1)
		on.exit(save)
		retryMe <- function(msg) {
			Message(message = msg, type = "error")
			QQPlot2()
		}
		switch(dist, norm = {
					args <- "dist=\"norm\""
				}, t = {
					df <- tclvalue(tDfVariable)
					df.num <- as.numeric(df)
					if (is.na(df.num) || df.num < 1) {
						retryMe(gettextRcmdr("df for t must be a positive number."))
						return()
					}
					args <- paste("dist=\"t\", df=", df, sep = "")
				}, chisq = {
					df <- tclvalue(chisqDfVariable)
					df.num <- as.numeric(df)
					if (is.na(df.num) || df.num < 1) {
						retryMe(gettextRcmdr("df for chi-square must be a positive number."))
						return()
					}
					args <- paste("dist=\"chisq\", df=", df, sep = "")
				}, f = {
					df1 <- tclvalue(FDf1Variable)
					df2 <- tclvalue(FDf2Variable)
					df.num1 <- as.numeric(df1)
					df.num2 <- as.numeric(df2)
					if (is.na(df.num1) || df.num1 < 1 || is.na(df.num2) || 
							df.num2 < 1) {
						retryMe(gettextRcmdr("numerator and denominator \ndf for F must be positive numbers."))
						return()
					}
					args <- paste("dist=\"f\", df1=", df1, ", df2=", 
							df2, sep = "")
				}, {
					dist <- tclvalue(otherNameVariable)
					params <- tclvalue(otherParamsVariable)
					args <- paste("dist=\"", dist, "\", ", params, sep = "")
				})
		.activeDataSet <- ActiveDataSet()
		if ("1" == tclvalue(identifyVariable)) {
			RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
							gettextRcmdr(if (MacOSXP()) 
												"esc key to exit."
											else "right button to exit."), sep = ""), icon = "info", 
					type = "ok")
			idtext <- paste(", labels=rownames(", .activeDataSet, 
					"), id.method=\"identify\"", sep = "")
		}
		else idtext <- ""
		command <- paste("qqPlot", "(", .activeDataSet, "$", 
				x, ", ", args, idtext, ",ylab= \"",x,"\")", sep = "")
		doItAndPrint(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "qqPlot", reset = "QQPlot2")
	distFrame <- tkframe(top)
	distVariable <- tclVar(dialog.values$initial.dist)
	normalButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "norm")
	tButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "t")
	chisqButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "chisq")
	FButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "f")
	otherButton <- ttkradiobutton(distFrame, variable = distVariable, 
			value = "other")
	tDfFrame <- tkframe(distFrame)
	tDfVariable <- tclVar(dialog.values$initial.df)
	tDfField <- ttkentry(tDfFrame, width = "6", textvariable = tDfVariable)
	chisqDfFrame <- tkframe(distFrame)
	chisqDfVariable <- tclVar(dialog.values$initial.chisqdf)
	chisqDfField <- ttkentry(chisqDfFrame, width = "6", textvariable = chisqDfVariable)
	FDfFrame <- tkframe(distFrame)
	FDf1Variable <- tclVar(dialog.values$initial.fdf1)
	FDf1Field <- ttkentry(FDfFrame, width = "6", textvariable = FDf1Variable)
	FDf2Variable <- tclVar(dialog.values$initial.fdf2)
	FDf2Field <- ttkentry(FDfFrame, width = "6", textvariable = FDf2Variable)
	otherParamsFrame <- tkframe(distFrame)
	otherParamsVariable <- tclVar(dialog.values$initial.otherparam)
	otherParamsField <- ttkentry(otherParamsFrame, width = "30", 
			textvariable = otherParamsVariable)
	otherNameVariable <- tclVar(dialog.values$initial.othername)
	otherNameField <- ttkentry(otherParamsFrame, width = "10", 
			textvariable = otherNameVariable)
	identifyVariable <- tclVar(dialog.values$initial.identify)
	identifyFrame <- tkframe(top)
	identifyCheckBox <- tkcheckbutton(identifyFrame, variable = identifyVariable)
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(identifyFrame, text = gettextRcmdr("Identify observations with mouse")), 
			identifyCheckBox, sticky = "w")
	tkgrid(identifyFrame, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Distribution"), 
					fg = "blue"), columnspan = 6, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Normal")), 
			normalButton, sticky = "w")
	tkgrid(labelRcmdr(tDfFrame, text = gettextRcmdr("df = ")), 
			tDfField, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = "t"), tButton, tDfFrame, 
			sticky = "w")
	tkgrid(labelRcmdr(chisqDfFrame, text = gettextRcmdr("df = ")), 
			chisqDfField, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Chi-square")), 
			chisqButton, chisqDfFrame, sticky = "w")
	tkgrid(labelRcmdr(FDfFrame, text = gettextRcmdr("Numerator df = ")), 
			FDf1Field, labelRcmdr(FDfFrame, text = gettextRcmdr("Denominator df = ")), 
			FDf2Field, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = "F"), FButton, FDfFrame, 
			sticky = "w")
	tkgrid(labelRcmdr(otherParamsFrame, text = gettextRcmdr("Specify: ")), 
			otherNameField, labelRcmdr(otherParamsFrame, text = gettextRcmdr("Parameters: ")), 
			otherParamsField, sticky = "w")
	tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Other")), 
			otherButton, otherParamsFrame, sticky = "w")
	tkgrid(distFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 5, columns = 1)
}


