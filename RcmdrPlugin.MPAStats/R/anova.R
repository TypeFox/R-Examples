#Modified on June 13, 2013 by Christa Schank

#Interpretation Function
wordsAnova <- function(x,group,response){
    #Word-wrap function
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }

    fstat <- x[[1]][["F value"]][1]
    pval <- x[[1]][["Pr(>F)"]][1]
    alpha <- .05
    
    text <- paste("Test Information: This test determines whether any level of ",group," affects ",response,"; 
      that is, whether the true mean ",response," levels of ",group," are all equivalent to each other, 
      or if the mean ",response," of at least one level of ",group," is significantly different from the others.
      \n The test assumes one categorical independent variable. Observations are sampled independently. Residuals in each category are distributed as a Normal,
      with equal variance across all levels of ",group,".
      \r **************************************************
    \n \n",sep="")
    wrapper(text)
    
    if(pval >= alpha){
       text1 <- paste("Test Results: There is no significant difference in the mean ",response," between levels of ",group,". (F=",round(fstat,3),", p=",round(pval,3),").",sep="")
       wrapper(text1)
    }
    else if(pval < alpha){
        text1 <- paste("Test Results: At least one mean ",response," among the levels of ",group," differs from the rest. (F=",round(fstat,3),",p=",round(pval,3),").",sep="")
        wrapper(text1)
    }
}

#Modified oneWayAnova function code from Rcmdr R Commander.
oneWayAnova2 <- function () {
	Library("multcomp")
	Library("abind")
	defaults <- list(initial.group = NULL, initial.response = NULL, initial.pairwise = 0)
	dialog.values <- getDialog("oneWayAnova2", defaults)
	initializeDialog(title = gettextRcmdr("One-Way Analysis of Variance"))
	UpdateModelNumber()
	modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), 
					sep = ""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
	groupBox <- variableListBox(top, Factors(), title = gettextRcmdr("Groups (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.group, "factor"))
	responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	optionsFrame <- tkframe(top)
	pairwiseVariable <- tclVar(dialog.values$initial.pairwise)
	pairwiseCheckBox <- tkcheckbutton(optionsFrame, variable = pairwiseVariable)
	onOK <- function() {
		modelValue <- trim.blanks(tclvalue(modelName))
		if (!is.valid.name(modelValue)) {
			UpdateModelNumber(-1)
			errorCondition(recall = oneWayAnova2, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
							modelValue))
			return()
		}
		if (is.element(modelValue, listAOVModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type = gettextRcmdr("Model")))) {
				UpdateModelNumber(-1)
				tkdestroy(top)
				oneWayAnova2()
				return()
			}
		}
		group <- getSelection(groupBox)
		response <- getSelection(responseBox)
		closeDialog()
		if (length(group) == 0) {
			errorCondition(recall = oneWayAnova2, message = gettextRcmdr("You must select a groups factor."))
			return()
		}
		if (length(response) == 0) {
			errorCondition(recall = oneWayAnova2, message = gettextRcmdr("You must select a response variable."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		command <- paste(modelValue, " <- aov(", response, " ~ ", 
				group, ", data=", .activeDataSet, ")", sep = "")
		justDoIt(command)
		logger(command)
                # added the  "anova1 <-"
		doItAndPrint(paste("anova1 <-summary(", modelValue, ")", sep = ""))
		# Inserted Code:
                doItAndPrint("anova1")
                doItAndPrint(paste("wordsAnova(anova1,",'"',group,'"',",",'"',response,'"',")",sep=""))
                # End Inserted Code
		doItAndPrint(paste("numSummary(", .activeDataSet, "$", 
						response, " , groups=", .activeDataSet, "$", group, 
						", statistics=c(\"mean\", \"sd\"))", sep = ""))
		activeModel(modelValue)
		pairwise <- tclvalue(pairwiseVariable)
		putDialog ("oneWayAnova2", list (initial.group = group, initial.response = response, initial.pairwise = pairwise))
		if (pairwise == 1) {
			if (eval(parse(text = paste("length(levels(", .activeDataSet, 
									"$", group, ")) < 3")))) 
				Message(message = gettextRcmdr("Factor has fewer than 3 levels; pairwise comparisons omitted."), 
						type = "warning")
			else {
				command <- paste(".Pairs <- glht(", modelValue, 
						", linfct = mcp(", group, " = \"Tukey\"))", 
						sep = "")
				justDoIt(command)
				logger(command)
				doItAndPrint("summary(.Pairs) # pairwise tests")
				doItAndPrint("confint(.Pairs) # confidence intervals")
				doItAndPrint("cld(.Pairs) # compact letter display")
				justDoIt("old.oma <- par(oma=c(0,5,0,0))")
				logger("old.oma <- par(oma=c(0,5,0,0))")
				justDoIt("plot(confint(.Pairs))")
				logger("plot(confint(.Pairs))")
				justDoIt("par(old.oma)")
				logger("par(old.oma)")
				logger("remove(.Pairs)")
			#	remove(.Pairs, envir = .GlobalEnv)
			}
		}
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "anova", model = TRUE, reset = "oneWayAnova2")
	tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model: ")), 
			model, sticky = "w")
	tkgrid(modelFrame, sticky = "w", columnspan = 2)
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Pairwise comparisons of means")), 
			pairwiseCheckBox, sticky = "w")
	tkgrid(optionsFrame, sticky = "w", columnspan = 2)
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 4, columns = 2)
}

