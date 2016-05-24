# Modified on March 19, 2014 by Jordan Gressel

#Interpretation Function
#dat.factor represents the vector of factor data we are testing
singleProportionTestWords <- function(varname,level,x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    } 

    pval <- x$p.value
    null.value <- x$null.value
    alpha <- 1-level

    up.down <- paste(x$alternative," than ",sep="")
    if(up.down == "two.sided than "){
      if (x$estimate < x$null.value){
        up.down="smaller than "}
      else{
        up.down="larger than "}
    }
    
    # text is the test assumption
    text <- paste("Test Information: This test determines whether the true proportion of ",varname," in the population is significantly ",up.down,null.value,".
                   \n The test assumes that data are randomly and independently sampled. Furthermore,
                    \r -The sample must include at least 30 observations (n >= 30)
                    \r -The size of the population must be at least ten times as large as the sample size (N >= 10n)
                    \r -The sample must include at least 5 successes and 5 failures (n * p >= 5 and n * (1-p) >= 5).
                    \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)
    
    # text1 is the test results
    if(pval >= alpha){
        text1 <- paste("Test Results:The proportion of ",varname ," (", round(x$estimate,2),") in the population is not significantly ",up.down," the hypothesized value ",null.value,". \n \n",sep="")
        wrapper(text1)
    }

    else if(pval < alpha){
        text1 <- paste("Test Results:The proportion of ",varname ," (", round(x$estimate,2),") in the population is significantly ",up.down," the hypothesized value ",null.value," (p=",round(pval,3),"). \n \n",sep="")
        wrapper(text1)
    }
}

#Modified singleProportionTest code from Rcmdr: R Commander
singleProportionTest2 <- function () {
	defaults <- list (initial.x = NULL, initial.alternative = "two.sided", initial.level = ".95", 
			initial.test = "normal" , initial.p = ".5")
	dialog.values <- getDialog ("singleProportionTest2", defaults)
	initializeDialog(title = gettextRcmdr("Single-Sample Proportion Test"))
	xBox <- variableListBox(top, TwoLevelFactors(), title = gettextRcmdr("Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.x,"factor"))
	onOK <- function() {
		x <- getSelection(xBox)
		if (length(x) == 0) {
			errorCondition(recall = singleProportionTest2, message = gettextRcmdr("You must select a variable."))
			return()
		}
	
		alternative <- as.character(tclvalue(alternativeVariable))
		level <- tclvalue(confidenceLevel)
		test <- as.character(tclvalue(testVariable))
		p <- tclvalue(pVariable)
		putDialog ("singleProportionTest2", list (initial.x = x, initial.alternative = alternative, 
						initial.level = level, initial.test = test , initial.p = p))
		closeDialog()
		command <- paste("xtabs(~", x, ", data=", ActiveDataSet(), 
				")")
		logger(paste(".Table <-", command))
		doItAndPrint(paste(".Table","<-",command))
		doItAndPrint(".Table")
    #Creates (and prints) the matrix "ntable" which has the levels and counts of the variable 
		doItAndPrint(paste("ntable<-(rbind(.Table))"))
    #varname takes the 1st column name (or 1st level) 
		varname <- colnames(ntable)[1]
    
                # Added the "model <-" to each of the pastes
		if (test == "normal") 
			doItAndPrint(paste("model <- prop.test(rbind(.Table), alternative='", 
							alternative, "', p=", p, ", conf.level=", level, 
							", correct=FALSE)", sep = ""))
		else if (test == "corrected") 
			doItAndPrint(paste("model <- prop.test(rbind(.Table), alternative='", 
							alternative, "', p=", p, ", conf.level=", level, 
							", correct=TRUE)", sep = ""))
		else doItAndPrint(paste("model <- binom.test(rbind(.Table), alternative='", 
							alternative, "', p=", p, ", conf.level=", level, 
							")", sep = ""))
        doItAndPrint("model")
        #Inserted Code
		    doItAndPrint(paste("singleProportionTestWords(",'"',varname,'",',level,",model)",sep=""))
        #End of Inserted Code
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "prop.test", reset = "singleProportionTest2")
	radioButtons(top, name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Population proportion != p0", 
							"Population proportion < p0", "Population proportion > p0")), 
			title = gettextRcmdr("Alternative Hypothesis"), initialValue = dialog.values$initial.alternative)
	rightFrame <- tkframe(top)
	confidenceFrame <- tkframe(rightFrame)
	confidenceLevel <- tclVar(dialog.values$initial.level)
	confidenceField <- ttkentry(confidenceFrame, width = "6", 
			textvariable = confidenceLevel)
	pFrame <- tkframe(rightFrame)
	pVariable <- tclVar(dialog.values$initial.p)
	pField <- ttkentry(pFrame, width = "6", textvariable = pVariable)
	radioButtons(name = "test", buttons = c("normal", "corrected", 
					"exact"), labels = gettextRcmdr(c("Normal approximation", 
							"Normal approximation with\ncontinuity correction", "Exact binomial")), 
			title = gettextRcmdr("Type of Test"), initialValue = dialog.values$initial.test)
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(pFrame, text = gettextRcmdr("Null hypothesis: p = "), 
					fg = "blue"), pField, sticky = "w")
	tkgrid(pFrame, sticky = "w")
	tkgrid(labelRcmdr(rightFrame, text = ""))
	tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: "), 
					fg = "blue"), confidenceField, sticky = "w")
	tkgrid(confidenceFrame, sticky = "w")
	tkgrid(alternativeFrame, rightFrame, sticky = "nw")
	tkgrid(testFrame, sticky = "w")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	tkgrid.configure(confidenceField, sticky = "e")
	dialogSuffix(rows = 4, columns = 2)
}


# Check for objects of class `sclm`
sclmP <- function() activeModelP() && inherits(get(ActiveModel()), 'sclm')


# Modify the built-in `ordinalRegressionModel()` to work with `clm()` in the ordinal library
ordinalRegressionModel.ordinal <- function(){
	defaults <- list(initial.type="logit")
	dialog.values <- getDialog("ordinalRegressionModel.ordinal", defaults)
	initializeDialog(title=gettextRcmdr("Ordinal Logisitic Regression Model"))
	.activeModel <- ActiveModel()
	.activeDataSet <- ActiveDataSet()
  currentModel <- if (!is.null(.activeModel))
        class(get(.activeModel, envir=.GlobalEnv))[1] == "sclm"
      else FALSE
	if (currentModel) {
		currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
		if (currentFields$data != .activeDataSet) currentModel <- FALSE
	}
	if (isTRUE(getRcmdr("reset.model"))) {
		currentModel <- FALSE
		putRcmdr("reset.model", FALSE)
	}
	UpdateModelNumber()
	modelName <- tclVar(paste("OrdinalLogit.", getRcmdr("modelNumber"), sep=""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width="20", textvariable=modelName)
	radioButtons(name="modelType",
			buttons=c("logit", "probit"), initialValue=dialog.values$initial.type,
			labels=gettextRcmdr(c("Proportional-odds logit", "Ordered probit")),
			title=gettextRcmdr("Type of Model"))
	onOK <- function(){
		modelValue <- trim.blanks(tclvalue(modelName))
		closeDialog()
		if (!is.valid.name(modelValue)){
			errorCondition(recall=proportionalOddsModel, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
			return()
		}
		subset <- tclvalue(subsetVariable)
		if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
			subset <- ""
			putRcmdr("modelWithSubset", FALSE)
		}
		else{
			subset <- paste(", subset=", subset, sep="")
			putRcmdr("modelWithSubset", TRUE)
		}
		check.empty <- gsub(" ", "", tclvalue(lhsVariable))
		if ("" == check.empty) {
			errorCondition(recall=proportionalOddsModel, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
			return()
		}
		check.empty <- gsub(" ", "", tclvalue(rhsVariable))
		if ("" == check.empty) {
			errorCondition(recall=proportionalOddsModel, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
			return()
		}
		if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=get(.activeDataSet, envir=.GlobalEnv)))){
#        if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=eval(parse(text=.activeDataSet), envir=.GlobalEnv)))){
			errorCondition(recall=proportionalOddsModel, message=gettextRcmdr("Response variable must be a factor"))
			return()
		}
		if (is.element(modelValue, listProportionalOddsModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
				UpdateModelNumber(-1)
				proportionalOddsModel()
				return()
			}
		}
		putDialog("ordinalRegressionModel.ordinal", list(initial.type = tclvalue(modelTypeVariable)))
		formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
		command <- paste("clm(", formula, ', link="', tclvalue(modelTypeVariable),
				'", data=', .activeDataSet, subset, ")", sep="")
		doItAndPrint(paste(modelValue, " <- ", command, sep=""))
		doItAndPrint(paste("summary(", modelValue, ")", sep=""))
		activeModel(modelValue)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="clm", model=TRUE, reset = "resetclm")
	tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
	tkgrid(modelFrame, sticky="w")
	modelFormula()
	subsetBox(model=TRUE)
	tkgrid(getFrame(xBox), sticky="w")
	tkgrid(outerOperatorsFrame, sticky="w")
	tkgrid(formulaFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(modelTypeFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}


# Reset the clm model dialog
resetclm <- function() {
	putRcmdr("reset.model", TRUE)
	putDialog("ordinalRegressionModel.ordinal", NULL)
	ordinalRegressionModel.ordinal()
}


#Interpretation Function
twoSampleProportionsTestWords <- function(x,groups,varname,table){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }            
    pval <- round(x$p.value,2)
    alpha <- 1-attr(x$conf.int,"conf.level")
    conf.level <- 100*attr(x$conf.ing,"conf.level")
    chisq <- round(x$statistic,2)

    grp1 <- rownames(table)[1]
    grp2 <- rownames(table)[2]

    prop1 <- x$estimate[1]
    prop2 <- x$estimate[2]
    
    # text is the test assumption
    text <- paste("Test Information: This test determines whether there is a difference in the true proportion of ",varname," between levels of ",groups," in the population.
        \n The test assumes that data are randomly and independently sampled. Furthermore, for each sample,
        \r -The sample must include at least 30 observations (n >= 30)
        \r -The size of the population must be at least ten times as large as the sample size (N >= 10n)
        \r -The sample must include at least 5 successes and 5 failures (n * p >= 5 and n * (1-p) >= 5).
        \r ****************************************************************
      \n \n",sep="")
    wrapper(text)

    if(pval >= alpha){
        text1 <- paste("Test Results: There is no significant difference in the proportion of ",varname, " between groups of ",groups,". (chi-square=",round(chisq,2),", p=",round(pval,3),")." ,sep="")
        wrapper(text1)
    }

    else if(pval < alpha){
        text1 <- paste("Test Results: There is a statistically significant difference in the proportion of ",varname," between groups of ",groups,". The proportion for ",grp1," was ",prop1," and the proportion of ",grp2," was ",prop2,". (chisquare=",round(chisq,2),", p=",round(pval,3),").",sep="")
        wrapper(text1)
    }
}

#Modified twoSampleProportionsTest from Rcmdr: R Commander
twoSampleProportionsTest2 <- function () {
	Library("abind")
	defaults <- list(initial.groups = NULL, initial.response = NULL, initial.alternative = "two.sided", 
			initial.confidenceLevel = ".95", initial.test = "normal", initial.label=NULL)
	dialog.values <- getDialog("twoSampleProportionsTest2", defaults)
	initializeDialog(title = gettextRcmdr("Two-Sample Proportions Test"))
	.twoLevelFactors <- TwoLevelFactors()
	groupsBox <- variableListBox(top, .twoLevelFactors, title = gettextRcmdr("Groups (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.groups, "twoLevelFactor"))
	xBox <- variableListBox(top, .twoLevelFactors, title = gettextRcmdr("Response Variable (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.response, "twoLevelFactor"))
	onOK <- function() {
		groups <- getSelection(groupsBox)
		if (length(groups) == 0) {
			errorCondition(recall = twoSampleProportionsTest2, 
					message = gettextRcmdr("You must select a groups variable."))
			return()
		}
		x <- getSelection(xBox)
		if (length(x) == 0) {
			errorCondition(recall = twoSampleProportionsTest2, 
					message = gettextRcmdr("You must select a response variable."))
			return()
		}
		if (x == groups) {
			errorCondition(recall = twoSampleProportionsTest2, 
					message = gettextRcmdr("Groups and response variables must be different."))
			return()
		}
		alternative <- as.character(tclvalue(alternativeVariable))
		level <- tclvalue(confidenceLevel)
		test <- as.character(tclvalue(testVariable))
		closeDialog()
		putDialog("twoSampleProportionsTest2", list(initial.groups = groups, initial.response = x, 
						initial.test = test, initial.alternative = alternative, initial.confidenceLevel = level,
						initial.label=.groupsLabel))
		command <- paste("xtabs(~", groups, "+", x, ", data=", 
				ActiveDataSet(), ")", sep = "")
		logger(paste(".Table <-", command))
		doItAndPrint(paste(".Table","<-",command))
		doItAndPrint("rowPercents(.Table)")
		if (test == "normal") 

                #Inserted "model1 <- " twice.
			doItAndPrint(paste("model1 <- prop.test(.Table, alternative='", 
							alternative, "', conf.level=", level, ", correct=FALSE)", 
							sep = ""))
		else doItAndPrint(paste("model1 <- prop.test(.Table, alternative='", 
							alternative, "', conf.level=", level, ", correct=TRUE)", 
							sep = ""))
        #Inserted Code:
        doItAndPrint("model1")
		#Creates (and prints) the matrix "ntable" which has the levels and counts of the variable 
		doItAndPrint(paste("ntable<-(rbind(.Table))"))
		#varname takes the 1st column name (or 1st level) 
		varname <- colnames(ntable)[1]
        doItAndPrint(paste("twoSampleProportionsTestWords(model1,",'"',groups,'"',", varname=",'"',varname,'"',",.Table)",sep=""))
        #End Inserted Code
		logger("remove(.Table)")
		#remove(.Table, envir = .GlobalEnv)  
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "prop.test", reset = "twoSampleProportionsTest2")
	radioButtons(name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
							"Difference > 0")), initialValue = dialog.values$initial.alternative, 
			title = gettextRcmdr("Alternative Hypothesis"))
	rightFrame <- tkframe(top)
	confidenceFrame <- tkframe(rightFrame)
	confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
	confidenceField <- ttkentry(confidenceFrame, width = "6", 
			textvariable = confidenceLevel)
	radioButtons(name = "test", buttons = c("normal", "corrected"), 
			labels = gettextRcmdr(c("Normal approximation", "Normal approximation with\ncontinuity correction")), 
			initialValue = dialog.values$initial.test, 
			title = gettextRcmdr("Type of Test"))
	tkgrid(getFrame(groupsBox), getFrame(xBox), sticky = "nw")
	groupsLabel(columnspan = 2, initialText=dialog.values$initial.label)
	tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: "), 
					fg = "blue"), confidenceField, sticky = "w")
	tkgrid(confidenceFrame, sticky = "w")
	tkgrid(alternativeFrame, rightFrame, sticky = "nw")
	tkgrid(testFrame, sticky = "w")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	tkgrid.configure(confidenceField, sticky = "e")
	dialogSuffix(rows = 5, columns = 2)
}

