# Modified on November 27, 2013 by Jordan Gressel

ologitWords <- function(x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
        cat("\n")
    }

    yname <- rownames(attr(x$terms,"factors"))[1]
    alpha <-.05 
    long <- nrow(x$coefficients)- length(x$beta)
    rev.coef <- x$coefficients[-c(1:long),]
    names <- rownames(rev.coef)
    # Identify coefficients which are factors and print the output
    # Put factor variable names into a list
    counter <- 1
    varnames <- NA
    for(i in 2:length(attr(x$terms,"dataClasses"))){
        if(attr(x$terms,"dataClasses")[i]=="factor"){
            varnames[counter] <- attr(attr(x$terms,"dataClasses")[i],"names")
            counter=counter+1
        }
    }
    # text is the test assumption
    text <- paste("Test Information: This test determines for each independent variable in the model whether that variable is a predicts a significant amount of the variance in ",yname,", all other variables held constant.
                   \n The test assumes proportional odds (the relationship between any pair of groups is the same as any other pair).
                  \n \n",sep="")
    wrapper(text)
    # Interpretations of coefficients that are factors
    if(is.na(varnames[1]) == "FALSE"){
        counter <- 1
        factor.coef <- NA
        for(i in 1:length(varnames)){
            for(j in 1:length(attr(rev.coef,"dimnames")[[1]])){
               if(substr(attr(rev.coef,"dimnames")[[1]][j],1,nchar(varnames[i]))==varnames[i]){
                   factor.coef[counter] <- j
                   counter <- counter+1
                   pval <- rev.coef[j,4]
                   factorname <- substring(names[j],nchar(varnames[i])+4,nchar(names[j])-1)
                   beta <- rev.coef[j,1]
                   if(beta >= 0){
                       sign <- "positive"
                       moreless <- "more" #alt. updown <- "increases"
                   }
                   else if(beta < 0){
                       sign <- "negative"
                       moreless <- "less" #alt. updown <- "decreases"
                   }
                   if(pval >= alpha){
                       text1 <- paste("Test Results: ",factorname," has no statistical relationship with ",yname," (alpha=.05). \n \n")
                       wrapper(text1)
                   }
                   else if(pval < alpha){
                       text1 <- paste("Test Results: Controlling for all other variables in the model, ",varnames[i]," has a statistically significant ",sign," relationship with ",yname,". On average, for a one unit increase in ",factorname,",",yname," is ",round(exp(beta),3),moreless," likely to increase to a higher level. (z=",round(rev.coef[j,3],3),", p=",round(rev.coef[j,4],3),"). \n \n",sep="")
                       wrapper(text1)
                   }
               }
            }
        }
    }
    # Non-factor interpretations
    if(exists("factor.coef")){
        factor.coef <- factor.coef
    }
    if(exists("factor.coef") == "FALSE"){
        factor.coef <- nrow(rev.coef)+1
    }
    for(i in seq(1,nrow(rev.coef))[-factor.coef]){
        pval <- rev.coef[i,4]
        xname <- names[i]
        beta <- rev.coef[i,1]
        if(beta >= 0){
            sign <- "positive"
            updown <- "increases"
        }
        else if(beta < 0){
            sign <- "negative"
            updown <- "decreases"
        }
        if(pval >=  alpha){
            text2 <- paste("Test Results: ",xname," has no statistical relationship with ",yname," (alpha=.05). \n \n",sep="")
            wrapper(text2)
        }
        else if(pval < alpha){
            text2 <- paste("Test Results: Controlling for all other variables in the model, ",xname," has a statistically significant ",sign," relationship with ",yname,". For a one unit increase in ",xname,", the likelihood of increasing to a higher level of ",yname," ",updown," by a factor of ",round(exp(beta),3),". (z=",round(rev.coef[i,3],3),", p=",round(rev.coef[i,4],3),"). \n \n",sep="")
            wrapper(text2)
        }
    }
}

# Modified ordinalRegressionModel.ordinal function from Rcmdr: R Commander
ordinalRegressionModelOrdinal2 <- function(){
	defaults <- list(initial.type="logit")
	dialog.values <- getDialog("ordinalRegressionModelOrdinal2", defaults)
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
		putDialog("ordinalRegressionModelOrdinal2", list(initial.type = tclvalue(modelTypeVariable)))
		formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
		command <- paste("clm(", formula, ', link="', tclvalue(modelTypeVariable),
				'", data=', .activeDataSet, subset, ")", sep="")
		# logger(paste(modelValue, " <- ", command, sep=""))
		### assign(modelValue, justDoIt(with(environment(),command)), envir=SUB)
		doItAndPrint(paste(modelValue,"<-",command))
                ### doItAndPrint(paste("with(SUB,mod  <- summary(", modelValue, "))", sep=""))
                doItAndPrint(paste("mod  <- summary(", modelValue, ")", sep=""))
                # Inserted Code:
                ### doItAndPrint("with(SUB,mod)")
                doItAndPrint("mod")
                ### doItAndPrint("with(SUB,ologitWords(mod))")
                doItAndPrint("ologitWords(mod)")
                # End Insertion
		### activeModel(with(SUB,modelValue))
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
	putDialog("ordinalRegressionModelOrdinal2", NULL)
	ordinalRegressionModelOrdinal2()
}
