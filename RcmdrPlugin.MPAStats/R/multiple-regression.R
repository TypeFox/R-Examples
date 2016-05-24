# Modified on March 19, 2014 by Jordan Gressel

# Interpretation function
multipleRegressionWords <- function(x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
        cat("\n")
    }

    Capitalize <- function(y){
        letter <- toupper(substr(y,1,1))
        substr(y,1,1) <- letter
        return(y)
    }

    yname <- rownames(attr(x$terms,"factors"))[1]
    alpha <-.05 
    names <- rownames(x$coefficients)
    # Identify coefficients which are factors and print output
    # Put factor variable names into a list
    counter <- 1
    varnames <- NA
    for(i in 1:length(attr(x$terms,"dataClasses"))){
        if(attr(x$terms,"dataClasses")[i]=="factor"){
            varnames[counter] <- attr(attr(x$terms,"dataClasses")[i],"names")
            counter=counter+1
        }
    }
    # Interpretations of coefficients that are factors
    if(is.na(varnames[1]) == "FALSE"){
    counter <- 1
    factor.coef <- NA
        for(i in 1:length(varnames)){
            for(j in 2:length(attr(x$coefficients,"dimnames")[[1]])){
               if(substr(attr(x$coefficients,"dimnames")[[1]][j],1,nchar(varnames[i]))==varnames[i]){
                   factor.coef[counter] <- j
                   counter <- counter+1
                   pval <- x$coefficients[j,4]
                   factorname <- substring(names[j],nchar(varnames[i])+4,nchar(names[j])-1)
                   beta <- x$coefficients[j,1]
                   if(beta >= 0){
                       sign <- "positive"
                       updown <- "increases"
                   }
                   else if(beta < 0){
                       sign <- "negative"
                       updown <- "decreases"
                   }
                   if(pval >= alpha){
                       factorname <- Capitalize(factorname)
                       text <- paste("Test Results: ",factorname," has no statistical relationship with ",yname," (alpha=.05). \n \n")
                       wrapper(text)
                   }
                   else if(pval < alpha){
                       text <- paste("Test Results: Holding all other variables constant, ",varnames[i]," has a statistically significant ",sign," relationship with ",yname,". ","On average, being ", factorname," ",updown," ",yname," by ",round(abs(beta),3)," units (t=",round(x$coefficients[j,3],3),", p=",round(x$coefficients[j,4],3),"). \n \n",sep="")
                       wrapper(text)
                   } 
               }
            }
        }
    }

    # Non-factor interpretations
    if(exists("factor.coef")){
        factor.coef <- factor.coef - 1
    }
    else if(exists("factor.coef") == "FALSE"){
        factor.coef <- nrow(x$coefficients)+1
    }
    for(i in seq(2,nrow(x$coefficients))[-factor.coef]){
        pval <- x$coefficients[i,4]
        xname <- names[i]
        beta <- x$coefficients[i,1]
            if(beta >= 0){
            sign <- "positive"
            updown <- "increases"
        }
        else if(beta < 0){
            sign <- "negative"
            updown <- "decreases"
        }
        # text is the test assumption
        text <- paste("Test Information: This test determines whether each independent variable predicts variation in ",yname,", all other variables held constant.
                   \n The test assumes that ",yname," is normally distributed with equal variance across all values of any independent variable. Each independent variable must have a linear relationship with ",yname,".
                  \r ****************************************************************
                  \n \n",sep="")
        wrapper(text)
        if(pval >=  alpha){
            xname <- Capitalize(xname)
            text1 <- paste("Test Results: ",xname," has no statistical relationship with ",yname," (alpha=.05). \n \n",sep="") 
            wrapper(text1)
        }
        else if(pval < alpha){
            text1 <- paste("Test Results: Holding all other variables constant, ",xname," has a statistically significant ",sign," relationship with ",yname,". For a one unit increase in ",xname,", ",yname," ",updown," by ",round(abs(beta),3)," units. (t=",round(x$coefficients[i,3],3),", p=",round(x$coefficients[i,4],3),"). \n \n",sep="")
            wrapper(text1)
        } 
    }
    # Interpret Model/R-squared
    text2<-paste("The R-squared value of ",round(x$r.squared,3)," (adjusted R-squared = ",round(x$adj.r.squared,3),") indicates that this formula explains about ",round(x$r.squared,3)*100," percent of the variation in ",yname,", based on results from the current data. \n
      To predict individual values of ",yname," based on the regression line, enter the appropriate values for each independent variable into the formula below.
      Predictions may be under- or over-estimates of the actual value of ",yname,". \n \n",sep="")
    wrapper(text2)

    coeff<-round(x$coef[,1],3)
    xx<-coeff[-1]
    stri<-paste(yname,"=",coeff[1])
    for (i in 1:length(xx)){
    if(xx[i]>0){
    stri<-cat(stri,"+")
        }
     stri<-cat(stri,xx[i],"*",names[i+1])
  if(i==length(xx)){
    stri<-cat(stri,"\n")
        }
    }
}

# Modified linearModel function in Rcmdr: R Commander
linearModel2 <- function(){
	initializeDialog(title=gettextRcmdr("Linear Model"))
	.activeModel <- ActiveModel()
	currentModel <- if (!is.null(.activeModel))
				class(get(.activeModel, envir=.GlobalEnv))[1] == "lm"
			else FALSE
	if (currentModel) {
		currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
		if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
	}
	if (isTRUE(getRcmdr("reset.model"))) {
		currentModel <- FALSE
		putRcmdr("reset.model", FALSE)
	}
	UpdateModelNumber()
	modelName <- tclVar(paste("LinearModel.", getRcmdr("modelNumber"), sep=""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width="20", textvariable=modelName)
	onOK <- function(){
		modelValue <- trim.blanks(tclvalue(modelName))
		closeDialog()
		if (!is.valid.name(modelValue)){
			errorCondition(recall=linearModel2, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
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
			errorCondition(recall=linearModel2, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
			return()
		}
		check.empty <- gsub(" ", "", tclvalue(rhsVariable))
		if ("" == check.empty) {
			errorCondition(recall=linearModel2, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
			return()
		}
		if (is.element(modelValue, listLinearModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
				UpdateModelNumber(-1)
				linearModel2()
				return()
			}
		}
		formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
		
                command <- paste("lm(", formula,
				", data=", ActiveDataSet(), subset, ")", sep="")
		# logger(paste(modelValue, " <- ", command, sep=""))
                doItAndPrint(paste(modelValue, " <- ", command, sep=""))
		###assign(modelValue, justDoIt(with(environment(),command)), envir=SUB)
		# assign(modelValue, justDoIt(command))
                # Inserted "model <-"
		###doItAndPrint(paste("with(SUB,model <- summary(", modelValue, "))", sep=""))               
		doItAndPrint(paste("model <- summary(", modelValue, ")", sep=""))               
                ###doItAndPrint("with(SUB,model)")
                doItAndPrint("model")
                # Inserted Code:
                ###doItAndPrint("with(SUB,multipleRegressionWords(model))")
		doItAndPrint(paste("plot(",ActiveDataSet(),"$",tclvalue(lhsVariable),",",ActiveDataSet(),"$",tclvalue(rhsVariable),")", sep = ""))
		doItAndPrint(paste("abline(",modelValue,")",sep=""))
                doItAndPrint("multipleRegressionWords(model)")
                # End Inserted Code
                activeModel(modelValue)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="linearModel2", model=TRUE, reset="resetLinearModel2")
	tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
	tkgrid(modelFrame, sticky="w")
	modelFormula()
	subsetBox(model=TRUE)
	tkgrid(getFrame(xBox), sticky="w")
	tkgrid(outerOperatorsFrame, sticky="w")
	tkgrid(formulaFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}

resetLinearModel <- function(){
	putRcmdr("reset.model", TRUE)
	linearModel2()
} 
