# Modified on Feb 13, 2014 by Jordan Gressel

# Interpretation Function
logitWords <- function(x){
    yname <- rownames(attr(x$terms,"factors"))[1]
    alpha <-.05 
    names <- rownames(x$coefficients)
    # Identify coefficient which are factors and print output
    # Put factor variable names into a list
    counter <- 1
    varnames <- NA
    for(i in 1:length(attr(x$terms,"dataClasses"))){
        if(attr(x$terms,"dataClasses")[i]=="factor"){
            varnames[counter] <- attr(attr(x$terms,"dataClasses")[i],"names")
            counter=counter+1
        }
    }
    # Test Assumption
    cat("Test Information: This test determines whether each independent variable predicts variation in ",yname,", all other variables held constant.
                   \n The test assumes that observations are independent, and that the model contains the appropriate set of variables.
                    \r ****************************************************************
                    \n",sep="")
    
    
    
    # Interpretation of coefficients if they are factors
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
			     moreless <- "more" #alt. updown <- "increases"
                   }
                   else if(beta < 0){
                       sign <- "negative"
			     moreless <- "less" #alt. updown <- "decreases"
                   }
                   if(pval >= alpha){
                       cat("Test Results: ",factorname," has no statistical relationship with ",yname," (alpha=.05). \n \n")
                   }
                   else if(pval < alpha){
                       cat("Test Results: Controlling for all other variables in the model ",varnames[i]," has a statistically significant ",sign," relationship with ",yname,". "," For every one unit increase in ",factorname,",",yname," is about ",1/round(exp(beta),3)," times ",moreless," likely to occur. (z=",round(x$coefficients[j,3],3),", p=",round(x$coefficients[j,4],3),"). \n \n",sep="")
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
        if(pval >=  alpha){
            cat("Test Results: ",xname," has no statistical relationship with ",yname," (alpha=.05). \n \n",sep="")
        }
        else if(pval < alpha){
            cat("Test Results: Controlling for all other variables in the model, ",xname," has a statistically significant ",sign," relationship with ",yname,". For a one unit increase in ",xname,", the likelihood of ",yname," ",updown," by a factor of ",round(exp(beta),3)," units. (z=",round(x$coefficients[i,3],3),", p=",round(x$coefficients[i,4],3),"). \n \n",sep="")
        }
    }
}


# Modified generalizedLinearModel function from Rcmdr: R Commander
generalizedLinearModel2  <- function(){
	families <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
			"quasibinomial", "quasipoisson")
	links <- c("identity", "inverse", "log", "logit", "probit",
			"cloglog", "sqrt", "1/mu^2")
	availableLinks <- matrix(c(
					TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
					FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, FALSE,
					TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE,  FALSE,
					TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
					TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, TRUE,
					FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, FALSE,
					TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE,  FALSE),
			7, 8, byrow=TRUE)
	rownames(availableLinks) <- families
	colnames(availableLinks) <- links
	canonicalLinks <- c("identity", "logit", "log", "inverse", "1/mu^2", "logit", "log")
	names(canonicalLinks) <- families
	initializeDialog(title=gettextRcmdr("Generalized Linear Model"))
	.activeModel <- ActiveModel()
	currentModel <- if (!is.null(.activeModel))
				class(get(.activeModel, envir=.GlobalEnv))[1] == "glm"
			else FALSE
	if (currentModel) {
		currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), glm=TRUE)
		if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
	}
	if (isTRUE(getRcmdr("reset.model"))) {
		currentModel <- FALSE
		putRcmdr("reset.model", FALSE)
	}
	modelFormula()
	UpdateModelNumber()
	modelName <- tclVar(paste("GLM.", getRcmdr("modelNumber"), sep=""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width="20", textvariable=modelName)
	linkFamilyFrame <- tkframe(top)
	familyFrame <- tkframe(linkFamilyFrame)
	familyBox <- tklistbox(familyFrame, height="4", exportselection="FALSE",
			selectmode="single", background="white")
	familyScroll <- ttkscrollbar(familyFrame,
			command=function(...) tkyview(familyBox, ...))
	tkconfigure(familyBox, yscrollcommand=function(...) tkset(familyScroll, ...))
	for (fam in families) tkinsert(familyBox, "end", fam)
	linkFrame <- tkframe(linkFamilyFrame)
	linkBox <- tklistbox(linkFrame, height="4", exportselection="FALSE",
			selectmode="single", background="white")
	subsetBox(model=TRUE)
	onFamilySelect <- function(){
		family <- families[as.numeric(tkcurselection(familyBox)) + 1]
		availLinks <- links[availableLinks[family,]]
		tkdelete(linkBox, "0", "end")
		for (lnk in availLinks) tkinsert(linkBox, "end", lnk)
		canLink <- canonicalLinks[family]
		tkconfigure(linkBox, height=length(availLinks))
		tkselection.set(linkBox, which(canLink == availLinks) - 1)
	}
	onOK <- function(){
		check.empty <- gsub(" ", "", tclvalue(lhsVariable))
		if ("" == check.empty) {
			errorCondition(recall=generalizedLinearModel2, model=TRUE, message=gettextRcmdr("Left-hand side of model empty."))
			return()
		}
		check.empty <- gsub(" ", "", tclvalue(rhsVariable))
		if ("" == check.empty) {
			errorCondition(recall=generalizedLinearModel2, model=TRUE, message=gettextRcmdr("Right-hand side of model empty."))
			return()
		}
		modelValue <- trim.blanks(tclvalue(modelName))
		if (!is.valid.name(modelValue)){
			errorCondition(recall=generalizedLinearModel2, model=TRUE, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
			return()
		}
		if (is.element(modelValue, listGeneralizedLinearModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
				UpdateModelNumber(-1)
				closeDialog()
				generalizedLinearModel2()
				return()
			}
		}
		formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
		family <- families[as.numeric(tkcurselection(familyBox)) + 1]
		availLinks <- links[availableLinks[family,]]
		link <- availLinks[as.numeric(tkcurselection(linkBox)) + 1]
		subset <- tclvalue(subsetVariable)
		closeDialog()
		if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
			subset <- ""
			putRcmdr("modelWithSubset", FALSE)
		}
		else{
			subset <- paste(", subset=", subset, sep="")
			putRcmdr("modelWithSubset", TRUE)
		}
		command <- paste("glm(", formula, ", family=", family, "(", link,
				"), data=", ActiveDataSet(), subset, ")", sep="")
		# logger(paste(modelValue, " <- ", command, sep=""))
		### assign(modelValue, justDoIt(with(environment(),command)), envir=SUB)
		doItAndPrint(paste(modelValue,"<-", command))
                # Modified code: (inserted the logit <-) 
		### doItAndPrint(paste("with(SUB,logit <- summary(", modelValue, "))", sep=""))
		doItAndPrint(paste("logit1 <- summary(", modelValue, ")", sep=""))
                # Inserted code:
                ### doItAndPrint("with(SUB,logit)")
                doItAndPrint("logit1")
                ### doItAndPrint("with(SUB,logitWords(logit))")
                doItAndPrint("logitWords(logit1)")
                # End insertion
		### activeModel(with(SUB,modelValue))
		activeModel(modelValue)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="generalizedLinearModel2", model=TRUE, reset="resetGLM")
	helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
	tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
	tkgrid(modelFrame, sticky="w")
	tkgrid(getFrame(xBox), sticky="w")
	tkgrid(outerOperatorsFrame, sticky="w")
	tkgrid(formulaFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(labelRcmdr(linkFamilyFrame, text=gettextRcmdr("Family (double-click to select)"), fg="blue"),
			labelRcmdr(linkFamilyFrame, text="   "), labelRcmdr(linkFamilyFrame, text=gettextRcmdr("Link function"), fg="blue"), sticky="w")
	tkgrid(familyBox, familyScroll, sticky="nw")
	tkgrid(linkBox, sticky="nw")
	tkgrid(familyFrame, labelRcmdr(linkFamilyFrame, text="   "), linkFrame, sticky="nw")
	tkgrid(linkFamilyFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	tkgrid.configure(familyScroll, sticky="ns")
	fam <- if (currentModel) which(currentFields$family == families) - 1
			else 1
	tkselection.set(familyBox, fam)
	availLinks <- links[availableLinks[fam + 1,]]
	for (lnk in availLinks) tkinsert(linkBox, "end", lnk)
	tkconfigure(linkBox, height=length(availLinks))
	lnk <- if (currentModel) which(currentFields$link == availLinks) - 1
			else 0
	tkselection.set(linkBox, lnk)
	tkbind(familyBox, "<Double-ButtonPress-1>", onFamilySelect)
	dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
}

resetGLM <- function(){
	putRcmdr("reset.model", TRUE)
	generalizedLinearModel2()
}
