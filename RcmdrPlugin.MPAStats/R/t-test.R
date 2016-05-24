# Modified on April 23, 2014 by Jordan Gressel

# Interpretation function
singleTTestWords <- function(x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }

    varcut <- nchar(ActiveDataSet())+2
    
    compare=paste(x$alternative," than ",sep="")
    if(compare == "two.sided than "){
        if (x$estimate < x$null.value){
          compare="smaller than "}
        else{
          compare="larger than "}
    }

    m <- x$estimate
    pval=x$p.value
    alpha=1-attr(x$conf.int,"conf.level")


    if(pval < alpha){
        .not=""
    }
    if(pval > alpha){
        .not="not "
    }

    # text is the test assumption
    text <- paste("Test Information: This test determines whether the true population mean is significantly ",compare," {greater than, less than, or different from} ",x$null.value,".
                   \n The test assumes that data are collected randomly from a normally-distributed population.
                    \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)
    
    # text1 is the test results
    text1 <- paste("Test Results: The ",tolower(substring(x$data.name,varcut)) ," mean ", round(m,4) ," is ",.not,"significantly ",compare," the hypothesized value ", x$null.value,". (t=",round(x$statistic,3)," p=",round(x$p.value,3),"). \n",sep="") 
    wrapper(text1)
}

# Modified singleSampleTTest function from Rcmdr: R Commander
singleSampleTTest2 <- function () {
	defaults <- list (initial.x = NULL, initial.alternative = "two.sided", initial.level = ".95", initial.plots="No",
			initial.mu = "0.0")
	dialog.values <- getDialog ("singleSampleTTest2", defaults)  
	initializeDialog(title = gettextRcmdr("Single-Sample t-Test"))
	xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.x, "numeric"))
	onOK <- function() {
		x <- getSelection(xBox)
		if (length(x) == 0) {
			errorCondition(recall = singleSampleTTest2, message = gettextRcmdr("You must select a variable."))
			return()
		}
		alternative <- as.character(tclvalue(alternativeVariable))
		level <- tclvalue(confidenceLevel)
    	mu <- tclvalue(muVariable)
		plots <- as.character(tclvalue(plotsVariable))
		putDialog ("singleSampleTTest2", list (initial.x = x, initial.alternative = alternative, 
						initial.level = level, initial.plots=plots, initial.mu = mu))
		closeDialog()
                # Inserted "t.test <-"
		doItAndPrint(paste("t.test1 <- t.test(", ActiveDataSet(), "$", x, 
						", alternative='", alternative, "', mu=", mu, ", conf.level=", 
						level, ")", sep = ""))
		doItAndPrint("t.test1")
    
    #added plot option
		if(plots == "Yes"){
		  doItAndPrint(paste("graphtest <- ", ActiveDataSet(), "$", x, sep = ""))  
		  doItAndPrint(paste("hist(graphtest, xlab='",ActiveDataSet(), "$", x,"', main='Histogram of ",ActiveDataSet(), "$", x,"')", sep = ""))
		}
	
        doItAndPrint("singleTTestWords(t.test1)")
                # End Inserted Code
		tkdestroy(top)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "t.test", reset = "singleSampleTTest2")
	radioButtons(top, name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Population mean != mu0", "Population mean < mu0", 
							"Population mean > mu0")), title = gettextRcmdr("Alternative Hypothesis"),
			initialValue = dialog.values$initial.alternative)
	rightFrame <- tkframe(top)
	radioButtons(top, name = "plots", buttons = c("Yes", 
	                                                       "No"), values = c("Yes", "No"), 
	             labels = gettextRcmdr(c("Yes",
	                                     "No")), title = gettextRcmdr("Plot?"),
	             initialValue = dialog.values$initial.plot)
	confidenceFrame <- tkframe(rightFrame)
	confidenceLevel <- tclVar(dialog.values$initial.level)
	confidenceField <- ttkentry(confidenceFrame, width = "6", 
			textvariable = confidenceLevel)
	muFrame <- tkframe(rightFrame)
	muVariable <- tclVar(dialog.values$initial.mu)
	muField <- ttkentry(muFrame, width = "8", textvariable = muVariable)
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
	tkgrid(labelRcmdr(muFrame, text = gettextRcmdr("Null hypothesis: mu = ")), 
			muField, sticky = "w")
	tkgrid(muFrame, sticky = "w")
	tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: ")), 
			confidenceField, sticky = "w")
	tkgrid(confidenceFrame, sticky = "w")
	tkgrid(alternativeFrame, rightFrame, sticky = "nw")
	tkgrid(plotsFrame, rightFrame, sticky = "w")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	tkgrid.configure(confidenceField, sticky = "e")
	dialogSuffix(rows = 4, columns = 2)
}


# Interpretation Function
pairedTTestWords=function(x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }


    prelim.names <- strsplit(x$data.name," ")[[1]][c(1,3)]
    one.pos <- which(strsplit(prelim.names,"")[[1]]=="$") + 1
    names <- substring(prelim.names[1],one.pos)

    two.pos <- which(strsplit(prelim.names,"")[[2]]=="$") + 1
    names[2] <- substring(prelim.names[2],two.pos)
        grp1 <- names[1]
        grp2 <- names[2]
    t.value <- x$statistic
    pval <- x$p.value
    conf.level <- 100*attr(x$conf.int,"conf.level")
    alpha <- 1-attr(x$conf.int,"conf.level")

    l.conf <- x$conf.int[1]
    u.conf <- x$conf.int[2]
    
    up.down <- paste(x$alternative," than ",sep="")
    if(x$alternative == "two.sided"){
        one.two <- "two"
        up.down <- "different from "
    }
    else if(x$alternative != "two.sided"){
        one.two <- "one"
    }
    
    # text is the test assumption
    text <- paste("Test Information: This test determines whether there is a difference between ",grp1," and ",grp2," in the population, where ",grp1," and ",grp2," are organized into connected pairs; that is, whether the true mean population difference ",grp1," - ",grp2," is ",up.down," 0.
                   \n The test assumes that data in each group are collected randomly and that the population of differences is normally-distributed.
                    \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)

    # text1 is the test results
    if(pval >= alpha){
        text1 <- paste("Test Results: The mean difference of ",grp1," - ",grp2," is not significantly ",up.down,"0. (t=",round(t.value,3),", p=",round(pval,3)," for a ",one.two,"-tailed test).",sep="")
        wrapper(text1)
    }
    else if(pval < alpha){
        text1 <- paste("Test Results: It appears that on average,",grp1," is different from ",grp2,". The true mean difference for ",grp1," -  ",grp2, " is likely between ",round(l.conf,2)," and ",round(u.conf,2),". (",conf.level,"% confidence, t=",round(t.value,3),", p=",round(pval,3)," for a ",one.two,"-tailed test). \n \n",sep="")
    wrapper(text1)
    }

}    

# Modifed pairedTTest from Rcmdr: R Commander
pairedTTest2 <- function () {
	defaults <- list(initial.x = NULL, initial.y = NULL, initial.alternative = "two.sided", initial.plots="No",
			initial.confidenceLevel = ".95")
	dialog.values <- getDialog("pairedTTest2", defaults)
	initializeDialog(title = gettextRcmdr("Paired t-Test"))
	.numeric <- Numeric()
	xBox <- variableListBox(top, .numeric, title = gettextRcmdr("First variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.x, "numeric"))
	yBox <- variableListBox(top, .numeric, title = gettextRcmdr("Second variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.y, "numeric"))
	onOK <- function() {
		x <- getSelection(xBox)
		y <- getSelection(yBox)
		if (length(x) == 0 | length(y) == 0) {
			errorCondition(recall = pairedTTest2, message = gettextRcmdr("You must select two variables."))
			return()
		}
		if (x == y) {
			errorCondition(recall = pairedTTest2, message = gettextRcmdr("Variables must be different."))
			return()
		}
		alternative <- as.character(tclvalue(alternativeVariable))
		level <- tclvalue(confidenceLevel)
		plots <- as.character(tclvalue(plotsVariable))
		putDialog ("pairedTTest2", list (initial.x = x, initial.y = y, initial.alternative = alternative, initial.plots=plots,
						initial.confidenceLevel = level))
		closeDialog()
		.activeDataSet <- ActiveDataSet()
                # Added "t.test2 <-"
                doItAndPrint(paste("t.test2 <-t.test(", .activeDataSet, "$", x, 
						", ", .activeDataSet, "$", y, ", alternative='", 
						alternative, "', conf.level=", level, ", paired=TRUE)", 
						sep = ""))
                # Inserted Code:
            doItAndPrint("t.test2")
    
		  if(plots == "Yes"){
		  doItAndPrint(paste("graphtest2 <- cbind(",.activeDataSet, "$", x," , ", .activeDataSet, "$", y,")", sep = ""))
		  doItAndPrint(paste("boxplot(graphtest2, names=c('",.activeDataSet, "$", x,"' , '", .activeDataSet, "$", y,"'))", sep = ""))
		  }
                doItAndPrint("pairedTTestWords(t.test2)")
                # End Inserted Code
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "t.test", reset = "pairedTTest2")
	radioButtons(top, name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
							"Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"), 
			initialValue = dialog.values$initial.alternative)
	radioButtons(top, name = "plots", buttons = c("Yes", 
	                                                       "No"), values = c("Yes", "No"), 
	             labels = gettextRcmdr(c("Yes",
	                                     "No")), title = gettextRcmdr("Plot?"),
	             initialValue = dialog.values$initial.plot)
	confidenceFrame <- tkframe(top)
	confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
	confidenceField <- ttkentry(confidenceFrame, width = "6", 
			textvariable = confidenceLevel)
	tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw")
	tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level"), 
					fg = "blue"))
	tkgrid(confidenceField, sticky = "w")
	tkgrid(alternativeFrame, confidenceFrame, sticky = "nw")
	tkgrid(plotsFrame, confidenceFrame, sticky = "w") 
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 3, columns = 2)
}



# Interpretation function
independentSamplesTTestWords <- function(x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }

    varname=strsplit(x$data.name," ")[[1]][1]
    if(x$estimate[1] > x$estimate[2]){
        grp1=strsplit(attr(x$estimate,"names")[1]," ")[[1]][4]
        grp2=strsplit(attr(x$estimate,"names")[2]," ")[[1]][4]
        meangrp1=x$estimate[1]
        meangrp2=x$estimate[2]
    }
    if(x$estimate[1] <= x$estimate[2]){
        grp1=strsplit(attr(x$estimate,"names")[2]," ")[[1]][4]
        grp2=strsplit(attr(x$estimate,"names")[1]," ")[[1]][4]
        meangrp1=x$estimate[2]
        meangrp2=x$estimate[1]
    }
    pval=x$p.value
    alpha <- 1 - attr(x$conf.int,"conf.level")
    t.value=x$statistic
    difference <-abs(x$estimate[1] - x$estimate[2])  
    one.two="one"
    if(x$alternative == "two.sided"){
        one.two="two"
    }
    
    # text is the test assumption
    text <- paste("Test Information: This test determines whether there is a difference in the true population means of ",grp1," and ",grp2," , respectively.
                   \n The test assumes that data in each group are collected randomly from a normally-distributed population.
                    \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)
    
    # text1 is the test results
    if(pval < alpha){
        text1 <- paste("Test Results: On average, the mean ",varname," for ",grp1, " is about ",
            round(difference,2)," more than the mean ",varname," for ",grp2,
            ". The average ",varname," among ",grp1," was about ",round(meangrp1,2),
            ", and the average ",varname," among ",grp2," was about ",
            round(meangrp2,2),". (t=",round(t.value,3),", p=",round(pval,3),", for a ",one.two,
            "-tailed test). \n \n",sep="")
        wrapper(text1)
    }

    else if(pval >= alpha){
        text1 <- paste("Test Results: There was no significant difference in the means of ",grp1," and ",grp2,". (alpha=",alpha,", ",one.two,"-tailed test).",sep="")
        wrapper(text1)
    }
}


# Modified from independentSamplesTTest from Rcmdr: R Commander
independentSamplesTTest2 <- function () {
	defaults <- list(initial.group = NULL, initial.response = NULL, initial.alternative = "two.sided", 
			initial.confidenceLevel = ".95", initial.variances = "FALSE", initial.plots="No", initial.label=NULL)
	dialog.values <- getDialog("independentSamplesTTest2", defaults)
	initializeDialog(title = gettextRcmdr("Independent Samples t-Test"))
	variablesFrame <- tkframe(top)
	groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), 
			title = gettextRcmdr("Groups (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.group, "twoLevelFactor"))
	responseBox <- variableListBox(variablesFrame, Numeric(), 
			title = gettextRcmdr("Response Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	onOK <- function() {
		group <- getSelection(groupBox)
		if (length(group) == 0) {
			errorCondition(recall = independentSamplesTTest2, 
					message = gettextRcmdr("You must select a groups variable."))
			return()
		}
		response <- getSelection(responseBox)
		if (length(response) == 0) {
			errorCondition(recall = independentSamplesTTest2, 
					message = gettextRcmdr("You must select a response variable."))
			return()
		}
		alternative <- as.character(tclvalue(alternativeVariable))
		level <- tclvalue(confidenceLevel)
		variances <- as.character(tclvalue(variancesVariable))
		plots <- as.character(tclvalue(plotsVariable))
		putDialog ("independentSamplesTTest2", list (initial.group = group, initial.response = response, initial.alternative = alternative, 
						initial.confidenceLevel = level, initial.variances = variances, initial.plots=plots, initial.label=.groupsLabel))        
		closeDialog()
                # Added  "t.test3 <-" 
		doItAndPrint(paste("t.test3 <- t.test(", response, "~", group, ", alternative='", 
						alternative, "', conf.level=", level, ", var.equal=", 
						variances, ", data=", ActiveDataSet(), ")", sep = ""))
                # Inserted Code:
                doItAndPrint("t.test3")
    
		doItAndPrint("independentSamplesTTestWords(t.test3)")
    
		if(plots == "Yes"){
      doItAndPrint(paste("boxplot(", response, "~", group, ", data= ",ActiveDataSet(), ", main='Independent Two Group T-test Boxplot')", sep = ""))
    }
	
                # End Insertion
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "t.test", reset = "independentSamplesTTest2")
	optionsFrame <- tkframe(top)
	radioButtons(optionsFrame, name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
							"Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"),
			initialValue = dialog.values$initial.alternative)
	#####ADDED PLOT BUTTONS
	radioButtons(optionsFrame, name = "plots", buttons = c("Yes", 
	                                                             "No"), values = c("Yes", "No"), 
	             labels = gettextRcmdr(c("Yes",
                                      "No")), title = gettextRcmdr("Plot?"),
	             initialValue = dialog.values$initial.plot)
	confidenceFrame <- tkframe(optionsFrame)
	confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
	confidenceField <- ttkentry(confidenceFrame, width = "6", 
			textvariable = confidenceLevel)
	radioButtons(optionsFrame, name = "variances", buttons = c("yes", 
					"no"), values = c("TRUE", "FALSE"),  
			labels = gettextRcmdr(c("Yes", "No")), title = gettextRcmdr("Assume equal variances?"),
			initialValue = dialog.values$initial.variances)
	tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text = "    "), 
			getFrame(responseBox), sticky = "nw")
	tkgrid(variablesFrame, sticky = "nw")
	tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level"), 
					fg = "blue"), sticky = "w")
	tkgrid(confidenceField, sticky = "w")
  groupsLabel(groupsBox = groupBox, initialText=dialog.values$initial.label)
	tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text = "    "), 
			confidenceFrame, labelRcmdr(optionsFrame, text = "    "), 
			plotsFrame, labelRcmdr(optionsFrame, text = "    "), 
			variancesFrame, sticky = "nw")
	tkgrid(optionsFrame, sticky = "nw")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 4, columns = 1)
} 



### Added function for Distribution Diagnostics
DistWords <- function(x){
  wrapper <- function(text){
    text2 <- strwrap(text)
    for(i in 1:length(text2)){
      cat(text2[i],"\n",sep="")
    }
  }

  
  if (x$p.value == 1){
    text <-"The two variances are equal so when doing a two-group t-test, check the equal variances box."
 }
  else { text <-"The two variances are NOT equal so when doing a two-group t-tests, do NOT check the equal variances box."
  }
  
  # text is the test assumption
  text1 <- paste("Test Information: This test determines whether the variances of the two samples are equal by using an F-test.
                \r ****************************************************************
                \n \n",sep="")
  
  wrapper(text1)
  wrapper(text)
}

distdiagnostics <- function(x){
  defaults <- list(initial.group = NULL, initial.response = NULL, initial.alternative = "two.sided", 
                   initial.confidenceLevel = ".95", initial.label=NULL)
  dialog.values <- getDialog("distdiagnostics", defaults)
  initializeDialog(title = gettextRcmdr("Distribution Diagnostics"))
 variablesFrame <- tkframe(top)
 groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), 
                             title = gettextRcmdr("Groups (pick one)"), 
                             initialSelection = varPosn(dialog.values$initial.group, "twoLevelFactor"))
 responseBox <- variableListBox(variablesFrame, Numeric(), 
                                title = gettextRcmdr("Response Variable (pick one)"),
                                initialSelection = varPosn(dialog.values$initial.response, "numeric"))
 onOK <- function() {
   group <- getSelection(groupBox)
   if (length(group) == 0) {
     errorCondition(recall = distdiagnostics, 
                    message = gettextRcmdr("You must select a groups variable."))
     return()
   }
   response <- getSelection(responseBox)
   if (length(response) == 0) {
     errorCondition(recall = distdiagnostics, 
                    message = gettextRcmdr("You must select a response variable."))
     return()
   }
   alternative <- as.character(tclvalue(alternativeVariable))
   level <- tclvalue(confidenceLevel)
   putDialog ("Distribution Diagnostics", list (initial.group = group, initial.response = response, initial.alternative = alternative, 
                                                initial.confidenceLevel = level, initial.label=.groupsLabel))  
    closeDialog()

   doItAndPrint(paste("variancetest <- var.test(", response, "~", group,", data=", ActiveDataSet(),",ratio=1,
                      alternative='",alternative,"', conf.level=",level,")", sep = ""))
   doItAndPrint("variancetest")
   doItAndPrint("DistWords(variancetest)")
   
   # End Insertion
    tkdestroy(top)
    tkfocus(CommanderWindow())
}

OKCancelHelp(helpSubject = "t.test", reset = "distdiagnostics")
optionsFrame <- tkframe(top)
radioButtons(optionsFrame, name = "alternative", buttons = c("twosided", 
                                                             "less", "greater"), values = c("two.sided", "less", "greater"), 
             labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
                                     "Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"),
             initialValue = dialog.values$initial.alternative)
confidenceFrame <- tkframe(optionsFrame)
confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
confidenceField <- ttkentry(confidenceFrame, width = "6", 
                            textvariable = confidenceLevel)
tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text = "    "), 
       getFrame(responseBox), sticky = "nw")
tkgrid(variablesFrame, sticky = "nw")
tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level"), 
                  fg = "blue"), sticky = "w")
tkgrid(confidenceField, sticky = "w")
groupsLabel(groupsBox = groupBox, initialText=dialog.values$initial.label)
tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text = "    "), 
       confidenceFrame, labelRcmdr(optionsFrame, text = "    "), sticky = "nw")
tkgrid(optionsFrame, sticky = "nw")
tkgrid(buttonsFrame, sticky = "w")
dialogSuffix(rows = 4, columns = 1)
}



### Another added function for Distribution Diagnostics
DistWords2 <- function(x,y,z){
  wrapper <- function(text){
    text2 <- strwrap(text)
    for(i in 1:length(text2)){
      cat(text2[i],"\n",sep="")
    }
  }
  clevel <- 1-y
  
  group1 <- z[1]
  group2 <- z[2]
  
  text1 <- paste("Test Information: This test determines whether the data of the two samples are normally distributed.
                  If the p-value is less than the alpha of ",clevel," then the data is not normally distributed. 
                 \r ****************************************************************
                 \n \n",sep="")
  
  wrapper(text1)
}

distnormal <- function(x){
  defaults <- list(initial.group = NULL, initial.response = NULL, initial.confidenceLevel = ".95", initial.label=NULL)
  dialog.values <- getDialog("distdiagnostics", defaults)
  initializeDialog(title = gettextRcmdr("Distribution Diagnostics"))
  variablesFrame <- tkframe(top)
  groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), 
                              title = gettextRcmdr("Groups (pick one)"), 
                              initialSelection = varPosn(dialog.values$initial.group, "twoLevelFactor"))
  responseBox <- variableListBox(variablesFrame, Numeric(), 
                                 title = gettextRcmdr("Response Variable (pick one)"),
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"))
  onOK <- function() {
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall = distdiagnostics, 
                     message = gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall = distdiagnostics, 
                     message = gettextRcmdr("You must select a response variable."))
      return()
    }
    level <- tclvalue(confidenceLevel)
    putDialog ("Distribution Diagnostics", list (initial.group = group, initial.response = response, initial.confidenceLevel = level, initial.label=.groupsLabel))  
    closeDialog()
    
    doItAndPrint(paste("normaltest <- tapply(", ActiveDataSet(),"$", response, ",", ActiveDataSet(),"$", group,", shapiro.test)", sep = ""))
    doItAndPrint("normaltest")
   doItAndPrint(paste("DistWords2(normaltest, ",level,",", ActiveDataSet(),"$", group, ")", sep = ""))
    # End Insertion
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject = "t.test", reset = "distnormal")
  optionsFrame <- tkframe(top)
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel)
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text = "    "), 
         getFrame(responseBox), sticky = "nw")
  tkgrid(variablesFrame, sticky = "nw")
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level"), 
                    fg = "blue"), sticky = "w")
  tkgrid(confidenceField, sticky = "w")
  groupsLabel(groupsBox = groupBox, initialText=dialog.values$initial.label)
  tkgrid( confidenceFrame, labelRcmdr(optionsFrame, text = "    "), sticky = "nw")
  tkgrid(optionsFrame, sticky = "nw")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix(rows = 4, columns = 1)
}
