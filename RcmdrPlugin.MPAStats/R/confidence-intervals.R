# Modified on March 22, 2015 by Jessica Peterson

# Extract and pretty print the confidence interval from an object of class `htest`
printConfint <- function(x) {
  if (class(x) != "htest") 
    stop("must supply an object of class htest")
  
  cat(format(100 * attr(x$conf.int, "conf.level")),
      "percent confidence interval:\n",
      format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
}


#Interpretation Function
confintContinuousWords <- function(x){
     varcut<-nchar(ActiveDataSet())+2

     cat("Our sample suggests that the true average (mean) of ",
       tolower(substring(x$data.name,varcut))," is between"," \n ",
       round(x$conf.int[1L],2)," and ",round(x$conf.int[2L],2), " (",
       100*attr(x$conf.int,"conf.level"),"% confidence level). \n",sep="" )
}


# Modified from confintContinuous function
confintContinuous2 <- function () {
  defaults <- list (initial.x = NULL, initial.level = ".95")
  dialog.values <- getDialog ("confintContinuous2", defaults)  
  initializeDialog(title = gettextRcmdr("Confidence intervals for continuous data"))
  xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"),
                          initialSelection = varPosn(dialog.values$initial.x, "numeric"))
  onOK <- function() {
    x <- getSelection(xBox)
    if (length(x) == 0) {
      errorCondition(recall = confintContinuous2, 
        message = gettextRcmdr("You must select a variable."))
      return()
    }
    level <- tclvalue(confidenceLevel)
    putDialog ("confintContinuous2", list (initial.x = x, initial.level = level))
    closeDialog()
    doItAndPrint(paste(".test <- t.test(", ActiveDataSet(), "$", x, 
                       ", conf.level=", level, ")", sep = ""))
    doItAndPrint("printConfint(.test)")
    #Insertion:
    doItAndPrint("confintContinuousWords(.test)")
    #End Insertion
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "t.test", reset = "confintContinuous2")
  
  # Create main frames
  leftFrame <- getFrame(xBox)
  rightFrame <- tkframe(top)

  # Confidence frame
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(dialog.values$initial.level)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel)

  # Labels
  tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: ")), 
                    confidenceField, sticky = "w")

  # Place frames
  tkgrid(leftFrame, rightFrame)
  tkgrid(confidenceFrame, sticky = "w")

  # Set anchoring options
  tkgrid.configure(confidenceField, sticky = "e")
  tkgrid.configure(rightFrame, sticky="nw")
  tkgrid.configure(leftFrame, sticky="nw")

  # Final buttons
  tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
  dialogSuffix(rows = 4, columns = 2)
}


# Dialog for calculating confidence intervals for Poisson data (NOT USED PRESENTLY?)
confintPoisson <- function () {
  defaults <- list (initial.x = NULL, initial.level = ".95")
  dialog.values <- getDialog ("confintPoisson", defaults)  
  initializeDialog(title = gettextRcmdr("Confidence intervals for Poisson data"))
  xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"),
                          initialSelection = varPosn(dialog.values$initial.x, "numeric"))
  onOK <- function() {
    x <- getSelection(xBox)
    if (length(x) == 0) {
      errorCondition(recall = confintPoisson, 
        message = gettextRcmdr("You must select a variable."))
      return()
    }
    level <- tclvalue(confidenceLevel)
    putDialog ("confintPoisson", list (initial.x = x, initial.level = level))
    closeDialog()
    doItAndPrint(paste(".test.poisson <- poisson.test(sum(", ActiveDataSet(), "$", x,
                       "), length(", ActiveDataSet(), "$", x,"), conf.level=", 
                       level, ")", sep = ""))
    doItAndPrint("printConfint(.test.poisson)")
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "t.test", reset = "confintPoisson")
  
  # Create main frames
  leftFrame <- getFrame(xBox)
  rightFrame <- tkframe(top)

  # Confidence frame
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(dialog.values$initial.level)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel)

  # Labels
  tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: ")), 
                    confidenceField, sticky = "w")

  # Place frames
  tkgrid(leftFrame, rightFrame)
  tkgrid(confidenceFrame, sticky = "w")

  # Set anchoring options
  tkgrid.configure(confidenceField, sticky = "e")
  tkgrid.configure(rightFrame, sticky="nw")
  tkgrid.configure(leftFrame, sticky="nw")

  # Final buttons
  tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
  dialogSuffix(rows = 4, columns = 2)
}

#Interpretation Function
confintBinomialWords <- function(level,varname,x){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }

    l.limit <- x$conf.int[1]
    u.limit <- x$conf.int[2]
    conf.level <- 100*level
    text <- paste("Our sample suggests that the proportion of ",varname," in the population is between ",round(l.limit,4)," and ",round(u.limit,4)," (",conf.level,"% confidence level, binomial corrected). \n \n",sep="")	
    wrapper(text)
}


# Dialog for calculating binary data (modified from confintBinomial function)
confintBinomial2 <- function () {
  defaults <- list (initial.x = NULL,  initial.level = ".95") 
  dialog.values <- getDialog ("confintBinomial2", defaults)  
  initializeDialog(title = gettextRcmdr("Confidence intervals for binomial data"))

 xBox <- variableListBox(top, TwoLevelFactors(), title = gettextRcmdr("Variable (pick one)"),
                        initialSelection = varPosn(dialog.values$initial.x,"factor"))
	
  onOK <- function() {
    x <- getSelection(xBox)
		if (length(x) == 0) {
			errorCondition(recall = singleProportionTest, message = gettextRcmdr("You must select a variable."))
			return()
		}


    level <- tclvalue(confidenceLevel)
    alternative <- "two.sided"
    p <- .5
    ### INSERTED CODE
    ### assign("x",x,envir=SUB)
    ### assign("level",level,envir=SUB)
    ### assign("alternative",alternative,envir=SUB)
    ### assign("p",p,envir=SUB)
    ### END INSERTION
 
  putDialog ("confintBinomial2", list (initial.x = x, initial.level = level))
 
    closeDialog()
    # Modified to interact with the SUB environment

var1 <- paste("table(",ActiveDataSet(),"$",x,")[1]")
var2 <- paste("table(",ActiveDataSet(),"$",x,")[2]")
command <- paste("c(",var2,",",var1,")")
#changed command from xtabs to c(var1,var2) in preparation to be able to choose which factor you want to build on
   # command <- paste("xtabs(~", x, ", data=", ActiveDataSet(), ")")
                # logger(paste(".Table <-", command))
		### assign(".Table", justDoIt(with(environment(),command)),envir=SUB)
		doItAndPrint(paste(".Table","<-",command))
		### doItAndPrint("with(SUB,.Table)")
		doItAndPrint(".Table")
    #Creates (and prints) the matrix "ntable" which has the levels and counts of the variable 
    doItAndPrint(paste("ntable<-(rbind(.Table))"))
    #varname takes the 1st column name (or 1st level) 
    varname <- colnames(ntable)[1]
    varname2 <- colnames(ntable)[2]
		### doItAndPrint(paste("with(SUB,.test.bi <- binom.test(rbind(.Table), alternative='", 
		###					alternative, "', p=", 
                ###                                        p, ", conf.level=", 
                ###                                        level, 
		###					"))", sep = "")) 
		doItAndPrint(paste(".test.bi <- binom.test(rbind(.Table), alternative='", 
							alternative, "', p=", 
                                                        p, ", conf.level=", 
                                                        level, 
							")", sep = "")) 
    ### doItAndPrint("with(SUB,printConfint(.test.bi))")
    doItAndPrint("printConfint(.test.bi)")
    # Inserted Code
    ### doItAndPrint(paste("with(SUB,confintBinomialWords(",level,",",'"',x,'"',",.test.bi))",sep=""))
    doItAndPrint(paste("confintBinomialWords(",level,",",'"',varname,'"',",.test.bi)",sep=""))
    doItAndPrint(paste("confintBinomialWords(",level,",",'"',varname2,'"',",.test.bi)",sep=""))
    # End Insertion  
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "t.test", reset = "confintBinomial2")
  
  # Create main frames
  leftFrame <- getFrame(xBox)
  rightFrame <- tkframe(top)
 
  # Confidence frame
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(dialog.values$initial.level)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel) 
  # Labels
  tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: ")), 
                    confidenceField, sticky = "w")


  # Place frames
  tkgrid(leftFrame, rightFrame)
  tkgrid(confidenceFrame, sticky = "w")

  # Set anchoring options
  tkgrid.configure(confidenceField, sticky = "e")
  tkgrid.configure(rightFrame, sticky="nw")
  tkgrid.configure(leftFrame, sticky="nw")

  # Final buttons
  tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
  dialogSuffix(rows = 4, columns = 2)
}
