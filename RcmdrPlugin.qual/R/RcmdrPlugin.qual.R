# Some Rcmdr dialogs for the qual package

# last modified: October 19 2010 by E. Hodgess

# Note: the following function (with contributions from Richard Heiberger) 
# can be included in any Rcmdr plug-in package to cause the package to load
# the Rcmdr if it is not already loaded

.onAttach <- function(libname, pkgname){
        if (!interactive()) return()
 #       Rcmdr <- options()$Rcmdr
 #       plugins <- Rcmdr$plugins
 #       if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
 #               Rcmdr$plugins <- c(plugins, pkgname)
 #               options(Rcmdr=Rcmdr)
 #               closeCommander(ask=FALSE, ask.save=TRUE)
 #               Commander()
 #       }

 Rcmdr <- options()$Rcmdr
  plugins <- Rcmdr$plugins

  if (!pkgname %in% plugins) {
      Rcmdr$plugins <- c(plugins, pkgname)
      Rcmdr$ask.on.exit <- FALSE
      options(Rcmdr=Rcmdr)

      if("package:Rcmdr" %in% search()) {
          if(!getRcmdr("autoRestart")) {
              closeCommander(ask=FALSE, ask.save=TRUE)
              Commander()
          }
      }
      else {
          if(packageVersion("Rcmdr") < "2.0-0")
              stop(.gettext("This package requires Rcmdr version 2.0-0 or higher."))

          Commander()
      }
}
}


    
   


xbara <- function() {
	
.activeDataSet <- ActiveDataSet()
        
 	command <- paste('qcc(',.activeDataSet,',"xbar")',sep="")
	justDoIt(command)
        logger(command)

        tkfocus(CommanderWindow())

            }
  

	


rchart <- function() {


	
.activeDataSet <- ActiveDataSet()
        
 	command <- paste('qcc(',.activeDataSet,',"R")',sep="")
	justDoIt(command)
        logger(command)

        tkfocus(CommanderWindow())

            }
  


schart <- function() {


	
.activeDataSet <- ActiveDataSet()
        
 	command <- paste('qcc(',.activeDataSet,',"S")',sep="")
	justDoIt(command)
        logger(command)

        tkfocus(CommanderWindow())

            }
  



npchart <- function() {
	initializeDialog(title=gettextRcmdr("np Chart"))

  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=npchart, message=gettextRcmdr("You must select a variable."))
            return()
            }

   
        freq1 <- tclvalue(freqVariable)
freq1 <- as.numeric(tclvalue(freqVariable))
        if (is.na(freq1) || freq1 <= 0){
            errorCondition(recall=npchart, message="Length must be a 
			positive integer.")
            return()
            }
	
         closeDialog()

 command <- paste('qcc(',ActiveDataSet(),'$',x,',"np",size=',freq1,')',sep="")
 

	justDoIt(command)
        logger(command)
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
 

    OKCancelHelp(helpSubject="qcc")
    freqFrame <- tkframe(top)
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
      tkgrid(getFrame(xBox), sticky="nw") 
   tkgrid(tklabel(freqFrame, text=gettextRcmdr("How many? ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(freqField, sticky="e")
 tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
 	
    dialogSuffix(rows=4, columns=2)
    }











pchart <- function() {
	initializeDialog(title=gettextRcmdr("p Chart"))


   xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=pchart, message=gettextRcmdr("You must select a variable."))
            return()
            }

   
        freq1 <- tclvalue(freqVariable)
freq1 <- as.numeric(tclvalue(freqVariable))
        if (is.na(freq1) || freq1 <= 0){
            errorCondition(recall=pchart, message="Length must be a 
			positive integer.")
            return()
            }
	
         closeDialog()

 command <- paste('qcc(',ActiveDataSet(),'$',x,',"p",size=',freq1,')',sep="")
 

	justDoIt(command)
        logger(command)
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="qcc")
    freqFrame <- tkframe(top)
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
      tkgrid(getFrame(xBox), sticky="nw") 
   tkgrid(tklabel(freqFrame, text=gettextRcmdr("How many? ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(freqField, sticky="e")
 tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
 	
    dialogSuffix(rows=4, columns=2)
    }






  

shapeMod <- function(){
    initializeDialog(title=gettextRcmdr("Shape a data frame"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))

    dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
    
    onOK <- function(){
  dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=shapeMod, 
                message=gettextRcmdr("You must enter the name of a data set."))  
            return()
            }  
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall=newDataSet,
                message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
                newDataSet()
                return()
	}
	}
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=shapeMod, message=gettextRcmdr("You must select a variable."))
            return()
            }

   
        freq1 <- tclvalue(freqVariable)
freq1 <- as.numeric(tclvalue(freqVariable))
        if (is.na(freq1) || freq1 <= 0){
            errorCondition(recall=pchart, message="Length must be a 
			positive integer.")
            return()
            }
	
	command <- paste(dsnameValue, "<- shape1(", ActiveDataSet(), "$", x,
         ",freq=",freq1,
            ")", sep="")
	justDoIt(command)
        logger(command)
activeDataSet(dsnameValue)

         closeDialog()
	}
 OKCancelHelp(helpSubject="qcc")
       freqFrame <- tkframe(top)
   tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
      tkgrid(getFrame(xBox), sticky="nw") 
   tkgrid(tklabel(freqFrame, text=gettextRcmdr("Number of rows: ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
      tkgrid(getFrame(xBox), sticky="nw") 

   tkgrid(buttonsFrame, columnspan=2, sticky="w")

 tkgrid.configure(entryDsname, sticky="w")
    tkgrid.configure(freqField, sticky="e")

  
 	
    dialogSuffix(rows=4, columns=2)
    }





shape1 <- function(x,freq) {
	reshape1<-  function(df,nrow,ncol,byrow=TRUE)data.frame(matrix(as.matrix(df),nrow,ncol,byrow=byrow))
	n1 <- length(x)
	nc1 <- floor(n1/freq)
	y <- reshape1(df=x,nrow=freq,ncol=nc1)
	}

paretoMod <- function(){
    initializeDialog(title=gettextRcmdr("Pareto Chart"))
    xBox <- variableListBox(top, Factors(), title=gettextRcmdr("Variable (pick one)"))
 
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=paretoMod, message=gettextRcmdr("You must select a variable."))
            return()
            }

   
	
          closeDialog()
	command <- paste("pareto1(", ActiveDataSet(), "$", x,
                   ")", sep="")
	justDoIt(command)
        logger(command)
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="barplot")
      tkgrid(getFrame(xBox), sticky="nw") 

   tkgrid(buttonsFrame, columnspan=2, sticky="w")
       dialogSuffix(rows=4, columns=2)
    }


pareto1 <- function(x) {
	y1 <- table(x)
	barplot(rev(y1[order(y1)]))
	abline(h=0)
}


 
xbaro <- function() {

.activeDataSet <- ActiveDataSet()
        
 	command <- paste('qcc(',.activeDataSet,',"xbar.one")',sep="")
	justDoIt(command)
        logger(command)

        tkfocus(CommanderWindow())

            }
  
 

 

ewMod <- function(){
  initializeDialog(title=gettextRcmdr("EWMA"))

    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
	print(x)
        if (length(x) == 0){
            errorCondition(recall=ewMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("ewma(", ActiveDataSet(), "$", x,")",sep=""))     
	
             tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ewma")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }



hist1Mod <- function() {
	
.activeDataSet <- ActiveDataSet()


      
 	command <- paste('hist(c(as.matrix(',.activeDataSet,')))',sep="")
	justDoIt(command)
        logger(command)

        tkfocus(CommanderWindow())
}





sum1Mod <- function() {
	
.activeDataSet <- ActiveDataSet()


    	command <- paste('sum1a(c(as.matrix(',.activeDataSet,')))',sep="")
        doItAndPrint(command)	
 
   

        tkfocus(CommanderWindow())
}


sum1a <- function(x) {
	z <- sd(x)
	names(z) <- "sd"
	y <- c(summary(x),z)
	return(y)
	}


pschart <- function() {
	initializeDialog(title=gettextRcmdr("p Chart for different size samples"))

     variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, 
        title=gettextRcmdr("Select attribute variable)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Select size variable"))
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (0 == length(y)) {
            errorCondition(recall=pschart, message=gettextRcmdr("You must select an attribute variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=pschart, message=gettextRcmdr("You must select a size variable."))
            return()
            }
        if (is.element(y, x)) {
            errorCondition(recall=pschart, message=gettextRcmdr("Attribute and size variables must be different."))
            return()
            }
    
        command <- paste('qcc(',ActiveDataSet(),'$',x,',"p",size=',ActiveDataSet(),'$',y,')',sep="")
 	justDoIt(command)
        logger(command)

        closeDialog()
  tkdestroy(top)
       tkfocus(CommanderWindow())

            }
    OKCancelHelp(helpSubject="qcc")
  
    tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
    dialogSuffix(rows=4, columns=2)
    }



npschart <- function() {
	initializeDialog(title=gettextRcmdr("np Chart for different size samples"))

     variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, 
        title=gettextRcmdr("Select attribute variable)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Select size variable"))
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (0 == length(y)) {
            errorCondition(recall=npschart, message=gettextRcmdr("You must select an attribute variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=npschart, message=gettextRcmdr("You must select a size variable."))
            return()
            }
        if (is.element(y, x)) {
            errorCondition(recall=npschart, message=gettextRcmdr("Attribute and size variables must be different."))
            return()
            }
    
        command <- paste('qcc(',ActiveDataSet(),'$',x,',"np",size=',ActiveDataSet(),'$',y,')',sep="")
 	justDoIt(command)
        logger(command)

        closeDialog()
  tkdestroy(top)
       tkfocus(CommanderWindow())

            }
    OKCancelHelp(helpSubject="qcc")
  
    tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
    dialogSuffix(rows=4, columns=2)
    }


cchart <- function() {
	initializeDialog(title=gettextRcmdr("c Chart"))

.activeDataSet <- ActiveDataSet()
  nVar <- tclVar("1")
    nEntry <- tkentry(top, width="6", textvariable=nVar)
  
  
    onOK <- function(){
  n <- round(as.numeric(tclvalue(nVar)))
        if (is.na(n) || n <= 0){
            errorCondition(recall=cchart, message="Length must be a 
			positive integer.")
            return()
            }
     
        command <- paste('qcc(',.activeDataSet,',"c",size=',n,')',sep="")
 	justDoIt(command)
        logger(command)

        closeDialog()
        tkfocus(CommanderWindow())

            }
    OKCancelHelp(helpSubject="qcc")
  
  rightFrame <- tkframe(top)
    tkgrid(tklabel(top, text="How many?"), nEntry, sticky="e")
  
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
 tkgrid.configure(nEntry, sticky="w")
    dialogSuffix(rows=10, columns=2, focus=nEntry)
    }


uchart <- function() {
	initializeDialog(title=gettextRcmdr("u Chart"))

.activeDataSet <- ActiveDataSet()
  nVar <- tclVar("1")
    nEntry <- tkentry(top, width="6", textvariable=nVar)
  
  
    onOK <- function(){
  n <- round(as.numeric(tclvalue(nVar)))
        if (is.na(n) || n <= 0){
            errorCondition(recall=uchart, message="Length must be a 
			positive integer.")
            return()
            }
     
        command <- paste('qcc(',.activeDataSet,',"u",size=',n,')',sep="")
 	justDoIt(command)
        logger(command)

        closeDialog()
        tkfocus(CommanderWindow())

            }
    OKCancelHelp(helpSubject="qcc")
  
  rightFrame <- tkframe(top)
    tkgrid(tklabel(top, text="How many?"), nEntry, sticky="e")
  
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
 tkgrid.configure(nEntry, sticky="w")
    dialogSuffix(rows=10, columns=2, focus=nEntry)
    }



uschart <- function() {
	initializeDialog(title=gettextRcmdr("u Chart for different size samples"))

     variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, 
        title=gettextRcmdr("Select attribute variable)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Select size variable"))
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (0 == length(y)) {
            errorCondition(recall=uschart, message=gettextRcmdr("You must select an attribute variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=uschart, message=gettextRcmdr("You must select a size variable."))
            return()
            }
        if (is.element(y, x)) {
            errorCondition(recall=uschart, message=gettextRcmdr("Attribute and size variables must be different."))
            return()
            }
    
        command <- paste('qcc(',ActiveDataSet(),'$',x,',"u",size=',ActiveDataSet(),'$',y,')',sep="")
 	justDoIt(command)
        logger(command)

        closeDialog()
  tkdestroy(top)
       tkfocus(CommanderWindow())

            }
    OKCancelHelp(helpSubject="qcc")
  
    tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
    dialogSuffix(rows=4, columns=2)
    }


cusumMod <- function(){
  initializeDialog(title=gettextRcmdr("CUSUM chart"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
	print(x)
        if (length(x) == 0){
            errorCondition(recall=cusumMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("cusum(", ActiveDataSet(), "$", x,",plot=TRUE)",sep=""))     
	
             tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="cusum")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }




movranMod <- function(){
  initializeDialog(title=gettextRcmdr("Moving Range chart"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
	print(x)
        if (length(x) == 0){
            errorCondition(recall=cusumMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("movrange1(", ActiveDataSet(), "$", x,")",sep=""))     
	
             tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="qcc")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }




movrange1 <- function(x) {
	   	  
	  n1 <- length(x)
	  y <- numeric(length=(n1-1))
	  for(i in 1:(n1-1)) {
	  	y[i] <- abs(x[i+1]-x[i])
		}
	  ym <- mean(y)
	  ysd <- ym/1.128
	  d3 <- qcc.options("se.R.unscaled")
	  ucl <- ym + 3*d3[n1]*ysd
	  lcl <- max(0,ym-3*d3[n1]*ysd)
	  cc <- ifelse(y>=ucl,"red","black")
	  plot(y,ylim=c(0,max(ucl,y)),type="b",pch=16,col=cc,main="Moving Range Chart")
	  abline(h=ym)
	  abline(h=ucl,lty=2)
	  abline(h=lcl,lty=2)
	  }

gofMod <- function(){
    initializeDialog(title=gettextRcmdr("Goodness of Fit for a Variable"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=gofMod, message=gettextRcmdr("You must select a variable."))
            return()
            }

  
          closeDialog()
	command <- paste("gof1(", ActiveDataSet(), "$", x,
            ")", sep="")

	     doItAndPrint(command)



        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="data.frame")
      tkgrid(getFrame(xBox), sticky="nw") 
  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
 tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
	
    dialogSuffix(rows=4, columns=2)
    }





gof1 <- function(x) {
        if(all(x>0)) {
	xw <- fitdistr(x,densfun=dweibull,start=list(scale=1,shape=2))
	xl <- fitdistr(x,densfun="lognormal")
	xkw <- ks.test(x,"pweibull",scale=xw[[1]][1],shape=xw[[1]][2])
	xkl <- ks.test(x,"plnorm",xl[[1]][1],xl[[1]][2])
	}
	xn <- fitdistr(x,densfun="normal")
	xkn <- ks.test(x,"pnorm",xn[[1]][1],xn[[1]][2])
	if(all(x>0)){
	p1 <- max(xkw$p,xkn$p,xkl$p)
	if(p1==xkw$p)z <- list(p=p1,dist="Weibull")
	if(p1==xkn$p)z <- list(p=p1,dist="Normal")
	if(p1==xkl$p)z <- list(p=p1,dist="Lognormal")
	}
	else {
	z <- list(p=xkn$p,dist="Normal")
	}
	return(z)

}
xbarnew  <- function () {
	 dataSets <- listDataSets()
	 	initializeDialog(title=gettextRcmdr("Select Data Sets"))
	dataSetsBox1 <- variableListBox(top, dataSets, title=gettextRcmdr("Data Set One  (pick one)"))
	dataSetsBox2 <- variableListBox(top, dataSets, title=gettextRcmdr("Data Set Two  (pick one)"))


	onOK <- function(){
	x <- getSelection(dataSetsBox1)
        y <- getSelection(dataSetsBox2)
		if (length(x) == 0 | length(y) == 0) {
			errorCondition(recall = xbarnew, message = gettextRcmdr("You must select two Sets."))
				return()
}
	if(x==y) {
		 errorCondition(recall=xbarnew,message= gettextRcmdr("The two sets must be different"))
		 x <- NULL
		 y <- NULL
		 }



	closeDialog()

 	command <- paste('qcc(',x,',"xbar",newdata=',y,')',sep="")
		doItAndPrint(command)
		tkdestroy(top)


		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="qcc")
	
	tkgrid(getFrame(dataSetsBox1), getFrame(dataSetsBox2),sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=3, columns=2)
}


cschart <- function() {
	initializeDialog(title=gettextRcmdr("c Chart for different size samples"))

     variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, 
        title=gettextRcmdr("Select attribute variable)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Select size variable"))
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (0 == length(y)) {
            errorCondition(recall=cschart, message=gettextRcmdr("You must select an attribute variable."))
            return()
            }
        if (0 == length(x)) {
            errorCondition(recall=cschart, message=gettextRcmdr("You must select a size variable."))
            return()
            }
        if (is.element(y, x)) {
            errorCondition(recall=cschart, message=gettextRcmdr("Attribute and size variables must be different."))
            return()
            }
    
        command <- paste('qcc(',ActiveDataSet(),'$',x,',"c",size=',ActiveDataSet(),'$',y,')',sep="")
 	justDoIt(command)
        logger(command)

        closeDialog()
  tkdestroy(top)
       tkfocus(CommanderWindow())

            }
    OKCancelHelp(helpSubject="qcc")
  
    tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
    dialogSuffix(rows=4, columns=2)
    }

