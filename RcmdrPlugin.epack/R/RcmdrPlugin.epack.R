# Some Rcmdr dialogs for the epack package

# last modified: March 31, 2012  by E. Hodgess

# Note: the following function (with contributions from Richard Heiberger) 
# can be included in any Rcmdr plug-in package to cause the package to load
# the Rcmdr if it is not already loaded


.onAttach <- function(libname, pkgname){
        if (!interactive()) return()
        Rcmdr <- options()$Rcmdr
        plugins <- Rcmdr$plugins
        if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
                Rcmdr$plugins <- c(plugins, pkgname)
                options(Rcmdr=Rcmdr)
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
        }
}


    
   




 


ArimaMod <- function(){
 initializeDialog(title=gettextRcmdr("ARIMA Models"))
   .activeModel <- ActiveModel()

    currentModel <- if (!is.null(.activeModel)) 
        eval(parse(text=paste("class(", .activeModel, ")[1] == 'Arima'", sep="")), 
            envir=.GlobalEnv) 
        else FALSE
#    if (currentModel) {
#        currentFields <- formulaFields(eval(parse(text=.activeModel), 
#            envir=.GlobalEnv))
#        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#        }
    UpdateModelNumber()
   modelName <- tclVar(paste("ArimaModel.", getRcmdr("modelNumber"), sep=""))
   modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=modelName)
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
          if (!is.valid.name(modelValue)){
            errorCondition(recall=ArimaMod, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
	         x <- getSelection(xBox)
       if (length(x) == 0){
            errorCondition(recall=ArimaMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
         AR1 <- tclvalue(ARLevel)
        MA1 <- tclvalue(MAVariable)
	DIF1 <- tclvalue(DIFVariable)
    
	SAR1 <- tclvalue(SARLevel)
        SMA1 <- tclvalue(SMAVariable)
	SDIF1 <- tclvalue(SDIFVariable)
 
   alternative <- as.character(tclvalue(alternativeVariable))
 	 closeDialog() 
    
       if (is.element(modelValue, listArimaModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                ArimaMod()
                return()
                }
            }
         
      pd1 <- paste("tsp(",ActiveDataSet(),"$",x,")[3]",sep="")

	assign("pd2",justDoIt(pd1),envir=.GlobalEnv)
    command <- paste("Arima(", ActiveDataSet(), "$", x,",order=c(",AR1,",",DIF1,",",MA1,"),
		include.mean=",alternative,",seasonal=list(order=c(",SAR1,",",SDIF1,",",SMA1
		,"),period=",pd2,"))",sep="")
    
 
  logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(modelValue)
           activeModel(modelValue)
    

         tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="arima")
   radioButtons(top, name="alternative", buttons=c("yes","no"), values=c(TRUE,FALSE),
        labels=gettextRcmdr(c("Mean Included","Mean Not Included")), 
        title=gettextRcmdr("Mean"))
      rightFrame <- tkframe(top)
    ARFrame <- tkframe(rightFrame)
    ARLevel <- tclVar("0")
    ARField <- tkentry(ARFrame, width="6", textvariable=ARLevel)
    MAFrame <- tkframe(rightFrame)
    MAVariable <- tclVar("0")
    MAField <- tkentry(MAFrame, width="6", textvariable=MAVariable)
    DIFFrame <- tkframe(rightFrame)
    DIFVariable <- tclVar("0")
    DIFField <- tkentry(DIFFrame,width="6",textvariable=DIFVariable)
 
    SARFrame <- tkframe(rightFrame)
    SARLevel <- tclVar("0")
    SARField <- tkentry(SARFrame, width="6", textvariable=SARLevel)
    SMAFrame <- tkframe(rightFrame)
    SMAVariable <- tclVar("0")
    SMAField <- tkentry(SMAFrame, width="6", textvariable=SMAVariable)
    SDIFFrame <- tkframe(rightFrame)
    SDIFVariable <- tclVar("0")
    SDIFField <- tkentry(SDIFFrame,width="6",textvariable=SDIFVariable)


 tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
     tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(rightFrame, text=""), sticky="w")   
   tkgrid(tklabel(ARFrame, text=gettextRcmdr("Reg. AR Order ")), ARField, sticky="w")
    tkgrid(ARFrame, sticky="w")
    tkgrid(tklabel(DIFFrame,text=gettextRcmdr("Reg. Diff Order ")),DIFField, sticky="w")
    tkgrid(DIFFrame,sticky="w")
    tkgrid(tklabel(MAFrame, text=gettextRcmdr("Reg. MA Order")), MAField, sticky="w")
    tkgrid(MAFrame, sticky="w")
  
   tkgrid(tklabel(SARFrame, text=gettextRcmdr("Sea. AR Order ")), SARField, sticky="w")
    tkgrid(SARFrame, sticky="w")
    tkgrid(tklabel(SDIFFrame,text=gettextRcmdr("Sea. Diff Order ")),SDIFField, sticky="w")
    tkgrid(SDIFFrame,sticky="w")
    tkgrid(tklabel(SMAFrame, text=gettextRcmdr("Sea. MA Order")), SMAField, sticky="w")
    tkgrid(SMAFrame, sticky="w")
 

   tkgrid(alternativeFrame, rightFrame, sticky="nw")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(ARField, sticky="e")
    dialogSuffix(rows=6, columns=2)
    }



bc2 <- function(x,lam1=seq(-2,2,.1)) {
  require(MASS)
  if(any(x)<=0)stop("Negative values present..no transformations permitted")
  t1 <- 1:length(x)
  xxx <-boxcox(x~t1,lambda=lam1,plot=FALSE)
  z <- xxx$x[xxx$y==max(xxx$y)]
return(z) 
}



newHistPrice <- function() {
    initializeDialog(title=gettextRcmdr("New Historical Price"))
    dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
    tick <- tclVar("ibm")
    ticker <- tkentry(top, width="8", textvariable=tick)
    starty <- tclVar("1998")
    stay1 <- tkentry(top,width="4",textvariable=starty)
    startm <- tclVar("01")
    stam1 <- tkentry(top,width="2",textvariable=startm)
    startd <- tclVar("01")
    stad1 <- tkentry(top,width="2",textvariable=startd)
    endy <- tclVar("1998")
    endy1 <- tkentry(top,width="4",textvariable=endy)
   
    endm <- tclVar("01")
    endm1 <- tkentry(top,width="2",textvariable=endm)
  
    endd <- tclVar("31")
   endd1 <- tkentry(top,width="2",textvariable=endd)
  
    onOK <- function(){
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=newHistPrice, 
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
       alternative <- as.character(tclvalue(alternativeVariable))
       startyVal <- tclvalue(starty)
       startmVal <- tclvalue(startm)
       startdVal <- tclvalue(startd)
   endyVal <- tclvalue(endy)
       endmVal <- tclvalue(endm)
       enddVal <- tclvalue(endd)
    
       TickValue <- trim.blanks(tclvalue(tick))
     if (TickValue == "") {
            errorCondition(recall=newHistPrice, 
                message=gettextRcmdr("You must enter a ticker symbol."))  
            return()
            }  
	tick2 <- paste("",TickValue,"",sep="")
	alt1 <- paste("",alternative,"",sep="")
        startb <- paste(startyVal,"-",startmVal,"-",startdVal,sep="")
        starta <- paste("",startb,"",sep="")
        endb <- paste(endyVal,"-",endmVal,"-",enddVal,sep="")
        enda <- paste("",endb,"",sep="")
	command <- paste('histprice2(inst1="',tick2,'",quot1="',alt1,'",start1="',starta,'",end1="',enda,'")',
	sep="")
	assign(dsnameValue, justDoIt(command), envir=.GlobalEnv)

        logger(paste(dsnameValue, "<-", command))
      if (eval(parse(text=paste("nrow(", dsnameValue, ")"))) == 0){
            errorCondition(recall=newHistPrice, message=gettextRcmdr("empty data set."))
            return()
            }
        activeDataSet(dsnameValue)
        closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="get.hist.quote")
    radioButtons(top, name="alternative", buttons=c("Open", "Close","High","Low","AdjClose"), 
	values=c("Open", "Close","High","Low","AdjClose"),
        labels=gettextRcmdr(c("Obtain Open Price", "Obtain Closing Price",
	"Obtain High Price","Obtain Low Price","Obtain Adj Closing Price")), 
        title=gettextRcmdr("Price Options"))
  
  rightFrame <- tkframe(top)
    tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
    tkgrid(alternativeFrame, rightFrame, sticky="nw")
    tkgrid(tklabel(top,text="Ticker"),ticker,sticky="w")
    tkgrid(tklabel(top,text="Start year"),stay1,sticky="w")
    tkgrid(tklabel(top,text="Start month"),stam1,sticky="w")
    tkgrid(tklabel(top,text="Start day"),stad1,sticky="w")
  tkgrid(tklabel(top,text="End year"),endy1,sticky="w")
    tkgrid(tklabel(top,text="End month"),endm1,sticky="w")
    tkgrid(tklabel(top,text="End day"),endd1,sticky="w")
   
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
    tkgrid.configure(entryDsname, sticky="w")
	tkgrid.configure(ticker,sticky="w")
tkgrid.configure(stay1,sticky="w")
tkgrid.configure(stam1,sticky="w")
tkgrid.configure(stad1,sticky="w")
tkgrid.configure(endy1,sticky="w")
tkgrid.configure(endm1,sticky="w")
tkgrid.configure(endd1,sticky="w")

    dialogSuffix(rows=10, columns=2, focus=entryDsname)
    }





bcMod <- function(){
  initializeDialog(title=gettextRcmdr("Box Cox Transformation"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
     if (length(x) == 0){
            errorCondition(recall=bcMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("bc2(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="boxcox")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }



opMod <- function(){
  initializeDialog(title=gettextRcmdr("Original Plots"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
	print(x)
        if (length(x) == 0){
            errorCondition(recall=opMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("ts.plot(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="boxcox")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


acfMod <- function(){
  initializeDialog(title=gettextRcmdr("ACF Plots"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
	print(x)
        if (length(x) == 0){
            errorCondition(recall=acfMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("acf(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="acf")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


pacfMod <- function(){
  initializeDialog(title=gettextRcmdr("PACF Plots"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
   if (length(x) == 0){
            errorCondition(recall=pacfMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("pacf(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="pacf")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


decomMod <- function(){
    initializeDialog(title=gettextRcmdr("Multiplicative Decomposition"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
  
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=decomMod, message=gettextRcmdr("You must select a variable."))
            return()
            }

           dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=editTSframe, 
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
    return() 
	}   
        freq1 <- tclvalue(freqVariable)
        fore1 <- tclvalue(foreVariable)
        closeDialog()
   
     comm1 <- paste(dsnameValue,"<-decom5(", ActiveDataSet(), "$", x,
            ",se1=",freq1,",fore1=",fore1,
            ")", sep="")
   	doItAndPrint(comm1)
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="t.test")
    freqFrame <- tkframe(top)
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
    foreFrame <- tkframe(top)
    foreVariable <- tclVar("0")
    foreField <- tkentry(foreFrame, width="8", textvariable=foreVariable)
    tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(foreFrame, text=gettextRcmdr("Number of Forecasts: ")), foreField, sticky="w")
    tkgrid(foreFrame, sticky="w")
    tkgrid(tklabel(freqFrame, text=gettextRcmdr("Freq: ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
   tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  
    tkgrid.configure(entryDsname, sticky="w")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(freqField, sticky="e")
    dialogSuffix(rows=4, columns=2)
    }



specMod <- function(){
  initializeDialog(title=gettextRcmdr("Spectral Density Plots"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
      if (length(x) == 0){
            errorCondition(recall=specMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("spectrum(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="spectrum")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


listArimaModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects, 
        function(.x) "Arima" == (class(eval(parse(text=.x), envir=envir))[1]))]
    }



listGARCHModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects, 
        function(.x) "garch" == (class(eval(parse(text=.x), envir=envir))[1]))]
    }


selectActiveARModel <- function(){
    models <- listAllModels()
    .activeModel <- ActiveModel()
   if ((length(models) == 1) && !is.null(.activeModel)) {
       Message(message=gettextRcmdr("There is only one model in memory."),
               type="warning")
       tkfocus(CommanderWindow())
       return()
       }
    if (length(models) == 0){
        Message(message=gettextRcmdr("There are no models from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
    initializeDialog(title=gettextRcmdr("Select ARIMA Model"))
    .activeDataSet <- ActiveDataSet()
    initial <- if (is.null(.activeModel)) NULL else which(.activeModel == models) - 1
    modelsBox <- variableListBox(top, models, title=gettextRcmdr("Models (pick one)"), 
        initialSelection=initial)
    onOK <- function(){
        model <- getSelection(modelsBox)
        closeDialog()
        if (length(model) == 0) {
            tkfocus(CommanderWindow())
            return()
	}
  #    dataSet <- eval(parse(text=paste("as.character(", model, "$call$data)")))
     dataSet <- eval(parse(text=paste("as.character(", model, "$series)")))
 	dataset1 <- unlist(strsplit(dataSet,"$",fixed=TRUE))[1]
	dataSet <- dataset1
  if (length(dataSet) == 0){
            errorCondition(message=gettextRcmdr("There is no dataset associated with this model."))
            return()
            }
        dataSets <- listDataSets()
        if (!is.element(dataSet, dataSets)){
            errorCondition(message=sprintf(gettextRcmdr("The dataset associated with this model, %s, is not in memory."), dataSet))
            return()
            }
        if (is.null(.activeDataSet) || (dataSet != .activeDataSet)) activeDataSet(dataSet)
      activeModel(model)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp()
    nameFrame <- tkframe(top)
    tkgrid(tklabel(nameFrame, fg="blue", text=gettextRcmdr("Current Model: ")), 
        tklabel(nameFrame, text=tclvalue(getRcmdr("modelName"))), sticky="w")
    tkgrid(nameFrame, sticky="w", columnspan="2")
    tkgrid(getFrame(modelsBox), columnspan="2", sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }



predARModel <- function(){
    models <- listAllModels()
     .activeModel <- ActiveModel()
    if ((length(models) == 1) && !is.null(.activeModel)) {
        Message(message=gettextRcmdr("There is only one model in memory."),
                type="warning")
        tkfocus(CommanderWindow())
        return()
        }
    if (length(models) == 0){
        Message(message=gettextRcmdr("There are no models from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
 dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
  
          dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=decomaMod, 
                message=gettextRcmdr("You must enter the name of a data set."))  
            return()
            }  
        if (!is.valid.name(dsnameValue)) {
            errorCondition(recall=newDataSet,
                message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
    initializeDialog(title=gettextRcmdr("ARIMA Forecasting"))
    .activeDataSet <- ActiveDataSet()
    initial <- if (is.null(.activeModel)) NULL else which(.activeModel == models) - 1
    modelsBox <- variableListBox(top, models, title=gettextRcmdr("Models (pick one)"), 
        initialSelection=initial)
    onOK <- function(){
        model <- getSelection(modelsBox)
        closeDialog()
        if (length(model) == 0) {
            tkfocus(CommanderWindow())
            return()
	}

     activeModel(model)
     fore1 <- tclvalue(foreVariable)
        closeDialog()
        doItAndPrint(paste("forecast.Arima(", model,",fore1=",fore1,
            ")", sep=""))
        tkdestroy(top)
      tkfocus(CommanderWindow())
        }
    OKCancelHelp()
    nameFrame <- tkframe(top)
  foreFrame <- tkframe(top)
    foreVariable <- tclVar("1")
    foreField <- tkentry(foreFrame, width="8", textvariable=foreVariable)
  
    tkgrid(tklabel(nameFrame, fg="blue", text=gettextRcmdr("Current Model: ")), 
        tklabel(nameFrame, text=tclvalue(getRcmdr("modelName"))), sticky="w")
    tkgrid(nameFrame, sticky="w", columnspan="2")
    tkgrid(getFrame(modelsBox), columnspan="2", sticky="w")
    tkgrid(tklabel(top,text=gettextRcmdr("New Data set name")),
	entryDsname,sticky="e")
   
   tkgrid(tklabel(foreFrame, text=gettextRcmdr("Number of Forecasts: ")), foreField, sticky="w")
    tkgrid(foreFrame, sticky="w")
  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }

 

pow2 <- function(x,lam1=1) {
  if(any(x)<=0)stop("Nonnegative values; transformation not possible")
  y <- exp(lam1*log(x))
  return(y)
}

decom5 <- function(x1,fore1=0,se1=1) {
#Changes made on 2/21/06
#Remove "its" section and remove options for irre.
#x1 is the original data series;
#fore1 is the number of forecast periods
#se1 is the annual frequency
  n1 <- length(x1)
if(is.ts(x1) !=TRUE){
  x <- ts(x1,start=1,frequency=1)
}
else {
  x <- x1
}
f1 <- tsp(x)[3]
f21 <- f1
ck1 <- n1/f1
if(se1 != 1)f21 <- 1
#if(ck1 != floor(ck1))stop("Need exact values for a year")
if(fore1 < 0)stop("Forecast value must be positive")
if(fore1 > n1)stop("Forecast value must be less than series length")
#Now start the seasonal process
#This is NOT done for annual data
if(f21 == 1) {
y <- filter(x,rep(1,f1))/f1
z <- filter(y,rep(1,2))/2
xx <- as.vector(z)
z1 <- c(NA,xx[-n1])
w1 <- x/z1
#w2 <- matrix(w1,nrow=f1)
#w3 <- apply(w2,1,function(x)mean(x,na.rm=TRUE))
xw <- vector("list",f1)
for(i in 1:n1) {
  j1 <- ifelse(i%%f1==0,f1,i%%f1)
  xw[[j1]] <- c(xw[[j1]],w1[i])
}
w3 <- unlist(lapply(xw,function(x)mean(x,na.rm=TRUE)))
w4 <- sum(w3)/f1
w3 <- w3/w4
sea1 <- rep(w3,length=n1)
sea1 <- ts(sea1,start=start(x),frequency=f1)
ab <- f1 - start(x)[2] +2
sea2 <- sea1[ab:(ab+f1-1)]
dy <- x/sea1
}
else {
sea1 <- rep(1,length=n1)
sea2 <- 1
dy <- x
}
#Begin fitting the trend
t1 <- 1:n1
trend.lm <- lm(dy ~ t1)
trend.ts <- ts(trend.lm$fitted.values,start=start(x),frequency=f1)
print(trend.lm$coef)
#Obtain Final Fitted series
#2/05/2006
#Make adjustments to set up cycle and irregular as time series
yhat <- trend.ts*sea1
#We will get cyclical and irregular values
cr1 <- x/yhat
cy1 <- ts(as.vector(filter(cr1,rep(1,3))/3),start=start(x),frequency=f1)
ir1 <- cr1/cy1
#Calculate forecasts if needed
if(fore1 != 0) {
new1 <- data.frame(t1=(n1+1):(n1+fore1))
pred1 <- predict(trend.lm,newdata=new1,interval="prediction")
pred2 <- (pred1[,3] - pred1[,2])/2
xs1 <- sea1[1:fore1]
pred4 <- pred1[,1]*xs1
pred5 <- pred4 - pred2
pred6 <- pred4 + pred2
pred.df <- data.frame(pred=pred4,lower=pred5,upper=pred6)
print(pred.df)
return(data.frame(mult=pred.df$pred))
}
#x1 <- data.frame(x1,deas=dy,
#		trend=trend.ts,seas=sea1,seay=sea2,cycle=cy1,irr=ir1)
return(data.frame(mult=pred.df$pred))
}

histprice2<- function(inst1,start1="1998-01-01",quot1="Close",end1) {

  library(tseries)
 z <- get.hist.quote(instrument=inst1, start=start1,end=end1,
                     quote=quot1,compression  = "m")
y <- as.ts(aggregate(z, as.yearmon, tail, 1))
y.df <- data.frame(y=y,time=time(y))
y.df$x <- ts(y.df[,1])
tsp(y.df$x) <- tsp(y.df[,2])
names(y.df) <- c("data","time","ts")
z.df <- data.frame(ts=y.df$ts)
 return(z.df)
}
              

foreMod <- function(){
    initializeDialog(title=gettextRcmdr("Forecast Accuracy Results"))
    .numeric <- Numeric()
    xBox <- variableListBox(top, .numeric, title=gettextRcmdr("First variable (pick one)"))
    yBox <- variableListBox(top, .numeric, title=gettextRcmdr("Second variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        if (length(x) == 0 | length(y) == 0){
            errorCondition(recall=foreMod, message=gettextRcmdr("You must select two variables."))
            return()
            }
        if (x == y){
            errorCondition(recall=foreMod, message=gettextRcmdr("Variables must be different."))
            return()
            }
       closeDialog()
        .activeDataSet <- ActiveDataSet()
        doItAndPrint(paste("forerr(", .activeDataSet, "$", x, ", ", 
            .activeDataSet, "$", y,
            ")", sep=""))
        tkfocus(CommanderWindow())
        }


    OKCancelHelp(helpSubject="t.test")
   tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")    
 
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }




garchMod <- function(){
 initializeDialog(title=gettextRcmdr("GARCH Models"))
   .activeModel <- ActiveModel()
	    currentModel <- if (!is.null(.activeModel)) 
        eval(parse(text=paste("class(", .activeModel, ")[1] == 'garch'", sep="")), 
            envir=.GlobalEnv) 
        else FALSE
#    if (currentModel) {
#        currentFields <- formulaFields(eval(parse(text=.activeModel), 
#            envir=.GlobalEnv))
#        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
#        }
    UpdateModelNumber()
   modelName <- tclVar(paste("GARCHModel.", getRcmdr("modelNumber"), sep=""))
   modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=modelName)
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
          if (!is.valid.name(modelValue)){
            errorCondition(recall=garchMod, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
	         x <- getSelection(xBox)
       if (length(x) == 0){
            errorCondition(recall=garchMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
         AR1 <- tclvalue(ARLevel)
        MA1 <- tclvalue(MAVariable)
	 closeDialog() 
    
       if (is.element(modelValue, listGARCHModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                garchMod()
                return()
                }
            }
  
    command <- paste("garch(", ActiveDataSet(), "$", x,",order=c(",AR1,",",MA1,"),
		trace=F)"
		,sep="")
    
 
  logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
      doItAndPrint(paste("summary(", modelValue, ")", sep=""))
  
	
           activeModel(modelValue)
    

         tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="garch")
       ARFrame <- tkframe(top)
    ARLevel <- tclVar("0")
    ARField <- tkentry(ARFrame, width="6", textvariable=ARLevel)
    MAFrame <- tkframe(top)
    MAVariable <- tclVar("0")
    MAField <- tkentry(MAFrame, width="6", textvariable=MAVariable)
 
 tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
     tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(top, text=""), sticky="w")   
   tkgrid(tklabel(ARFrame, text=gettextRcmdr("AR Order ")), ARField, sticky="w")
    tkgrid(ARFrame, sticky="w")
    tkgrid(tklabel(MAFrame, text=gettextRcmdr("MA Order")), MAField, sticky="w")
    tkgrid(MAFrame, sticky="w")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(ARField, sticky="e")
    dialogSuffix(rows=6, columns=2)
    }



tsConv <- function(){
    initializeDialog(title=gettextRcmdr("Time Series Conversion"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
 
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=tsConv, message=gettextRcmdr("You must select a variable."))
            return()
            }

        year1 <- tclvalue(yearVariable)
        per1 <- tclvalue(perVariable)
        freq1 <- tclvalue(freqVariable)
        closeDialog()

    comm1 <- paste(ActiveDataSet(),"$",x,sep="")
    command <- paste("ts(", ActiveDataSet(), "$", x,",start=c(",year1,",",per1,"),frequency=",
		freq1,")",sep="")
    
 
 # logger(paste(comm1, " <- ", command, sep=""))
	comm2 <- paste(comm1, " <- ", command, sep="")
     #  assign(comm1, justDoIt(command), envir=.GlobalEnv)
   doItAndPrint(comm2)
      
 
      tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ts")
    yearFrame <- tkframe(top)
    yearVariable <- tclVar("2000")
    yearField <- tkentry(yearFrame, width="6", textvariable=yearVariable)
 
    perFrame <- tkframe(top)
    perVariable <- tclVar("1")
    perField <- tkentry(perFrame, width="6", textvariable=perVariable)
    freqFrame <- tkframe(top)
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
    tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(yearFrame, text=gettextRcmdr("Starting Year: ")), yearField, sticky="w")
    tkgrid(yearFrame, sticky="w")
   tkgrid(tklabel(perFrame, text=gettextRcmdr("Starting Period: ")), perField, sticky="w")
    tkgrid(perFrame, sticky="w")
     tkgrid(tklabel(freqFrame, text=gettextRcmdr("Freq: ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(yearField, sticky="e")
    dialogSuffix(rows=4, columns=2)
    }


decom6 <- function(x1,fore1=0,se1=1) {
#Additive decomposition pgs 338 - 340
#Bowerman O'Connell Koehler
#x1 is the original data series;
#fore1 is the number of forecast periods
#se1 is the annual frequency
  n1 <- length(x1)
if(is.ts(x1) !=TRUE){
  x <- ts(x1,start=1,frequency=1)
}
else {
  x <- x1
}
f1 <- tsp(x)[3]
f21 <- f1
ck1 <- n1/f1
if(se1 != 1)f21 <- 1
#if(ck1 != floor(ck1))stop("Need exact values for a year")
if(fore1 < 0)stop("Forecast value must be positive")
if(fore1 > n1)stop("Forecast value must be less than series length")
#Now start the seasonal process
#This is NOT done for annual data
if(f21 == 1) {
y <- filter(x,rep(1,f1))/f1
z <- filter(y,rep(1,2))/2
xx <- as.vector(z)
z1 <- c(NA,xx[-n1])
w1 <- x/z1
#w2 <- matrix(w1,nrow=f1)
#w3 <- apply(w2,1,function(x)mean(x,na.rm=TRUE))
xw <- vector("list",f1)
for(i in 1:n1) {
  j1 <- ifelse(i%%f1==0,f1,i%%f1)
  xw[[j1]] <- c(xw[[j1]],w1[i])
}
w3 <- unlist(lapply(xw,function(x)mean(x,na.rm=TRUE)))
w4 <- sum(w3)/f1
w3 <- w3/w4
sea1 <- rep(w3,length=n1)
sea1 <- ts(sea1,start=start(x),frequency=f1)
ab <- f1 - start(x)[2] +2
sea2 <- sea1[ab:(ab+f1-1)]
dy <- x-sea1
}
else {
sea1 <- rep(0,length=n1)
sea2 <- 0
dy <- x
}
#Begin fitting the trend
t1 <- 1:n1
trend.lm <- lm(dy ~ t1)
trend.ts <- ts(trend.lm$fitted.values,start=start(x),frequency=f1)
print(trend.lm$coef)
#Obtain Final Fitted series
#2/05/2006
#Make adjustments to set up cycle and irregular as time series
yhat <- trend.ts+sea1
#We will get cyclical and irregular values
cr1 <- x-yhat
cy1 <- ts(as.vector(filter(cr1,rep(1,3))/3),start=start(x),frequency=f1)
ir1 <- cr1-cy1
#Calculate forecasts if needed
if(fore1 != 0) {
new1 <- data.frame(t1=(n1+1):(n1+fore1))
pred1 <- predict(trend.lm,newdata=new1,interval="prediction")
pred2 <- (pred1[,3] - pred1[,2])/2
xs1 <- sea1[1:fore1]
pred4 <- pred1[,1]+xs1
pred5 <- pred4 - pred2
pred6 <- pred4 + pred2
pred.df <- data.frame(pred=pred4,lower=pred5,upper=pred6)

print(pred.df)
return(data.frame(add=pred.df$pred))
}
#x1 <- data.frame(x1,deas=dy,
#		trend=trend.ts,seas=sea1,seay=sea2,cycle=cy1,irr=ir1)
return(data.frame(add=pred.df$pred))
}

decomaMod <- function(){
    initializeDialog(title=gettextRcmdr("Additive Decomposition"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
   
 
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=decomaMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
          dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=decomaMod, 
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
    return() 
	}   
   
   
        freq1 <- tclvalue(freqVariable)
        fore1 <- tclvalue(foreVariable)
        closeDialog()
     
 comm1 <- paste(dsnameValue,"<-decom6(", ActiveDataSet(), "$", x,
            ",se1=",freq1,",fore1=",fore1,
            ")", sep="")
   	doItAndPrint(comm1)

        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="t.test")
    freqFrame <- tkframe(top)
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
    foreFrame <- tkframe(top)
    foreVariable <- tclVar("0")
    foreField <- tkentry(foreFrame, width="8", textvariable=foreVariable)
    tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(foreFrame, text=gettextRcmdr("Number of Forecasts: ")), foreField, sticky="w")
    tkgrid(foreFrame, sticky="w")
    tkgrid(tklabel(freqFrame, text=gettextRcmdr("Freq: ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
   tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  
    tkgrid.configure(entryDsname, sticky="w")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(freqField, sticky="e")
    dialogSuffix(rows=4, columns=2)
    }


hwseasMod <- function(){
  initializeDialog(title=gettextRcmdr("Seasonal Holt Winters Filter"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
     if (length(x) == 0){
            errorCondition(recall=bcMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("HoltWinters(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="HoltWinters")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


hwnonseasMod <- function(){
  initializeDialog(title=gettextRcmdr("Non-Seasonal Holt Winters Filter"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
     if (length(x) == 0){
            errorCondition(recall=bcMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("HoltWinters(", ActiveDataSet(), "$", x,",gamma=0)",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="HoltWinters")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


hwexpoMod <- function(){
  initializeDialog(title=gettextRcmdr("Exponential Smoothing"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
     if (length(x) == 0){
            errorCondition(recall=bcMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("HoltWinters(", ActiveDataSet(), "$", x,",gamma=0,beta=0)",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="HoltWinters")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }

HoltWintersMod <- function(){
 initializeDialog(title=gettextRcmdr("Seasonal HoltWinters Models"))
   .activeModel <- ActiveModel()
	print(.activeModel)
    currentModel <- if (!is.null(.activeModel)) 
        eval(parse(text=paste("class(", .activeModel, ")[1] == 'HoltWinters'", sep="")), 
            envir=.GlobalEnv) 
        else FALSE
   UpdateModelNumber()
   modelName <- tclVar(paste("HoltWintersModel.", getRcmdr("modelNumber"), sep=""))
   modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=modelName)
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
          if (!is.valid.name(modelValue)){
            errorCondition(recall=HoltWintersMod, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
	         x <- getSelection(xBox)
       if (length(x) == 0){
            errorCondition(recall=HoltWintersMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
    alternative <- as.character(tclvalue(alternativeVariable))
 	 closeDialog() 
    
       if (is.element(modelValue, listHoltWintersModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                HoltWintersMod
                return()
                }
            }
         
  pd1 <- paste("tsp(",ActiveDataSet(),"$",x,")[3]",sep="")
	
	assign("pd2",justDoIt(pd1),envir=.GlobalEnv)
		sea2 <- as.character(ifelse(alternative==1,"add","mult"))
	
   command <- paste("HoltWinters(", ActiveDataSet(), "$", x,",
	   seasonal='", alternative, "')", sep="")
       	
			
 	   
 
  logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(modelValue)
           activeModel(modelValue)
    

         tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="HoltWinters")
   radioButtons(top, name="alternative", buttons=c("yes","no"), values=c("add","mult"),
        labels=gettextRcmdr(c("Additive", "Multiplicative")), 
        title=gettextRcmdr("Seasonal Coefficients"))
      rightFrame <- tkframe(top)
  
 

 tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
     tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(rightFrame, text=""), sticky="w")   


   tkgrid(alternativeFrame, rightFrame, sticky="nw")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=6, columns=2)
    }

listHoltWintersModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects, 
        function(.x) "HoltWinters" == (class(eval(parse(text=.x), envir=envir))[1]))]
    }

HoltWintersNonMod <- function(){
 initializeDialog(title=gettextRcmdr("HoltWinters NonSeasonal Models"))
   .activeModel <- ActiveModel()
	print(.activeModel)
    currentModel <- if (!is.null(.activeModel)) 
        eval(parse(text=paste("class(", .activeModel, ")[1] == 'HoltWinters'", sep="")), 
            envir=.GlobalEnv) 
        else FALSE
   UpdateModelNumber()
   modelName <- tclVar(paste("HoltWintersModel.", getRcmdr("modelNumber"), sep=""))
   modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=modelName)
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
          if (!is.valid.name(modelValue)){
            errorCondition(recall=HoltWintersMod, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
	         x <- getSelection(xBox)
       if (length(x) == 0){
            errorCondition(recall=HoltWintersMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
  	 closeDialog() 
    
       if (is.element(modelValue, listHoltWintersModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                HoltWintersMod
                return()
                }
            }
         
 	
   command <- paste("HoltWinters(", ActiveDataSet(), "$", x,",
	   gamma=0)", sep="")
       	
			
 	   
 
  logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(modelValue)
           activeModel(modelValue)
    

         tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="HoltWinters")
     rightFrame <- tkframe(top)
  
 

 tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
     tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(rightFrame, text=""), sticky="w")   


  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=6, columns=2)
    }


HoltWintersExpoMod <- function(){
 initializeDialog(title=gettextRcmdr("Exponential Smoothing Models"))
   .activeModel <- ActiveModel()
	print(.activeModel)
    currentModel <- if (!is.null(.activeModel)) 
        eval(parse(text=paste("class(", .activeModel, ")[1] == 'HoltWinters'", sep="")), 
            envir=.GlobalEnv) 
        else FALSE
   UpdateModelNumber()
   modelName <- tclVar(paste("HoltWintersModel.", getRcmdr("modelNumber"), sep=""))
   modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=modelName)
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
          if (!is.valid.name(modelValue)){
            errorCondition(recall=HoltWintersMod, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
	         x <- getSelection(xBox)
       if (length(x) == 0){
            errorCondition(recall=HoltWintersMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
  	 closeDialog() 
    
       if (is.element(modelValue, listHoltWintersModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                HoltWintersMod
                return()
                }
            }
         
 	
   command <- paste("HoltWinters(", ActiveDataSet(), "$", x,",
	   gamma=0,beta=0)", sep="")
       	
			
 	   
 
  logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(modelValue)
           activeModel(modelValue)
    

         tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="HoltWinters")
     rightFrame <- tkframe(top)
  
 

 tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
     tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(rightFrame, text=""), sticky="w")   


  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=6, columns=2)
    }


predar3 <- function(x,fore1=1) {

z <- forecast.Arima(x,h=fore1)

zlow <- z$lower
zup <- z$upper
zz <- data.frame(pred=z$mean,lower=zlow,upper=zup) 
plot(z)
return(as.data.frame(zz$pred))
}


predhw <- function(x,fore1=1) {

z <- predict(x,prediction.interval=TRUE,n.ahead=fore1)



print(z)
return(as.data.frame(z[,1]))
}



predAllModel <- function(){
    models <- listAllModels()
    .activeModel <- ActiveModel()
    if ((length(models) == 1) && !is.null(.activeModel)) {
        Message(message=gettextRcmdr("There is only one model in memory."),
                type="warning")
        tkfocus(CommanderWindow())
        return()
        }
    if (length(models) == 0){
        Message(message=gettextRcmdr("There are no models from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
    initializeDialog(title=gettextRcmdr("Forecasting"))
    .activeDataSet <- ActiveDataSet()
    initial <- if (is.null(.activeModel)) NULL else which(.activeModel == models) - 1
    modelsBox <- variableListBox(top, models, title=gettextRcmdr("Models (pick one)"), 
        initialSelection=initial)
 dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)

    onOK <- function(){
        model <- getSelection(modelsBox)
        closeDialog()
        if (length(model) == 0) {
            tkfocus(CommanderWindow())
            return()
	}

        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=editTSframe, 
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
	return()
}

     activeModel(model)
     fore1 <- tclvalue(foreVariable)
        closeDialog()
 if(  eval(parse(text=paste("class(", model, ")[1] == 'HoltWinters'", sep="")), 
            envir=.GlobalEnv)) {
     
	   comm1 <- paste(dsnameValue, " <- predhw(", model,",fore1=",fore1,
            ")", sep="")
	doItAndPrint(comm1)




	}
	else {
	   comm1 <- paste(dsnameValue, " <- predar3(", model,",fore1=",fore1,
            ")", sep="")
	doItAndPrint(comm1)

  	}
        tkdestroy(top)
      tkfocus(CommanderWindow())
        }
    OKCancelHelp()
    nameFrame <- tkframe(top)
  foreFrame <- tkframe(top)
    foreVariable <- tclVar("1")
    foreField <- tkentry(foreFrame, width="8", textvariable=foreVariable)
  
    tkgrid(tklabel(nameFrame, fg="blue", text=gettextRcmdr("Current Model: ")), 
        tklabel(nameFrame, text=tclvalue(getRcmdr("modelName"))), sticky="w")
    tkgrid(nameFrame, sticky="w", columnspan="2")
    tkgrid(getFrame(modelsBox), columnspan="2", sticky="w")
   tkgrid(tklabel(foreFrame, text=gettextRcmdr("Number of Forecasts: ")), foreField, sticky="w")
    tkgrid(foreFrame, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  
    tkgrid.configure(entryDsname, sticky="w")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }


fixa <- function(x)
{
tsa <- tsp(x[,1])
y <- edit(as.matrix(x))
y <- data.frame(y)
for(i in 1:ncol(y)) {
y[,i] <- ts(y[,i],start=tsa[1],frequency=tsa[3])
}
return(y)
}




editTSframe <- function() {
    initializeDialog(title=gettextRcmdr("Edit TS Data"))
    dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
    onOK <- function(){
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=editTSframe, 
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
  	command <- paste('fixa(',dsnameValue,')',sep="")
	assign(dsnameValue, justDoIt(command), envir=.GlobalEnv)

        logger(paste(dsnameValue, "<-", command))
      if (eval(parse(text=paste("nrow(", dsnameValue, ")"))) == 0){
            errorCondition(recall=editTSframe, message=gettextRcmdr("empty data set."))
            return()
            }
        activeDataSet(dsnameValue)
        closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="edit")
  
  rightFrame <- tkframe(top)
    tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
  
    tkgrid.configure(entryDsname, sticky="w")
    dialogSuffix(rows=10, columns=2, focus=entryDsname)
    }


bulkfit <- function(x) {
w <- matrix(0,nrow=27,ncol=4)
ii <- 0
for(i in 0:2) {
	for(k in 0:2) {
	for(j in 0:2) {
		ii <- ii + 1
		fit <- try(arima(x,order=c(i,k,j)),silent=TRUE)
	
			if(inherits(fit,"try-error")) {
				w[ii,4] <- 99999 	
				}
			else {
			w[ii,4] <- fit$aic
			w[ii,1] <- i
			w[ii,2] <- k	
			w[ii,3] <- j
		
		}
		}
		}
	}
	
	dimnames(w) <- list(NULL,c("ar","d","ma","AIC"))
	xxx <- which(w[,4]==min(w[,4]))[1]
return(list(res=w,min=w[xxx,]))

}




decomModcom <- function(){
    initializeDialog(title=gettextRcmdr("Decomposition: Both Versions"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
 dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
 
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=decomModcom, message=gettextRcmdr("You must select a variable."))
            return()
            }
           dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=editTSframe, 
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
     return()
            }
        freq1 <- tclvalue(freqVariable)
        fore1 <- tclvalue(foreVariable)
        closeDialog()
	command <- print("Multiplicative")
	logger(command)
  comm1 <- paste(".One <-decom5(", ActiveDataSet(), "$", x,
            ",se1=",freq1,",fore1=",fore1,
            ")", sep="")
   	doItAndPrint(comm1)
   
     	command <- print("Additive")
	logger(command)
   comm1 <- paste(".Two <-decom6(", ActiveDataSet(), "$", x,
            ",se1=",freq1,",fore1=",fore1,
            ")", sep="")
   	doItAndPrint(comm1)
     
	comm2 <- paste(dsnameValue,"<-cbind(.One,.Two)",sep="")
	doItAndPrint(comm2)
	rm(.One,.Two)
   tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="t.test")
    freqFrame <- tkframe(top)
    freqVariable <- tclVar("1")
    freqField <- tkentry(freqFrame, width="6", textvariable=freqVariable)
    foreFrame <- tkframe(top)
    foreVariable <- tclVar("0")
    foreField <- tkentry(foreFrame, width="8", textvariable=foreVariable)
    tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(foreFrame, text=gettextRcmdr("Number of Forecasts: ")), foreField, sticky="w")
    tkgrid(foreFrame, sticky="w")
    tkgrid(tklabel(freqFrame, text=gettextRcmdr("Freq: ")), freqField, sticky="w")
    tkgrid(freqFrame, sticky="w")
  tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  
    tkgrid.configure(entryDsname, sticky="w")
  
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(freqField, sticky="e")
    dialogSuffix(rows=4, columns=2)
    }






runbulk <- function(){
  initializeDialog(title=gettextRcmdr("Run Multiple ARIMA models"))
   xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
     if (length(x) == 0){
            errorCondition(recall=runbulk, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("bulkfit(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="arima")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }


runbulkg <- function(){
  initializeDialog(title=gettextRcmdr("Run Multiple GARCH models"))
   xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
     if (length(x) == 0){
            errorCondition(recall=runbulkg, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("bulkfitg(", ActiveDataSet(), "$", x,")",sep=""))     
		
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="garch")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }



bulkfitg <- function(x) {
w <- matrix(0,nrow=8,ncol=3)
oldopts <- options()
options(show.error.messages=FALSE)
ii <- 0
for(i in 0:2) {
	for(j in 0:2) {
		if(i==0 & j==0)next()
		ii <- ii + 1
		cat(i,j,ii,"\n")
		fit <- try(garch(x,order=c(i,j)),silent=TRUE)
		if(inherits(fit,"try-error")) {
				w[ii,3] <- 99999 	
				}
			else {
				if(any(is.nan(fit[[5]])))w[ii,3] <- 99999
				else {
				w[ii,3] <- AIC(fit)
				w[ii,1] <- i
				w[ii,2] <- j
			}
		
		}
		}
	}
	
	dimnames(w) <- list(NULL,c("ar","ma","AIC"))
	xxx <- which(w[,3]==min(w[,3]))[1]
	options(oldopts)
return(list(res=w,min=w[xxx,]))

}



mse1 <- function(x,y) {
     n1 <- length(x)
     z <- sum((x-y)^2)/n1
     return(z)
}

mad1 <- function(x,y) {
     n1 <- length(x)
     z <- sum(abs(x-y))/n1
     return(z)
}




 fore1 <- function(x,w,y) {
        z <- which(names(y) == as.character(x))
	u <- y[,-z,drop=FALSE]
	print("MSE")
	a <- apply(u,2,mse1,y=w)
	print(a)
	print("MAD")
	b <- apply(u,2,mad1,y=w)
	print(b)
	return(data.frame(mse=a,mad=b))
	
	}







forebMod <- function(){
    initializeDialog(title=gettextRcmdr("Multiple Forecast Accuracy Results"))
    .numeric <- Numeric()
    xBox <- variableListBox(top, .numeric, title=gettextRcmdr("First variable (pick one)"))
 dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
   
     onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=forebMod, message=gettextRcmdr("You must select a variable."))
		}
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=editTSframe, 
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
     return()
            }
       closeDialog()
        .activeDataSet <- ActiveDataSet()
	x1 <- as.character(x)
 #       doItAndPrint(paste("forerr(", .activeDataSet, "$", x, ", ", 
 #           .activeDataSet, "$", y,
 #           ")", sep=""))

     comm1 <- paste(dsnameValue, "<- fore1('",x1,"',", ActiveDataSet(), "$", x, ", ",ActiveDataSet(), 
            ")", sep="")
      doItAndPrint(comm1)


        tkfocus(CommanderWindow())
        }


    OKCancelHelp(helpSubject="t.test")
   tkgrid(getFrame(xBox), sticky="nw")    
   tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  
    tkgrid.configure(entryDsname, sticky="w")
  
 
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }



subMod <- function(){
  initializeDialog(title=gettextRcmdr("Plot Plus Subset"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
	print(x)
        if (length(x) == 0){
            errorCondition(recall=subMod, message=gettextRcmdr("You must select a variable."))
            return()
            }
        closeDialog()
    doItAndPrint(paste("ts.plot(", ActiveDataSet(), "$", x,")",sep=""))     
	 comm1 <- paste("plotsub(", ActiveDataSet(), "$", x,
              ")", sep="")
   	doItAndPrint(comm1)
     
	
        tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ts.plot")
   tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
     dialogSuffix(rows=4, columns=2)
    }

plotsub <- function(x) {
plot(x)
z <- locator(2)
plot(window(x,start=min(z$x),end=max(z$x)))
}

aggConv <- function(){
    initializeDialog(title=gettextRcmdr("Low Frequency Aggregation"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  dsname <- tclVar(gettextRcmdr("Dataset"))
   
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
 
 
    onOK <- function(){
        x <- getSelection(xBox)
        if (length(x) == 0){
            errorCondition(recall=aggConv, message=gettextRcmdr("You must select a variable."))
            return()
            }
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "") {
            errorCondition(recall=aggConv, 
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
  alternative <- as.character(tclvalue(alternativeVariable))
      
 alternativea <- as.character(tclvalue(alternativeaVariable))
  

     alt1 <- paste("",alternative,"",sep="")
          alt2 <- paste("",alternativea,"",sep="")
        closeDialog()

  	comma <- paste(ActiveDataSet(),"$",x,sep="")
	comm3 <- paste("ep <- endpoints(",comma,",on='",alt2,"')",sep="")
     	doItAndPrint(comm3)
 
 # logger(paste(comm1, " <- ", command, sep=""))
	
	comm2 <- paste(dsnameValue, " <- data.frame(ts=period.apply(as.xts(", comma,"),ep,",alt1,"))", sep="")
     
      comm4 <- paste(".One <- period.apply(as.xts(", comma,"),ep,",alt1,")", sep="")
      doItAndPrint(comm4)
      comm5 <- paste(dsnameValue, "<- data.frame(x=1:(length(ep)-1))",sep="")
     doItAndPrint(comm5)
      comm6 <- paste(dsnameValue,"$x <- zoo(.One)",sep="")
     doItAndPrint(comm6)
    
     #  assign(comm1, justDoIt(command), envir=.GlobalEnv)
     
 
      tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="period.apply")
    radioButtons(top, name="alternative", buttons=c("Sum", "Mean","Min","Max"), 
	values=c("sum", "mean","min","max"),
        labels=gettextRcmdr(c("Sum", "Mean",
	"Minimum","Maximum")), 
        title=gettextRcmdr("Aggregation Options"))

 rightFrame <- tkframe(top)
   tkgrid(alternativeFrame, rightFrame, sticky="nw")
   
     tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
 
    radioButtons(top, name="alternativea", buttons=c("Months", "Quarters","Years"), 
	values=c("months", "quarters","years"),
        labels=gettextRcmdr(c("Months",
	"Quarters","Years")), 
        title=gettextRcmdr("Time  Options"))

 rightFrame1 <- tkframe(top)
   tkgrid(alternativeaFrame, rightFrame1, sticky="nw")



    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }

ArimaModIn <- function(){
 initializeDialog(title=gettextRcmdr("Interactive ARIMA Models"))
 
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  
    onOK <- function(){
	         x <- getSelection(xBox)
       if (length(x) == 0){
            errorCondition(recall=ArimaModIn, message=gettextRcmdr("You must select a variable."))
            return()
            }
         AR1 <- tclvalue(ARLevel)
        MA1 <- tclvalue(MAVariable)
	DIF1 <- tclvalue(DIFVariable)
    
	SAR1 <- tclvalue(SARLevel)
        SMA1 <- tclvalue(SMAVariable)
	SDIF1 <- tclvalue(SDIFVariable)
 
   alternative <- as.character(tclvalue(alternativeVariable))
 	 closeDialog() 
    
         
      pd1 <- paste("tsp(",ActiveDataSet(),"$",x,")[3]",sep="")

	assign("pd2",justDoIt(pd1),envir=.GlobalEnv)
    command <- paste("fit <- arima(", ActiveDataSet(), "$", x,",order=c(",AR1,",",DIF1,",",MA1,"),
		include.mean=",alternative,",seasonal=list(order=c(",SAR1,",",SDIF1,",",SMA1
		,"),period=",pd2,"))",sep="")
    
 
  
#doItAndPrint(paste("ts.plot(", ActiveDataSet(), "$", x,")",sep=""))     
doItAndPrint(paste("fit1 <- predict(fit,n.ahead=7)",sep=""))        

         tkdestroy(top)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="arima")
   radioButtons(top, name="alternative", buttons=c("yes","no"), values=c(TRUE,FALSE),
        labels=gettextRcmdr(c("Mean Included","Mean Not Included")), 
        title=gettextRcmdr("Mean"))
      rightFrame <- tkframe(top)
    ARFrame <- tkframe(rightFrame)
    ARLevel <- tclVar("0")
    ARField <- tkentry(ARFrame, width="6", textvariable=ARLevel)
    MAFrame <- tkframe(rightFrame)
    MAVariable <- tclVar("0")
    MAField <- tkentry(MAFrame, width="6", textvariable=MAVariable)
    DIFFrame <- tkframe(rightFrame)
    DIFVariable <- tclVar("0")
    DIFField <- tkentry(DIFFrame,width="6",textvariable=DIFVariable)
 
    SARFrame <- tkframe(rightFrame)
    SARLevel <- tclVar("0")
    SARField <- tkentry(SARFrame, width="6", textvariable=SARLevel)
    SMAFrame <- tkframe(rightFrame)
    SMAVariable <- tclVar("0")
    SMAField <- tkentry(SMAFrame, width="6", textvariable=SMAVariable)
    SDIFFrame <- tkframe(rightFrame)
    SDIFVariable <- tclVar("0")
    SDIFField <- tkentry(SDIFFrame,width="6",textvariable=SDIFVariable)


    tkgrid(getFrame(xBox), sticky="nw") 
    tkgrid(tklabel(rightFrame, text=""), sticky="w")   
   tkgrid(tklabel(ARFrame, text=gettextRcmdr("Reg. AR Order ")), ARField, sticky="w")
    tkgrid(ARFrame, sticky="w")
    tkgrid(tklabel(DIFFrame,text=gettextRcmdr("Reg. Diff Order ")),DIFField, sticky="w")
    tkgrid(DIFFrame,sticky="w")
    tkgrid(tklabel(MAFrame, text=gettextRcmdr("Reg. MA Order")), MAField, sticky="w")
    tkgrid(MAFrame, sticky="w")
  
   tkgrid(tklabel(SARFrame, text=gettextRcmdr("Sea. AR Order ")), SARField, sticky="w")
    tkgrid(SARFrame, sticky="w")
    tkgrid(tklabel(SDIFFrame,text=gettextRcmdr("Sea. Diff Order ")),SDIFField, sticky="w")
    tkgrid(SDIFFrame,sticky="w")
    tkgrid(tklabel(SMAFrame, text=gettextRcmdr("Sea. MA Order")), SMAField, sticky="w")
    tkgrid(SMAFrame, sticky="w")
 

   tkgrid(alternativeFrame, rightFrame, sticky="nw")
 
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    tkgrid.configure(ARField, sticky="e")
    dialogSuffix(rows=6, columns=2)
    }


foreAllModel <- function(){
    models <- listAllModels()
    .activeModel <- ActiveModel()
    if ((length(models) == 1) && !is.null(.activeModel)) {
        Message(message=gettextRcmdr("There is only one model in memory."),
                type="warning")
        tkfocus(CommanderWindow())
        return()
        }
    if (length(models) == 0){
        Message(message=gettextRcmdr("There are no models from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
    initializeDialog(title=gettextRcmdr("Plotting Forecasts"))
    .activeDataSet <- ActiveDataSet()
    initial <- if (is.null(.activeModel)) NULL else which(.activeModel == models) - 1
    modelsBox <- variableListBox(top, models, title=gettextRcmdr("Models (pick one)"), 
        initialSelection=initial)
 
    onOK <- function(){
        model <- getSelection(modelsBox)
        closeDialog()
        if (length(model) == 0) {
            tkfocus(CommanderWindow())
            return()
	}


           
     activeModel(model)
	

	   comm1 <- paste("plot.forecast(", model,")", sep="")
	doItAndPrint(comm1)

          tkdestroy(top)
      tkfocus(CommanderWindow())
	}
    OKCancelHelp()
    nameFrame <- tkframe(top)
 
    tkgrid(tklabel(nameFrame, fg="blue", text=gettextRcmdr("Current Model: ")), 
        tklabel(nameFrame, text=tclvalue(getRcmdr("modelName"))), sticky="w")
    tkgrid(nameFrame, sticky="w", columnspan="2")
    tkgrid(getFrame(modelsBox), columnspan="2", sticky="w")

    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }


