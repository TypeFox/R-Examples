# Distributions menu dialogs for plots

# last modified 2014-07-30 by J. Fox

#   many distributions added (and some other changes) by Miroslav Ristic  (20 July 06)
# modified by Miroslav M. Ristic (15 January 11)

normalDistributionPlot <- function() {distributionPlot("normal")}
tDistributionPlot <- function() {distributionPlot("t")}
chisqDistributionPlot <- function() {distributionPlot("chisq")}
FDistributionPlot <- function() {distributionPlot("F")}
exponentialDistributionPlot <- function() {distributionPlot("exponential")}
uniformDistributionPlot <- function() {distributionPlot("uniform")}
betaDistributionPlot <- function(){distributionPlot("beta")}
CauchyDistributionPlot <- function(){distributionPlot("Cauchy")}
logisticDistributionPlot <- function(){distributionPlot("logistic")}
lognormalDistributionPlot <- function(){distributionPlot("lognormal")}
gammaDistributionPlot <- function(){distributionPlot("gamma")}
WeibullDistributionPlot <- function(){distributionPlot("Weibull")}
GumbelDistributionPlot <- function(){distributionPlot("Gumbel")}
binomialDistributionPlot <- function(){discreteDistributionPlot("binomial")}
PoissonDistributionPlot <- function(){discreteDistributionPlot("Poisson")}

# the following functions were contributed by G. Jay Kerns, Andy Chang, and  Theophilius Boye
#  last modified 26 July 06 by J. Fox

geomDistributionPlot  <- function(){discreteDistributionPlot("geom")}
hyperDistributionPlot  <- function(){discreteDistributionPlot("hyper")}
negbinomialDistributionPlot  <- function(){discreteDistributionPlot("negbinomial")}

distributionPlot <- function(nameVar){
    fVar<-get(paste(nameVar,"Distribution",sep=""))
    nnVar<-length(fVar$params)
    dialogName <- paste(nameVar,"DistributionPlot", sep="")
    defaults <- list(initialValues=fVar$initialValues, type="Density")
    initial <- getDialog(dialogName, defaults=defaults)
    initializeDialog(title=gettextRcmdr(paste(fVar$titleName,"Distribution",sep=" ")))
    entriesFrame <- tkframe(top)
    paramsVar<-paste(fVar$params,"Var",sep="")
    paramsEntry<-paste(fVar$params,"Entry",sep="")
    for (i in 1:nnVar) {
        eval(parse(text=paste(paramsVar[i],"<-tclVar('",initial$initialValues[i],"')",sep="")))
        eval(parse(text=paste(paramsEntry[i],"<-ttkentry(entriesFrame, width='6', textvariable=",paramsVar[i],")",sep="")))
    }
    functionVar <- tclVar(initial$type)
    buttonFrame <- tkframe(top)
    densityButton <- ttkradiobutton(buttonFrame, variable=functionVar, value="Density")
    distributionButton <- ttkradiobutton(buttonFrame, variable=functionVar, value="Cumulative Probability")
    onOK <- function(){
        nameVarF<-get(paste(nameVar,"DistributionPlot",sep=""),mode="function")
        closeDialog()
        warn <- options(warn=-1)
        vars<-numeric(nnVar)
        for (i in 1:nnVar) {
            vars[i]<-as.numeric(tclvalue(get(paramsVar[i])))
        }
        if (length(fVar$paramsRound)>0) {
            for (j in fVar$paramsRound) {
                vars[j]<-round(vars[j])
            }
        }
        options(warn)
        for (i in 1:length(fVar$errorConds)) {
            if (eval(parse(text=fVar$errorConds[i]))) {
                errorCondition(recall=nameVarF, message=gettextRcmdr(fVar$errorTexts[i]))
                return()
            }
        }
        fun <- tclvalue(functionVar)
        fn <- if (fun == "Density") paste("d",fVar$funName,sep="") else paste("p",fVar$funName,sep="")
        dist.arg <- if (fun == "Density") "FALSE" else "TRUE"
        pasteVar<-""
        for (i in 1:nnVar) {
            pasteVar<-paste(pasteVar,", ",fVar$params[i],"=",vars[i],sep="")
        }
        mainVar<-""
        for (i in 1:nnVar) {
            mainVar<-paste(mainVar,", ",fVar$paramsLabels[i],"=",vars[i],sep="")
        }
        if (nameVar=="Gumbel") {
            min <- eval(parse(text=paste("round(log(q",fVar$funName,"(.0005",pasteVar,")),3)",sep="")))
            max <- eval(parse(text=paste("round(log(q",fVar$funName,"(.9995",pasteVar,")),3)",sep="")))
        } else {
            min <- eval(parse(text=paste("round(q",fVar$funName,"(.0005",pasteVar,"),3)",sep="")))
            max <- eval(parse(text=paste("round(q",fVar$funName,"(.9995",pasteVar,"),3)",sep="")))
        }
        if (nameVar=="Gumbel") {
            command <- paste("exp(seq(", min, ", ", max, ", length.out=1000))", sep="")
        } else {
            command <- paste("seq(", min, ", ", max, ", length.out=1000)", sep="")
        }
        command <- paste("local({\n  .x <- ", command, sep="")
        doVar<-"\n  plotDistr(.x, "
        if (nameVar=="Gumbel") {doVar<-"\n  plotDistr(log(.x), "}
        if (nameVar=="F") {mainVar<-paste(", Numerator df = ",vars[1],", Denominator df = ",vars[2],sep="")}
        command <- paste(command, "  ", doVar, fn, "(.x", pasteVar,'), cdf=', dist.arg, ', xlab="x", ylab="', fun, 
            '", main=paste("',fVar$titleName,' Distribution: ',substr(mainVar,2,nchar(mainVar)),'"))\n})', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        putDialog(dialogName, list(initialValues=vars, type=fun), resettable=FALSE)
    }
    OKCancelHelp(helpSubject=paste("d",fVar$funName,sep=""), reset=dialogName, apply=dialogName)
    for (i in 1:nnVar) {
        tkgrid(labelRcmdr(entriesFrame, text=gettextRcmdr(fVar$paramsLabels[i])), get(paramsEntry[i]), sticky="w", padx=6)
    }
    tkgrid(entriesFrame, sticky="w")
    tkgrid(buttonFrame, sticky="w")
    tkgrid(densityButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Plot density function")), sticky="w")
    tkgrid(distributionButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Plot distribution function")), sticky="w")
    tkgrid(buttonsFrame, sticky="ew")
    for (i in 1:nnVar) {
        tkgrid.configure(get(paramsEntry[i]), sticky="w")
    }
    dialogSuffix(focus=get(paramsEntry[1]))
}

discreteDistributionPlot <- function(nameVar){
    fVar<-get(paste(nameVar,"Distribution",sep=""))
    nnVar<-length(fVar$params)
    dialogName <- paste(nameVar,"DistributionPlot", sep="")
    defaults <- list(initialValues=fVar$initialValues, type="Probability")
    initial <- getDialog(dialogName, defaults=defaults)
    initializeDialog(title=gettextRcmdr(paste(fVar$titleName,"Distribution",sep=" ")))
    entriesFrame <- tkframe(top)
    paramsVar<-paste(fVar$params,"Var",sep="")
    paramsEntry<-paste(fVar$params,"Entry",sep="")
    for (i in 1:nnVar) {
        eval(parse(text=paste(paramsVar[i],"<-tclVar('",initial$initialValues[i],"')",sep="")))
        eval(parse(text=paste(paramsEntry[i],"<-ttkentry(entriesFrame, width='6', textvariable=",paramsVar[i],")",sep="")))
    }
    functionVar <- tclVar(initial$type)
    buttonFrame <- tkframe(top)
    densityButton <- ttkradiobutton(buttonFrame, variable=functionVar, value="Probability")
    distributionButton <- ttkradiobutton(buttonFrame, variable=functionVar, value="Cumulative Probability")
    onOK <- function(){
        nameVarF<-get(paste(nameVar,"DistributionPlot",sep=""),mode="function")
        closeDialog()
        warn <- options(warn=-1)
        vars<-numeric(nnVar)
        for (i in 1:nnVar) {
            vars[i]<-as.numeric(tclvalue(get(paramsVar[i])))
        }
        if (length(fVar$paramsRound)>0) {
            for (j in fVar$paramsRound) {
                vars[j]<-round(vars[j])
            }
        }
        options(warn)
        for (i in 1:length(fVar$errorConds)) {
            if (eval(parse(text=fVar$errorConds[i]))) {
                errorCondition(recall=nameVarF, message=gettextRcmdr(fVar$errorTexts[i]))
                return()
            }
        }
        fun <- tclvalue(functionVar)
        pasteVar<-""
        for (i in 1:nnVar) {
            pasteVar<-paste(pasteVar,", ",fVar$params[i],"=",vars[i],sep="")
        }
        min <- eval(parse(text=paste("q",fVar$funName,"(.0005",pasteVar,")",sep="")))
        max <- eval(parse(text=paste("q",fVar$funName,"(.9995",pasteVar,")",sep="")))
        command <- paste(min, ":", max, sep="")
        command <- paste("local({\n  .x <- ", command, sep="")
        switch(nameVar,
            "binomial" = xlabVar<-"Number of Successes",
            "Poisson" = xlabVar<-"x",
            "geom" = xlabVar<-"Number of Failures until Success",
            "hyper" = xlabVar<-"Number of White Balls in Sample",
            "negbinomial" = xlabVar <-"Number of Failures Until Target Successes"
        )
        mainVar<-""
        if (nameVar=="negbinomial") {
            mainVar<-paste(", Trials=",vars[1],", Prob=",vars[2],sep="")
        } else if (nameVar=="hyper") {
            mainVar<-paste(", m=",vars[1],", n=",vars[2],", k=",vars[3],sep="")
        } else {
            for (i in 1:nnVar) {
                mainVar<-paste(mainVar,", ", fVar$paramsLabels[i],"=",vars[i],sep="")
            }   
        }
        command <- if (fun == "Probability"){
           paste(command, "\n  plotDistr(.x, d",fVar$funName,"(.x", pasteVar,
                '), xlab="',xlabVar,'", ylab="Probability Mass", main="',fVar$titleName,
                ' Distribution: ',substr(mainVar,2,nchar(mainVar)),'", discrete=TRUE)', sep="")
        }
        else {
            paste(command, "\n  plotDistr(.x, p",fVar$funName,"(.x",
                pasteVar,'), xlab="',xlabVar,
                '",ylab="Cumulative Probability", main="',
                fVar$titleName,' Distribution: ',substr(mainVar,2,nchar(mainVar)),'", discrete=TRUE, cdf=TRUE)', 
                sep="")
        }
        command <- paste(command, "\n})", sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        putDialog(dialogName, list(initialValues=vars, type=fun), resettable=FALSE)
    }
    OKCancelHelp(helpSubject=paste("d",fVar$funName,sep=""), reset=dialogName, apply=dialogName)
    for (i in 1:nnVar) {
        tkgrid(labelRcmdr(entriesFrame, text=gettextRcmdr(fVar$paramsLabels[i])), get(paramsEntry[i]), sticky="w", padx=6)
    }
    tkgrid(densityButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Plot probability mass function")), sticky="w")
    tkgrid(distributionButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Plot distribution function")), sticky="w")
    tkgrid(entriesFrame, sticky="w")
    tkgrid(buttonFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    for (i in 1:nnVar) {
        tkgrid.configure(get(paramsEntry[i]), sticky="w")
    }
    dialogSuffix(focus=get(paramsEntry[1]))
}
