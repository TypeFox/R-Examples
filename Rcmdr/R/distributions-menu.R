# Distributions menu dialogs

# last modified 2014-10-08 by J. Fox

#   many distributions added (and some other changes) by Miroslav Ristic (20 July 06)
# Modified by Miroslav Ristic (15 January 11)

normalDistribution <- list(titleName="Normal",
		funName="norm",
		params=c("mean","sd"),
		paramsLabels=c("Mean","Standard deviation"),
		initialValues=c("0","1"),
		errorTexts=c("Mean not specified.",
				"Standard deviation must be positive."),
		errorConds=c("is.na(vars[1])",
				"is.na(vars[2]) || (vars[2] <= 0)"),
		paramsRound=c()
)

normalQuantiles <- function() {distributionQuantiles("normal")}
normalProbabilities <-function() {distributionProbabilities("normal")}

tDistribution <- list(titleName="t",
		funName="t",
		params=c("df"),
		paramsLabels=c("Degrees of freedom"),
		initialValues=c(""),
		errorTexts=c("Degrees of freedom not specified.",
				"Degrees of freedom must be positive."),
		errorConds=c("is.na(vars[1])",
				"(vars[1] <= 0)"),
		paramsRound=c()
)

tQuantiles <- function() {distributionQuantiles("t")}
tProbabilities <-function() {distributionProbabilities("t")}

chisqDistribution <- list(titleName="ChiSquared",
		funName="chisq",
		params=c("df"),
		paramsLabels=c("Degrees of freedom"),
		initialValues=c(""),
		errorTexts=c("Degrees of freedom not specified.",
				"Degrees of freedom must be positive."),
		errorConds=c("is.na(vars[1])",
				"(vars[1] <= 0)"),
		paramsRound=c()
)

chisqQuantiles <- function() {distributionQuantiles("chisq")}
chisqProbabilities <-function() {distributionProbabilities("chisq")}

FDistribution <- list(titleName="F",
		funName="f",
		params=c("df1","df2"),
		paramsLabels=c("Numerator degrees of freedom",
				"Denominator degrees of freedom"),
		initialValues=c("",""),
		errorTexts=c("Degrees of freedom not specified.",
				"Degrees of freedom must be positive."),
		errorConds=c("is.na(vars[1]) || is.na(vars[2])",
				"(vars[1] <= 0 || vars[2] <= 0)"),
		paramsRound=c()
)

FQuantiles <- function() {distributionQuantiles("F")}
FProbabilities <-function() {distributionProbabilities("F")}

exponentialDistribution <- list(titleName="Exponential",
		funName="exp",
		params=c("rate"),
		paramsLabels=c("Rate"),
		initialValues=c("1"),
		errorTexts=c("Rate must be positive."),
		errorConds=c("is.na(vars[1]) || vars[1] <= 0"),
		paramsRound=c()
)

exponentialQuantiles <- function() {distributionQuantiles("exponential")}
exponentialProbabilities <-function() {distributionProbabilities("exponential")}

uniformDistribution <- list(titleName="Uniform",
		funName="unif",
		params=c("min","max"),
		paramsLabels=c("Minimum","Maximum"),
		initialValues=c("0","1"),
		errorTexts=c("Lower limit must be less than upper limit."),
		errorConds=c("is.na(vars[1]) || is.na(vars[2]) || vars[1] >= vars[2]"),
		paramsRound=c()
)

uniformQuantiles <- function() {distributionQuantiles("uniform")}
uniformProbabilities <-function() {distributionProbabilities("uniform")}

betaDistribution <- list(titleName="Beta",
		funName="beta",
		params=c("shape1","shape2"),
		paramsLabels=c("Shape 1","Shape 2"),
		initialValues=c("",""),
		errorTexts=c("Shapes not specified.",
				"Shapes must be positive."),
		errorConds=c("is.na(vars[1]) || is.na(vars[2])",
				"vars[1] <= 0 || vars[2] <= 0"),
		paramsRound=c()
)

betaQuantiles <- function() {distributionQuantiles("beta")}
betaProbabilities <-function() {distributionProbabilities("beta")}

CauchyDistribution <- list(titleName="Cauchy",
		funName="cauchy",
		params=c("location","scale"),
		paramsLabels=c("Location","Scale"),
		initialValues=c("0","1"),
		errorTexts=c("Location not specified.",
				"Scale must be positive."),
		errorConds=c("is.na(vars[1])",
				"is.na(vars[2]) || vars[2] <= 0"),
		paramsRound=c()
)

CauchyQuantiles <- function() {distributionQuantiles("Cauchy")}
CauchyProbabilities <-function() {distributionProbabilities("Cauchy")}

logisticDistribution <- list(titleName="Logistic",
		funName="logis",
		params=c("location","scale"),
		paramsLabels=c("Location","Scale"),
		initialValues=c("0","1"),
		errorTexts=c("Location not specified.",
				"Scale must be positive."),
		errorConds=c("is.na(vars[1])",
				"is.na(vars[2]) || vars[2] <= 0"),
		paramsRound=c()
)

logisticQuantiles <- function() {distributionQuantiles("logistic")}
logisticProbabilities <-function() {distributionProbabilities("logistic")}

lognormalDistribution <- list(titleName="Lognormal",
		funName="lnorm",
		params=c("meanlog","sdlog"),
		paramsLabels=c("Mean (log scale)","Standard deviation (log scale)"),
		initialValues=c("0","1"),
		errorTexts=c("Mean not specified.",
				"Standard deviation must be positive."),
		errorConds=c("is.na(vars[1])",
				"is.na(vars[2]) || vars[2] <= 0"),
		paramsRound=c()
)

lognormalQuantiles <- function() {distributionQuantiles("lognormal")}
lognormalProbabilities <-function() {distributionProbabilities("lognormal")}

gammaDistribution <- list(titleName="Gamma",
		funName="gamma",
		params=c("shape","scale"),
		paramsLabels=c("Shape","Scale (inverse rate)"),
		initialValues=c("","1"),
		errorTexts=c("Shape not specified.",
				"Shape must be positive.",
				"Scale must be positive."),
		errorConds=c("is.na(vars[1])",
				"vars[1]<=0",
				"is.na(vars[2]) || vars[2] <= 0"),
		paramsRound=c()
)

gammaQuantiles <- function() {distributionQuantiles("gamma")}
gammaProbabilities <-function() {distributionProbabilities("gamma")}

WeibullDistribution <- list(titleName="Weibull",
		funName="weibull",
		params=c("shape","scale"),
		paramsLabels=c("Shape","Scale"),
		initialValues=c("","1"),
		errorTexts=c("Shape not specified.",
				"Shape must be positive.",
				"Scale must be positive."),
		errorConds=c("is.na(vars[1])",
				"vars[1]<=0",
				"is.na(vars[2]) || vars[2] <= 0"),
		paramsRound=c()
)

WeibullQuantiles <- function() {distributionQuantiles("Weibull")}
WeibullProbabilities <-function() {distributionProbabilities("Weibull")}

GumbelDistribution <- list(titleName="Gumbel",
		funName="weibull",
		params=c("shape","scale"),
		paramsLabels=c("Shape (log scale)","Scale (log scale)"),
		initialValues=c("","1"),
		errorTexts=c("Shape not specified.",
				"Shape must be positive.",
				"Scale must be positive."),
		errorConds=c("is.na(vars[1])",
				"vars[1]<=0",
				"is.na(vars[2]) || vars[2] <= 0"),
		paramsRound=c()
)

GumbelQuantiles <- function() {distributionQuantiles("Gumbel")}
GumbelProbabilities <-function() {distributionProbabilities("Gumbel")}

binomialDistribution <- list(titleName="Binomial",
		funName="binom",
		params=c("size","prob"),
		paramsLabels=c("Binomial trials","Probability of success"),
		initialValues=c("","0.5"),
		errorTexts=c("Binomial trials not specified.",
				"Binomial trials must be positive.",
				"Probability of success not specified.",
				"Probability of success must be between 0 and 1."),
		errorConds=c("is.na(vars[1])",
				"vars[1]<=0",
				"is.na(vars[2])",
				"vars[2]<0 || vars[2]>1"),
		paramsRound=c(1)
)

binomialQuantiles <- function() {distributionQuantiles("binomial")}
binomialProbabilities <-function() {distributionProbabilities("binomial")}
binomialMass <- function() {distributionMass("binomial")}

# the following functions were contributed by G. Jay Kerns, Andy Chang, and  Theophilius Boye
#  modified by J. Fox
#  modified by Miroslav Ristic (15 January 2011)

PoissonDistribution <- list(titleName="Poisson",
		funName="pois",
		params=c("lambda"),
		paramsLabels=c("Mean"),
		initialValues=c("1"),
		errorTexts=c("Poisson mean not specified.",
				"Poisson mean cannot be negative."),
		errorConds=c("is.na(vars[1])",
				"vars[1]<0"),
		paramsRound=c()
)

PoissonQuantiles <- function() {distributionQuantiles("Poisson")}
PoissonProbabilities <-function() {distributionProbabilities("Poisson")}
PoissonMass <-function() {distributionMass("Poisson")}

geomDistribution <- list(titleName="Geometric",
		funName="geom",
		params=c("prob"),
		paramsLabels=c("Probability of success"),
		initialValues=c("0.5"),
		errorTexts=c("Probability of success not specified.",
				"Probability of success must be between 0 and 1."),
		errorConds=c("is.na(vars[1])",
				"vars[1] < 0 || vars[1] > 1"),
		paramsRound=c()
)

geomQuantiles <- function() {distributionQuantiles("geom")}
geomProbabilities <-function() {distributionProbabilities("geom")}
geomMass <-function() {distributionMass("geom")}

hyperDistribution <- list(titleName="Hypergeometric",
		funName="hyper",
		params=c("m","n","k"),
		paramsLabels=c("m (number of white balls in the urn)",
				"n (number of black balls in the urn)",
				"k (number of balls drawn from the urn)"),
		initialValues=c("1","1","1"),
		errorTexts=c("The m parameter was not specified.",
				"The m parameter cannot be negative.",
				"The n parameter was not specified.",
				"The n parameter cannot be negative.",
				"The k parameter was not specified.",
				"The k parameter cannot be negative.",
				"The k parameter cannot be greater than m + n."),
		errorConds=c("is.na(vars[1])",
				"vars[1] < 0",
				"is.na(vars[2])",
				"vars[2] <0",
				"is.na(vars[3])",
				"vars[3] <0",
				"vars[3]>(vars[1]+vars[2])"),
		paramsRound=c(1,2,3)
)

hyperQuantiles <- function() {distributionQuantiles("hyper")}
hyperProbabilities <-function() {distributionProbabilities("hyper")}
hyperMass <-function() {distributionMass("hyper")}

negbinomialDistribution <- list(titleName="NegativeBinomial",
		funName="nbinom",
		params=c("size","prob"),
		paramsLabels=c("Target number of successes",
				"Probability of success"),
		initialValues=c("1","0.5"),
		errorTexts=c("Target number of successes not specified.",
				"Target number of successes cannot be negative.",
				"Probability of success not specified.",
				"Probability of success must be between 0 and 1."),
		errorConds=c("is.na(vars[1])",
				"vars[1]<0",
				"is.na(vars[2])",
				"vars[2] < 0 || vars[2] >1 "),
		paramsRound=c(1)
)

negbinomialQuantiles <- function() {distributionQuantiles("negbinomial")}
negbinomialProbabilities <-function() {distributionProbabilities("negbinomial")}
negbinomialMass <-function() {distributionMass("negbinomial")}

distributionQuantiles <- function(nameVar){
	fVar<-get(paste(nameVar,"Distribution",sep=""))
	nnVar<-length(fVar$params)
	dialogName <- paste(nameVar,"Quantiles", sep="")
	defaults <- list(initialValues=fVar$initialValues, tail="lower", quantiles="")
	initial <- getDialog(dialogName, defaults=defaults)
	initializeDialog(title=gettextRcmdr(paste(fVar$titleName,"Quantiles",sep=" ")))
    entryFrame <- tkframe(top)
	quantilesVar <- tclVar(initial$quantiles)
	quantilesEntry <- ttkentry(entryFrame, width="30", textvariable=quantilesVar)
	paramsVar<-paste(fVar$params,"Var",sep="")
	paramsEntry<-paste(fVar$params,"Entry",sep="")
	for (i in 1:nnVar) {
		eval(parse(text=paste(paramsVar[i],"<-tclVar('",initial$initialValues[i],"')",sep="")))
		eval(parse(text=paste(paramsEntry[i],"<-ttkentry(entryFrame, width='6', textvariable=",paramsVar[i],")",sep="")))
	}
	tailVar <- tclVar(initial$tail)
    buttonFrame <- tkframe(top)
	lowerTailButton <- ttkradiobutton(buttonFrame, variable=tailVar, value="lower")
	upperTailButton <- ttkradiobutton(buttonFrame, variable=tailVar, value="upper")
	onOK <- function(){
		nameVarF<-get(paste(nameVar,"Quantiles",sep=""),mode="function")
		closeDialog()
		quantiles <- gsub(" +", ",", gsub(",", " ", tclvalue(quantilesVar)))
		if ("" == quantiles) {
			errorCondition(recall=nameVarF, message=gettextRcmdr("No probabilities specified."))
			return()
		}
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
		tail <- tclvalue(tailVar)
		pasteVar<-""
		for (i in 1:nnVar) {
			pasteVar<-paste(pasteVar,fVar$params[i],"=",vars[i],", ",sep="")
		}
		if (nameVar=="Gumbel") {
			doItAndPrint(paste("log(q",fVar$funName,"(c(", quantiles, "), ",pasteVar,
							"lower.tail=", tail == "lower",")) # Gumbel distribution", sep=""))
		} else {
			doItAndPrint(paste("q",fVar$funName,"(c(", quantiles, "), ",
							pasteVar,"lower.tail=", tail == "lower",")", sep=""))
		}
		putDialog(dialogName, list(initialValues=vars, tail=tclvalue(tailVar), quantiles=tclvalue(quantilesVar)), resettable=FALSE)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject=paste("q",fVar$funName,sep=""), reset = dialogName, apply=dialogName)
	tkgrid(labelRcmdr(entryFrame, text=gettextRcmdr("Probabilities")), quantilesEntry, sticky="w", padx=6)
	for (i in 1:nnVar) {
		tkgrid(labelRcmdr(entryFrame, text=gettextRcmdr(fVar$paramsLabels[i])), get(paramsEntry[i]), sticky="w", padx=6)
	}
	tkgrid(lowerTailButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Lower tail")), sticky="w")
	tkgrid(upperTailButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Upper tail")), sticky="w")
    tkgrid(entryFrame, sticky="w")
	tkgrid(buttonFrame, sticky="w")
	tkgrid.configure(quantilesEntry, sticky="w")
	for (i in 1:nnVar) {
		tkgrid.configure(get(paramsEntry[i]), sticky="w")
	}
    tkgrid(buttonsFrame, sticky="ew")
	dialogSuffix(focus=quantilesEntry)
}

distributionProbabilities <- function(nameVar){
	fVar<-get(paste(nameVar,"Distribution",sep=""))
	nnVar<-length(fVar$params)
	dialogName <- paste(nameVar,"Probabilities", sep="")
	defaults <- list(initialValues=fVar$initialValues, tail="lower", probabilities="")
	initial <- getDialog(dialogName, defaults=defaults)
	initializeDialog(title=gettextRcmdr(paste(fVar$titleName,"Probabilities",sep=" ")))
    entryFrame <- tkframe(top)
	probabilitiesVar <- tclVar(initial$probabilities)
	probabilitiesEntry <- ttkentry(entryFrame, width="30", textvariable=probabilitiesVar)
	paramsVar<-paste(fVar$params,"Var",sep="")
	paramsEntry<-paste(fVar$params,"Entry",sep="")
	for (i in 1:nnVar) {
		eval(parse(text=paste(paramsVar[i],"<-tclVar('", initial$initialValues[i],"')", sep="")))
		eval(parse(text=paste(paramsEntry[i],"<-ttkentry(entryFrame, width='6', textvariable=",paramsVar[i],")",sep="")))
	}
	tailVar <- tclVar(initial$tail)
    buttonFrame <- tkframe(top)
	lowerTailButton <- ttkradiobutton(buttonFrame, variable=tailVar, value="lower")
	upperTailButton <- ttkradiobutton(buttonFrame, variable=tailVar, value="upper")
	onOK <- function(){
		nameVarF<-get(paste(nameVar,"Probabilities",sep=""),mode="function")
		closeDialog()
		probabilities <-  gsub(" +", ",", gsub(",", " ", tclvalue(probabilitiesVar)))
		if ("" == probabilities) {
			errorCondition(recall=nameVarF, message=gettextRcmdr("No values specified."))
			return()
		}
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
		tail <- tclvalue(tailVar)
		pasteVar<-""
		for (i in 1:nnVar) {
			pasteVar<-paste(pasteVar,fVar$params[i],"=",vars[i],", ",sep="")
		}
		if (nameVar=="Gumbel") {
			doItAndPrint(paste("p",fVar$funName,"(exp(c(", probabilities, ")), ",
							pasteVar,"lower.tail=", tail == "lower",") #Gumbel Distribution", sep=""))
		} else {
			doItAndPrint(paste("p",fVar$funName,"(c(", probabilities, "), ",
							pasteVar,"lower.tail=", tail == "lower",")", sep=""))
		}
		tkfocus(CommanderWindow())
		putDialog(dialogName, list(initialValues=vars, tail=tclvalue(tailVar), probabilities=tclvalue(probabilitiesVar)), resettable=FALSE)
	}
	OKCancelHelp(helpSubject=paste("p",fVar$funName,sep=""), reset = dialogName, apply = dialogName)
	tkgrid(labelRcmdr(entryFrame, text=gettextRcmdr("Variable value(s)")), probabilitiesEntry, sticky="w", padx=6)
	for (i in 1:nnVar) {
		tkgrid(labelRcmdr(entryFrame, text=gettextRcmdr(fVar$paramsLabels[i])), get(paramsEntry[i]), sticky="w", padx=6)
	}
	tkgrid(lowerTailButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Lower tail")), sticky="w")
	tkgrid(upperTailButton, labelRcmdr(buttonFrame, text=gettextRcmdr("Upper tail")), sticky="w")
    tkgrid(entryFrame, sticky="w")
    tkgrid(buttonFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="ew")
	tkgrid.configure(probabilitiesEntry, sticky="w")
	for (i in 1:nnVar) {
		tkgrid.configure(get(paramsEntry[i]), sticky="w")
	}
	dialogSuffix(focus=probabilitiesEntry)
}

distributionMass  <- function(nameVar) {
  fVar<-get(paste(nameVar,"Distribution",sep=""))
  nnVar<-length(fVar$params)
  dialogName <- paste(nameVar,"Mass", sep="")
  defaults <- list(initialValues=fVar$initialValues)
  initial <- getDialog(dialogName, defaults=defaults)
  checkRange <- function(range){
    if (nameVar=="binomial") {
      messageVar<-"Number of trials, %d, is large.\nCreate long output?"
    } else {
      messageVar<-"Range of values over which to plot, %d, is large.\nCreate long output?"
    }
    RcmdrTkmessageBox(message=sprintf(gettextRcmdr(messageVar),range),
                      icon="warning", type="yesno", default="no")
  }
  initializeDialog(title=gettextRcmdr(paste(fVar$titleName,"Probabilities",sep=" ")))
  entryFrame <- tkframe(top)
  paramsVar<-paste(fVar$params,"Var",sep="")
  paramsEntry<-paste(fVar$params,"Entry",sep="")
  for (i in 1:nnVar) {
    eval(parse(text=paste(paramsVar[i],"<-tclVar('",initial$initialValues[i],"')",sep="")))
    eval(parse(text=paste(paramsEntry[i],"<-ttkentry(entryFrame, width='6', textvariable=",paramsVar[i],")",sep="")))
  }
  onOK <- function(){
    nameVarF<-get(paste(nameVar,"Mass",sep=""),mode="function")
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
    if (nameVar=="binomial") {
      if (vars[1] > 50){
        if ("no" == tclvalue(checkRange(vars[1]))){
          if (getRcmdr("grab.focus")) tkgrab.release(top)
          tkdestroy(top)
          nameVarF()
          return()
        }
      }
    } else {
      pasteVar<-""
      for (i in 1:nnVar) {
        pasteVar<-paste(pasteVar,", ",fVar$params[i],"=",vars[i], sep="")
      }
      xmin <- eval(parse(text=paste("q",fVar$funName,"(.0005",pasteVar,")",sep="")))
      xmax <- eval(parse(text=paste("q",fVar$funName,"(.9995",pasteVar,")",sep="")))
      range <- xmax-xmin
      if (xmax - xmin > 50){
        if ("no" == tclvalue(checkRange(range))){
          if (getRcmdr("grab.focus")) tkgrab.release(top)
          tkdestroy(top)
          nameVarF()
          return()
        }
      }
    }
    if (nameVar=="binomial") {
      command <- paste("local({\n  .Table <- data.frame(Probability=dbinom(0:", vars[1], 
                       ", size=", vars[1], 
                       ", prob=", vars[2], "))", sep="")
      command <- paste(command, "\n  rownames(.Table) <- 0:", vars[1], sep="")
    } else {
      command <- paste("local({\n  .Table <- data.frame(Probability=d",fVar$funName,
                       "(", xmin, ":", xmax, pasteVar, "))", sep="")
      command <- paste(command, "\n  rownames(.Table) <- ", xmin, ":", xmax, sep="")
    }
    command <- paste(command, "\n  print(.Table)\n})")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
    putDialog(dialogName, list(initialValues=vars), resettable=FALSE)
  }
  OKCancelHelp(helpSubject=paste("d",fVar$funName,sep=""), reset = dialogName, apply = dialogName)
  for (i in 1:nnVar) {
    tkgrid(labelRcmdr(entryFrame, text=gettextRcmdr(fVar$paramsLabels[i])), get(paramsEntry[i]), sticky="w", padx=6)
  }
  tkgrid(entryFrame, sticky="w")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  for (i in 1:nnVar) {
    tkgrid.configure(get(paramsEntry[i]), sticky="w")
  }
  dialogSuffix(focus=get(paramsEntry[1]))
}

setRandomSeed <- function(){
  initializeDialog(title = gettextRcmdr("Set Random Number Generator Seed"))
  seed <- sample(1e5, 1)
  seedValue <- tclVar(seed)
  seedSlider <- tkscale(top, from = 1, to = 1e5, showvalue = TRUE,
                    variable = seedValue, resolution = 1, orient = "horizontal")
  onOK <- function(){
    seed <- tclvalue(seedValue)
    doItAndPrint(paste("set.seed(", seed, ")", sep=""))
    closeDialog()                     
  }
  tkgrid(seedSlider, sticky="ew")
  OKCancelHelp(helpSubject = "set.seed")
  tkgrid(buttonsFrame, stick = "ew")
  dialogSuffix()
}
