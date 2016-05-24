##########################################################################
## SAVE structure class
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

## ----------------
## CLASS definition
## ----------------

## SAVE Class

setClass("SAVE", 		
         representation(
			#data
			 responsename="character",      ##name of the response
			 controllablenames="character", ##names of the controllable inputs
			 constant.controllables="logical", ##whether the controllable inputs are constant T or not F. Default in the Code is F.
			 calibrationnames="character",  ##names of the calibration inputs
			 df="matrix",                   ##field design (no replicates) 
			 dm="matrix",                   ##model design
			 ym="numeric",                  ##model response associated with dm
			 yf="numeric",                  ##field response associated with df
			#mle fit
			 meanformula="formula",         ##linear model for the mean of the GP
			 mle="list",                    ##three components: thetaL, thetaM, thetaF
			 bestguess="numeric",           ##bestguess
			 xm="matrix",				   ##design matrix associated with mean formula
			 xf="matrix",				   ##design matrix associated with
			                               ##=model.matrix(mean.formula,x.unique)
			#bayesian options
			 prior="matrix",
			 method="numeric",
			 mcmcMultmle="numeric",
			#bayesian fit
			 mcmcsample="matrix",
			#other
			 wd="character",
			 call="call",
			 bayesfitcall="call"
			)
		)
		
validSAVEobject <- function(object) {
		if (is.null(object@responsename) || length(object@responsename)==0)
		{return("Response names are not provided\n")}
		
		if (is.null(object@controllablenames) || length(object@controllablenames)==0)
		{return("Controllable names are not provided\n")}
		
		if (is.null(object@yf) || length(object@yf)==0)
		{return("field data are not provided\n")}
		
		if (is.null(object@ym) || length(object@ym)==0)
		{return("Model data are not provided\n")}
		TRUE
	}
setValidity("SAVE",validSAVEobject)


setClass("summary.SAVE", representation(
		#only STAGEI:
		stage2 = "logical",
		call = "call",
		mle = "list",
		#STAGEII also:
		bayesfitcall = "call",
		mcmcsum = "matrix"
		)
	)

setClass("predictcode.SAVE", representation(
		predictcodecall = "call",
		mle = "list",
		samples = "data.frame",
		newdesign = "data.frame",
		modelmean = "vector",
		covmat = "matrix"
		)
	)

setClass("summary.predictcode.SAVE", representation(
    		call = "call",
		summariesmodelpred = "matrix",
		modelmean = "vector"
		)
	)

setClass("predictreality.SAVE", representation(
		#first half of samples is prediction and second half is bias
		predictrealitycall = "call",
		modelpred = "data.frame",
		biaspred = "data.frame",
		newdesign = "data.frame"
		)
	)

setClass("summary.predictreality.SAVE", representation(
		call = "call",
		biascorr = "matrix",
		biaspred = "matrix"
		)
	)


setClass("validate.SAVE", representation(
		call = "call",
		bayesfitcall = "call",
		validatecall = "call",
		newdesign = "data.frame",
		validate = "matrix"
		)
	)
	
setClass("summary.validate.SAVE", representation(
		callSAVE = "call",
		callbayes = "call",
		callvalidate = "call",
		summaries = "matrix"
		)
	)
