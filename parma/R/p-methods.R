#################################################################################
##
##   R package parma
##   Alexios Ghalanos Copyright (C) 2012-2013 (<=Aug)
##   Alexios Ghalanos and Bernhard Pfaff Copyright (C) 2013- (>Aug)
##   This file is part of the R package parma.
##
##   The R package parma is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package parma is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
parmautility = function(U = c("CARA", "Power"), method = c("moment", "scenario"), 
		scenario = NULL, M1 = NULL, M2 =  NULL, M3 = NULL, M4 = NULL, RA = 1, 
		budget = 1, LB = rep(0, length(M1)), UB = rep(1, length(M1)))
{
	UseMethod("parmautility")
}

setMethod("parmautility", definition = .parmautility)

parmaspec = function(scenario = NULL, probability = NULL, S = NULL, Q = NULL, qB = NULL,
		benchmark = NULL, benchmarkS = NULL, forecast = NULL, target = NULL, 
		targetType =  c("inequality", "equality"), 
		risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM", "LPMUPM"), 
		riskType = c("minrisk", "optimal", "maxreward"), riskB = NULL,
		options = list(alpha = 0.05, threshold = 999, moment = 1, 
				lmoment = 1, umoment = 1, lthreshold = -0.01, uthreshold = 0.01), 		
		LB = NULL, UB = NULL, budget = 1, leverage = NULL, 
		ineqfun = NULL, ineqgrad = NULL, eqfun = NULL, eqgrad = NULL, 
		uservars = list(), ineq.mat = NULL, ineq.LB = NULL, 
		ineq.UB = NULL, eq.mat = NULL, eqB = NULL, max.pos = NULL, 
		asset.names = NULL, ...)
{
	UseMethod("parmaspec")
}

setMethod("parmaspec", definition = .parmaspec)


parmasolve = function(spec, type = NULL, solver = NULL, solver.control = list(), 
		x0 = NULL, w0 = NULL, parma.control = list(ubounds = 1e4, 
				mbounds = 1e5, penalty = 1e4, eqSlack = 1e-5), ...)
{
	UseMethod("parmasolve")
}

setMethod("parmasolve", signature(spec = "parmaSpec"), .parmasolve)

parmafrontier = function(spec, n.points = 100, miny = NULL, maxy = NULL, 
		type = NULL, solver = NULL, solver.control = list(), 
		parma.control = list(ubounds = 10000, mbounds = 1e+05, penalty = 10000), 
		cluster = NULL)
{
	UseMethod("parmafrontier")
}

setMethod("parmafrontier", signature(spec = "parmaSpec"), .parmafrontier)


.parmaweights = function(object, ...)
{
	w = object@solution$weights
	names(w) = object@model$asset.names
	return( w )
}
setMethod("weights", signature(object = "parmaPort"), .parmaweights)

tictoc = function(object, ...){
	UseMethod("tictoc")
}

.tictoc = function(object, ...){
	return(object@model$elapsed)
}
setMethod("tictoc", signature(object = "parmaPort"), .tictoc)

.parmaweights = function(object, ...)
{
	w = object@solution$weights
	names(w) = object@model$asset.names
	return( w )
}
setMethod("weights", signature(object = "parmaPort"), .parmaweights)



checkarbitrage = function(object){
	UseMethod("checkarbitrage")
}
# the function
.arbcheck = function(weights, S, options, risk){
	AC = c(0, 0)
	if(tolower(risk) == "lpmupm"){
		f = fun.lpm(weights, S, options$lthreshold, options$lmoment)
	} else{
		f = fun.risk(weights, S, options, risk)
	}
	if(f<=0) AC[1] = 1
	mx = min(as.numeric(S %*% weights))
	if(mx>=0) AC[2] = 1
	# should also look at correlation extemes (-0.9, 0.9)
	return(AC)
}
# the extractor
.checkarbitrage = function(object){
	x = object@solution$arbitrage
	x = as.logical(x)
	names(x) = c("Risk(<=0)", "ALL[Port(>=0)]")
	return(x)
}

setMethod("checkarbitrage", signature(object = "parmaPort"), .checkarbitrage)

# show method

parmarisk = function(object, ...)
{
	UseMethod("parmarisk")
}

.parmarisk = function(object, ...){
	ans = object@solution$risk
	names(ans) = object@model$risk
	return(ans)
}

setMethod("parmarisk", signature(object = "parmaPort"), .parmarisk)

parmareward = function(object, ...)
{
	UseMethod("parmareward")
}

.parmareward = function(object, ...){
	ans = object@solution$reward
	return(ans)
}

setMethod("parmareward", signature(object = "parmaPort"), .parmareward)

setMethod("show",
		signature(object = "parmaSpec"),
		function(object){
			cat(paste("\n+---------------------------------+", sep = ""))
			cat(paste("\n|       PARMA Specification       |", sep = ""))
			cat(paste("\n+---------------------------------+ ", sep = ""))
			cat(paste("\nNo.Assets\t\t: ", object@model$indx[8], sep = ""))
			tmp = c("LP", "MILP", "QP", "MIQP", "SOCP", "NLP", "MINLP", "GNLP")[which(object@model$type==1)]
			cat("\nProblem\t\t\t:", paste(tmp, sep=" ", collapse = ","))
			cat(paste("\nInput\t\t\t: ", if(object@model$indx[1]==1) "Scenario" else "Covariance", sep = ""))
			cat(paste("\nRisk Measure\t: ", object@model$risk, sep = ""))
			cat(paste("\nObjective\t\t: ", object@model$riskType, sep = ""))
			cat("\n\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "parmaPort"),
		function(object){
			if(!is.null(object@model$utility)){
				w = weights(object) 
				cat(paste("\n+---------------------------------+", sep = ""))
				cat(paste("\n|        PARMA Portfolio          |", sep = ""))
				cat(paste("\n+---------------------------------+ ", sep = ""))
				cat(paste("\nNo.Assets\t\t: ", length(w), sep = ""))
				cat(paste("\nProblem\t\t\t: NLP"))
				cat(paste("\nUtility\t\t\t: CARA"))
				cat(paste("\nn.moments\t\t: ", object@model$moments, sep = ""))
				cat(paste("\nObjective\t\t: ", round(object@solution$utility,4), sep = ""))
				cat(paste("\nReward\t\t\t: ", round(object@solution$reward,7), sep = ""))
				cat("\n\n")
				w = weights(object)
				cnames = object@model$asset.names
				if(is.null(cnames)) cnames = paste("A_",1:length(w), sep = "")
				names(w) = cnames
				wx = sort(w, decreasing = TRUE)
				wx = wx[which(abs(wx)>0.001)]
				w = data.frame(Optimal_Weights = wx)
				print(round(w, 4))
			} else{
				cat(paste("\n+---------------------------------+", sep = ""))
				cat(paste("\n|        PARMA Portfolio          |", sep = ""))
				cat(paste("\n+---------------------------------+ ", sep = ""))
				cat(paste("\nNo.Assets\t\t: ", object@model$indx[8], sep = ""))
				cat(paste("\nProblem\t\t\t: ", object@model$type, sep = ""))
				cat(paste("\nRisk Measure\t: ", object@model$risk, sep = ""))
				cat(paste("\nObjective\t\t: ", object@model$riskType, sep = ""))
				cat(paste("\nRisk\t\t\t: ", round(object@solution$risk,7), sep = ""))
				cat(paste("\nReward\t\t\t: ", round(object@solution$reward,7), sep = ""))
				cat("\n\n")
				w = weights(object)
				cnames = object@model$asset.names
				if(is.null(cnames)) cnames = paste("A_",1:length(w), sep = "")
				names(w) = cnames
				wx = sort(w, decreasing = TRUE)
				wx = wx[which(abs(wx)>0.001)]
				w = data.frame(Optimal_Weights = wx)
				print(round(w, 4))
			}
			invisible(object)
		})

setGeneric("parmaset<-", function(object, value){standardGeneric("parmaset<-")})

parmaset<-function(object, value){
	UseMethod("parmaset")
}


.parmaset<-function(object, value){
	parnames = names(value)
	valid.choices = c("scenario", "probability", "S", "benchmark", "benchmarkS", 
			"forecast", "target", "targetType", "risk", "riskType", "options", 
			"LB", "UB", "budget", "leverage", "ineqfun", "ineqgrad", "eqfun", 
			"eqgrad", "uservars", "ineq.mat", "ineq.LB", 
			"ineq.UB", "eq.mat", "eqB", "max.pos", "asset.names")
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], valid.choices))){
			warning( (paste("Unrecognized Value: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	indx = object@model$indx
	scenario = object@modeldata$scenario
	probability = object@modeldata$probability
	S = object@modeldata$S
	# need to do this else benchmarkS will be mistakenly added (spec sets it to zero
	# which is NOT NULL)
	if(indx[2]==1) benchmark = object@modeldata$benchmark else benchmark = NULL
	if(indx[2]==1) benchmarkS = object@modeldata$benchmarkS else benchmarkS = NULL
	forecast = object@modeldata$forecast
	target = object@modeldata$target
	targetType = object@model$targetType
	risk = object@model$risk
	riskType = object@model$riskType
	options = object@model$options
	LB = object@constraints$LB
	UB = object@constraints$UB
	budget = object@constraints$budget
	leverage = object@constraints$leverage
	ineqfun = object@constraints$ineqfun
	ineqgrad = object@constraints$ineqgrad
	eqfun = object@constraints$eqfun
	eqgrad = object@constraints$eqgrad
	uservars = object@constraints$uservars
	ineq.mat = object@constraints$ineq.mat
	ineq.LB = object@constraints$ineq.LB
	ineq.UB = object@constraints$ineq.UB
	eq.mat = object@constraints$eq.mat
	eqB = object@constraints$eqB
	max.pos = object@constraints$max.pos
	asset.names = object@modeldata$asset.names
	
	if(length(inc)>0){
		for(i in 1:length(inc)){
			eval(parse(text = paste(parnames[inc[i]],"=unname(value[[inc[i]]])",sep="")))
		}
	}
	newspec = parmaspec(scenario = scenario, probability = probability, 
			S = S, benchmark = benchmark, benchmarkS = benchmarkS, 
			forecast = forecast, target = target, targetType =  targetType, 
			risk = risk, riskType = riskType, options = options, 
			LB = LB, UB = UB, budget = budget, leverage = leverage, 
			ineqfun = ineqfun, ineqgrad = ineqgrad, eqfun = eqfun, eqgrad = eqgrad, 
			uservars = uservars, ineq.mat = ineq.mat, ineq.LB = ineq.LB, 
			ineq.UB = ineq.UB , eq.mat = eq.mat, eqB = eqB, maxp.pos = max.pos, 
			asset.names =  asset.names)
	return(newspec)
}

setReplaceMethod(f="parmaset", signature= c(object = "parmaSpec", value = "vector"), definition = .parmaset)


parmaget = function(object, value){
	UseMethod("parmaget")
}

.parmaget<-function(object, value){
	valid.choices = c("scenario", "probability", "S", "benchmark", "benchmarkS", 
			"forecast", "target", "targetType", "risk", "riskType", 
			"options", "LB", "UB", "budget", "leverage", "ineqfun", "ineqgrad", 
			"eqfun", "eqgrad", "uservars", "ineq.mat", "ineq.LB", 
			"ineq.UB", "eq.mat", "eqB", "max.pos", "asset.names")
	inc = match.arg(value[1], valid.choices)
	if(is.na(inc)) stop("\nunrecognized value...try again (case sensitive)")
	scenario = object@modeldata$scenario
	probability = object@modeldata$probability
	
	S = object@modeldata$S
	# need to do this else benchmarkS will be mistakenly added (spec sets it to zero
	# which is NOT NULL)
	if(object@model$indx[2]==1) benchmark = object@modeldata$benchmark else benchmark = NULL
	if(object@model$indx[2]==1) benchmarkS = object@modeldata$benchmarkS else benchmarkS = NULL
	forecast = object@modeldata$forecast
	target = object@modeldata$target
	targetType = object@model$targetType
	risk = object@model$risk
	riskType = object@model$riskType
	options = object@model$options
	LB = object@constraints$LB
	UB = object@constraints$UB
	budget = object@constraints$budget
	leverage = object@constraints$leverage
	ineqfun = object@constraints$ineqfun
	ineqgrad = object@constraints$ineqgrad
	eqfun = object@constraints$eqfun
	eqgrad = object@constraints$eqgrad
	uservars = object@constraints$uservars
	ineq.mat = object@constraints$ineq.mat
	ineq.LB = object@constraints$ineq.LB
	ineq.UB = object@constraints$ineq.UB
	eq.mat = object@constraints$eq.mat
	eqB = object@constraints$eqB
	asset.names = object@modeldata$asset.names
	max.pos = object@constraints$max.pos
	ans = NA
	eval(parse(text = paste("ans=",inc,sep="")))
	return(ans)
}
setMethod("parmaget", signature(object = "parmaSpec"), definition = .parmaget)


parmastatus = function(object){
	setMethod("parmastatus")
}

.pstatus = function(object){
	object@solution$status
}

setMethod("parmastatus", signature(object = "parmaPort"), definition = .pstatus)
