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
.parmaspec = function(scenario = NULL, probability = NULL, S = NULL, Q = NULL,
		qB = NULL, benchmark = NULL, benchmarkS = NULL, forecast = NULL, target = NULL, 
		targetType =  c("inequality", "equality"), 
		risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM", "LPMUPM"), 
		riskType = c("minrisk", "optimal", "maxreward"), riskB = NULL,
		options = list(alpha = 0.05, threshold = 999, moment = 1, 
				lmoment = 1, umoment = 1, lthreshold = -0.01, uthreshold = 0.01), 
		LB = NULL, UB = NULL, budget = 1, leverage = NULL, 
		ineqfun = NULL, ineqgrad = NULL, eqfun = NULL, eqgrad = NULL, 
		uservars = list(), ineq.mat = NULL, ineq.LB = NULL, 
		ineq.UB = NULL, eq.mat = NULL, eqB = NULL, max.pos = NULL, 
		asset.names = NULL, ...){
	
	# for future upgrades:
	hasminlp = 0
	hasqcqp = 0
	hasmiqp = 0
	##### determine problem type #####################
	
	### a bunch of checks
	indx = rep(0, 8)
	names(indx) = c("datatype", "benchmark", "targettype", "risk", "risktype", "leverage", "aux1", "aux2")
	# datatype: 1 = scenario, 2 = S
	# benchmark: 1 (TRUE), 0 (FALSE)
	# targettype: 1 = inequality, 2 = equality
	# risk: 1 (MAD), 2 (MiniMax), 3 (CVaR), 4 (CDaR), 5(EV), 6 (LPM), 7 (LPMUPM)
	# risktype: 1 = minrisk, 2 = optimal (fractional), 3 = maxreward (only for EV)
	# leverage: !=0
	type = rep(0, 8)
	# problem types c("LP", "MILP", "QP", "MIQP", "SOCP", "NLP", "MINLP", "GNLP")
	
	if(is.null(scenario) && is.null(S)) stop("\nparma: You cannot have both the scenario AND covariance matrix (S) NULL!")
	if(!is.null(scenario)){
		if(tolower(riskType)=="maxreward") stop("\nparma: maxreward type only supported for covariance matrix (S) at present.")
		S = NULL #  set this to NULL in case both are not NULL!
		scenario = as.matrix(scenario)
		indx[1] = 1
		m = NCOL(scenario)
		n = NROW(scenario)
		indx[8] = m
		if(is.null(asset.names)) asset.names = colnames(scenario)
		if(is.null(forecast)){
			if(!is.null(benchmark)){
				warning("\nparma: no forecast provided but benchmark is included...using means from scenario - benchmark instead")
				forecast = as.numeric( colMeans(scenario) - mean(as.vector(benchmark)) )
			} else{
				forecast = as.numeric( colMeans(scenario) )
				warning("\nparma: no forecast provided...using means from scenario instead.")
			}
		} else{
			forecast = as.numeric(forecast)[1:m]
		}
		if(!is.null(benchmark)){		
			if(length(benchmark)<n){
				benchmark = rep(benchmark[1], n)
			} else{
				nb = length(as.numeric(benchmark))
				if(nb!=n) stop("\nparma: benchmark length not equal to scenario length.")			
				benchmark = as.numeric(benchmark[1:n])
			}
			#mbenchmark = mean(benchmark)
			indx[2] = 1
		} else{
			benchmark = rep(0, n)
			#mbenchmark = 0
		}

		if(is.null(target)){
			target = 0
			if(tolower(riskType[1]) == "minrisk") warning("\nparma: no target provided...setting target reward to zero.")
		} else{
			target = as.numeric(target)[1]
		}
		# Probability only currently support for LP problems
		if(is.null(probability)){
			probability = rep(1/n, n)
		} else{
			probability = probability[1:n]
			if(sum(probability)!=1) warning("\nProbability does not sum to 1!")
		}
		
		xtmp =  match(tolower(riskType[1]), c("minrisk", "optimal"))
		if(is.na(xtmp)) stop("\nparma: riskType not recognized") else riskType = c("minrisk", "optimal")[xtmp]
		indx[5] = xtmp
		
		tmp =  match(tolower(targetType[1]), c("inequality", "equality"))
		if(is.na(tmp)) stop("\nparma: targetType not recognized") else targetType = c("inequality", "equality")[tmp]
		indx[3] = tmp
		tmp = match(tolower(risk[1]), tolower(c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM", "LPMUPM")))
		if(is.na(tmp)) stop("\nparma: risk not recognized") else risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM", "LPMUPM")[tmp]
		indx[4] = tmp

		if(!is.null(leverage)) indx[6] = as.numeric( leverage )
		
		# Some initial checks
		if(!is.null(ineqfun) && (!is.null(ineq.mat) | !is.null(eq.mat))){
			stop("\nparma: you cannot mix ineqfun with linear custom constraints")
		}
		if(!is.null(eqfun) && (!is.null(ineq.mat) | !is.null(eq.mat))){
			stop("\nparma: you cannot mix eqfun with linear custom constraints")
		}
		if(!is.null(leverage) && (!is.null(ineq.mat) | !is.null(eq.mat))){
			stop("\nparma: you cannot mix leverage (NLP) with linear custom constraints")
		}
		if(!is.null(ineqfun)){
			chk = .checkconsfun(ineqfun, name = "ineqfun")
		}
		if(!is.null(ineqgrad)){
			chk = .checkconsfun(ineqgrad, name = "ineqgrad")
		}
		if(!is.null(eqfun)){
			chk = .checkconsfun(eqfun, name = "eqfun")
		}
		if(!is.null(eqgrad)){
			chk = .checkconsfun(eqgrad, name = "eqgrad")
		}
		
		
		# Go trough each problem in turn (simpler but redundant)
		midx = widx = vidx = NA
		# MAD [ weights[m] ... ]
		if(tmp == 1){
			widx = 1:m
			if(riskType == "optimal") midx=m+1
			if(!is.null(max.pos)){
				if(riskType == "minrisk"){
					type[c(2, 8)] = 1
					if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
						type[2] = 0
					}
					if(!is.null(ineq.mat) | !is.null(eq.mat)){
						type[8] = 0
					}
					if(hasminlp && type[2]>0) type[7] = 1
				} else{
					if(!is.null(ineq.mat) | !is.null(eq.mat))
						stop("\nparma: you cannot have ineq.mat or eq.mat (LP constraints) with cardinality constraints in optimal risk case (MINLP/GNLP problem)!")
					if(hasminlp) type[c(7,8)] = 1 else type[8] = 1
				}
			} else{
				type[c(1, 6, 8)] = 1
				if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
					type[1] = 0
				}
				if(!is.null(ineq.mat) | !is.null(eq.mat)){
					type[c(6,8)] = 0
				}
			}
		}
		# MiniMax [ weights[m] max.loss[1] ... ]
		if(tmp == 2){
			if(riskType=="optimal"){
				vidx = 1
				widx = 2:(m+1)
				midx = m+2
			} else{
				vidx = 1
				widx = 2:(m+1)
				midx = 0
			}
			if(!is.null(max.pos)){
				if(riskType == "minrisk"){
					type[c(2, 8)] = 1
					if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
						type[2] = 0
					}
					if(!is.null(ineq.mat) | !is.null(eq.mat)){
						type[8] = 0
					}
					if(hasminlp && type[2]>0) type[7] = 1
				} else{
					if(!is.null(ineq.mat) | !is.null(eq.mat))
						stop("\nparma: you cannot have ineq.mat or eq.mat (LP constraints) with cardinality constraints in optimal risk case (MINLP/GNLP problem)!")
					if(hasminlp) type[c(7,8)] = 1 else type[8] = 1
				}
			} else{
				type[c(1, 6, 8)] = 1
				if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
					type[1] = 0
				}
				if(!is.null(ineq.mat) | !is.null(eq.mat)){
					type[c(6,8)] = 0
				}
			}
		}
		
		# CVaR [ VaR[1] weights[m] ... ]
		if(tmp == 3){
			widx = 2:(m+1)
			vidx = 1
			if(riskType == "optimal") midx=m+2
			if(!is.null(max.pos)){
				if(riskType == "minrisk"){
					type[c(2, 8)] = 1
					if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
						type[2] = 0
					}
					if(!is.null(ineq.mat) | !is.null(eq.mat)){
						type[8] = 0
					}
					if(hasminlp && type[2]>0) type[7] = 1
				} else{
					if(!is.null(ineq.mat) | !is.null(eq.mat))
						stop("\nparma: you cannot have ineq.mat or eq.mat (LP constraints) with cardinality constraints in optimal risk case (MINLP/GNLP problem)!")
					if(hasminlp) type[c(7,8)] = 1 else type[8] = 1
				}
			} else{
				type[c(1, 6, 8)] = 1
				if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
					type[1] = 0
				}
				if(!is.null(ineq.mat) | !is.null(eq.mat)){
					type[c(6,8)] = 0
				}
			}
		}
		
		# CDaR [ DaR[1] weights[m] ... ]
		if(tmp == 4){
			widx = 2:(m+1)
			vidx = 1
			if(riskType == "optimal") midx=m+2
			# scenario && CDaR = LP (no NLP yet)
			if(!is.null(max.pos)){
				if(riskType == "minrisk") type[2] = 1 else stop("\nparma: CDaR with cardinality constraints and optimal risk type not supported.")
			} else{
				type[1] = 1
			}
			if(!is.null(leverage)) stop("\nparma: CDaR with leverage requires NLP formulation (not supported)")
			if(!is.null(ineqfun)) stop("\nparma: CDaR with custom NLP (ineqfun) constraints requires NLP formulation (not supported)")
			if(!is.null(eqfun)) stop("\nparma: CDaR with custom NLP (eqfun) constraints requires NLP formulation (not supported)")
		}
	
		# EV [ weights[m] ... ]
		if(tmp == 5){
			widx = 1:m
			if(!is.null(max.pos)){
				if(hasminlp) type[c(7,8)] = 1 else type[8] = 1
			} else{
				type[c(6, 8)] = 1
			}
			if(riskType == "optimal") midx = m+1
			if(!is.null(ineq.mat) | !is.null(eq.mat)){
				stop("\nparma: EV scenario problem is NLP. Does not accept linear constraints (use QP with S matrix instead).")
			}
		}
		# LPM [ weights[m] ... ]
		if(tmp == 6){
			# NO benchmark for LPM
			benchmark = 0
			widx = 1:m
			if(riskType == "optimal") midx=m+1
			# LPM[1] has both LP and NLP
			if(options$moment == 1){
				if(!is.null(max.pos)){
					if(riskType == "minrisk"){
						type[c(2, 8)] = 1
						if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
							type[2] = 0
						}
						if(!is.null(ineq.mat) | !is.null(eq.mat)){
							type[8] = 0
						}
						if(hasminlp && type[2]>0) type[7] = 1
					} else{
						if(!is.null(ineq.mat) | !is.null(eq.mat))
							stop("\nparma: you cannot have ineq.mat or eq.mat (LP constraints) with cardinality constraints in optimal risk case (MINLP/GNLP problem)!")
						if(hasminlp) type[c(7,8)] = 1 else type[8] = 1
					}
				} else{
					type[c(1, 6, 8)] = 1
					if(!is.null(leverage) | !is.null(ineqfun) | !is.null(eqfun)){
						type[1] = 0
					}
					if(!is.null(ineq.mat) | !is.null(eq.mat)){
						type[c(6,8)] = 0
					}
				}
			} else{
				# LPM[!=1] only NLP (though NLP[2] also has a QP not implemented here)
				if(!is.null(max.pos)){
					if(hasminlp) type[c(7,8)]<-1 else type[8]<-1
				} else{
					type[c(6, 8)] = 1
				}
				if(!is.null(ineq.mat) | !is.null(eq.mat)){
					stop("\nparma: LPM scenario problem with moment!=1 is NLP. Does not accept linear constraints.")
				}
			}
		}
		
		if(tmp == 7){
			widx = 1:m
			midx = m+1
			type[8]<-1
			if(!is.null(ineq.mat) | !is.null(eq.mat)){
				stop("\nparma: LPMUPM is GNLP. Does not accept linear constraints.")
			}
			riskType = "optimal"
			indx[5] = 2
		}
				
		if(is.null(LB)){
			LB = rep(0, m)
			warning("\nparma: no LB provided...setting Lower Bounds to zero.")
		}
		if(is.null(UB)){
			UB = rep(1, m)
			warning("\nparma: no UB provided...setting Upper Bounds to 1.")
		}
		if(any(UB<LB)) stop("\nparma: UB must be greater than LB.")
		if(is.null(leverage)){
			if(is.null(budget)){
				budget = 1
				warning("\nparma: no budget (or leverage) provided...setting budget constraint to 1.")
			}
		} else{
			indx[6] = 1
		}
		if(!is.null(ineq.mat)){
			ineq.mat = as.matrix(ineq.mat)
			nb = dim(ineq.mat)[2]
			if(nb!=m) stop("\nparma: ineq.mat columns not equal to number of assets")
			nb = dim(ineq.mat)[1]
			if(nb!=length(ineq.LB)) stop("\nparma: ineq.mat rows not equal to length of ineq.LB")
			if(nb!=length(ineq.UB)) stop("\nparma: ineq.mat rows not equal to length of ineq.UB")
		}
		if(!is.null(eq.mat)){
			ineq.mat = as.matrix(eq.mat)
			nb = dim(eq.mat)[2]
			if(nb!=m) stop("\nparma: eq.mat columns not equal to number of assets")
			nb = dim(eq.mat)[1]
			if(nb!=length(eqB)) stop("\nparma: eq.mat rows not equal to length of eqB")
		}
		# Enforce GNLP for custom constraints without jacobians
		if(!is.null(ineqfun)){
			if(is.null(ineqgrad)){
				type = rep(0,8)
				type[8] = 1
			}
		}
		if(!is.null(eqfun)){
			if(is.null(eqgrad)){
				type = rep(0,8)
				type[8] = 1
			}
		}
	} else{
		# QP [ v[1] weights[m] ... ]
		m = NCOL(S)
		benchmark = ineqfun = eqfun = ineqgrad = eqgrad = NULL
		midx = 1
		vidx = 0
		widx = 2:(m+1)
		
		if( !is.null(max.pos) ){
			if(hasmiqp){
				type[4] = 1
			} else{
				stop("\nparma: max.pos NOT NULL but MIQP not available!")
			}
		} else{
			# Q is a list of matrices
			if(!is.null(Q)){
				type[5] = 1
				if(riskType == "optimal") stop("\nparma: QCQP not yet implemented for optimal risk problem")
				#if(ncol(Q)!=m | nrow(Q)!=m) stop("\nparma: Q matrix dimensions must be the same as S.")
			} else{
				type[3] = 1
				type[5] = 1
			}
		}
		if(!is.null(leverage)){
			if(leverage<=0) stop("\nparma: leverage must be strictly positive!")
			# leverage only allowed in SOCP formulation
			type[3] = 0
			type[5] = 1
			indx[6] = as.numeric( leverage )
		}
		if(tolower(riskType)=="maxreward"){
			# maxreward only supported by SOCP
			type[3] = 0
			type[5] = 1
		}
		if(riskType == "optimal") midx=1
		
		if(tolower(risk[1])!="ev") stop("\nparma: only EV risk type allowed with covariance matrix (S)")
		indx[1] = 2
		indx[8] = m
		if(is.null(asset.names)) asset.names = colnames(S)
		if(!is.null(benchmarkS)){
			benchmarkS = matrix(benchmarkS, nrow = 1)
			indx[2] = 1
			nb = NCOL(benchmarkS)
			if((nb-1)!=m) stop("\nparma: benchmarkS length must be equal to ncol S + 1.")
		} else{
			benchmarkS = matrix(0, nrow = 1, ncol = m+1)
		}
		if(is.null(forecast)){
			forecast = rep(0, m)
			warning("\nparma: no forecast provided...setting to zero.")
		} else{
			forecast = as.numeric(forecast)[1:m]
		}
		if(is.null(target)){
			target = 0
			if(tolower(riskType[1]) == "minrisk"){
				warning("\nparma: no target provided...setting target reward to zero.")
				indx[5] = 1
			} else if(tolower(riskType[1])=="optimal"){
				indx[5] = 2
			} else{
				indx[5] = 3
				if(is.null(riskB)) stop("\nparma: maxreward option chosen but riskB not provided!")
				riskB = as.numeric(riskB[1])
			}
		} else{
			target = as.numeric(target)[1]
			if(tolower(riskType[1]) == "minrisk"){
				indx[5] = 1
			} else if(tolower(riskType[1])=="optimal"){
				indx[5] = 2
			} else{
				indx[5] = 3
				if(is.null(riskB)) stop("\nparma: maxreward option chosen but riskB not provided!")
				riskB = as.numeric(riskB[1])
				target = NULL
				warning("\nparma: maxreward chosen AND target given...setting target to NULL.")
			}
		}
		tmp =  match(tolower(targetType[1]), c("inequality", "equality"))
		if(is.na(tmp)) stop("\nparma: targetType not recognized") else targetType = c("inequality", "equality")[tmp]
		indx[3] = tmp
		if(is.null(LB)){
			LB = rep(0, m)
			warning("\nparma: no LB provided...setting Lower Bounds to zero.")
		}
		if(is.null(UB)){
			UB = rep(1, m)
			warning("\nparma: no UB provided...setting Upper Bounds to 1.")
		}
		if(any(UB<LB)) stop("\nparma: UB must be greater than LB.")
		if(is.null(leverage)){
			if(is.null(budget)){
				budget = 1
				warning("\nparma: no budget (or leverage) provided...setting budget constraint to 1.")
			}
		}
		if(!is.null(ineq.mat)){
			ineq.mat = as.matrix(ineq.mat)
			nb = dim(ineq.mat)[2]
			if(nb!=m) stop("\nparma: ineq.mat columns not equal to number of assets")
			nb = dim(ineq.mat)[1]
			if(nb!=length(ineq.LB)) stop("\nparma: ineq.mat rows not equal to length of ineq.LB")
			if(nb!=length(ineq.UB)) stop("\nparma: ineq.mat rows not equal to length of ineq.UB")
		}
		if(!is.null(eq.mat)){
			eq.mat = as.matrix(eq.mat)
			nb = dim(eq.mat)[2]
			if(nb!=m) stop("\nparma: eq.mat columns not equal to number of assets")
			nb = dim(eq.mat)[1]
			if(nb!=length(eqB)) stop("\nparma: eq.mat rows not equal to length of eqB")
		}
	}
	model = list(indx = indx, risk = risk, riskType = riskType, targetType = targetType,
			options = options, type = type, widx = widx, midx = midx, vidx = vidx)
	
	modeldata = list(scenario = scenario, probability = probability, S = S, Q = Q,
			qB = qB, benchmark = benchmark, benchmarkS = benchmarkS, forecast = forecast, 
			target = target, riskB = riskB, asset.names = asset.names, uservars = uservars)
	
	constraints = list(LB = LB, UB = UB, budget = budget, leverage = leverage, 
			ineqfun = ineqfun, ineqgrad = ineqgrad, eqfun = eqfun, eqgrad = eqgrad, 
			ineq.mat = ineq.mat, ineq.LB = ineq.LB, ineq.UB = ineq.UB, 
			eq.mat = eq.mat, eqB = eqB, uservars = uservars, max.pos = max.pos)
	
	ans = new("parmaSpec", 
			model = model, 
			modeldata = modeldata, 
			constraints = constraints)
	return(ans)
}

# REM index = c("datatype", "benchmark", "targettype", "risk", "risktype", "leverage", "aux1", "aux2")
# targettype: 1 = inequality, 2 = equality

.spec2minNLP = function(spec){
	optvars = list()
	optvars$widx = spec@model$widx
	optvars$midx = spec@model$midx
	optvars$vidx = spec@model$vidx
	optvars$index = spec@model$indx
	# MAD and EV are deviations from the mean. In the presence of a benchmark however,
	# this is not subtracted
	if(spec@model$indx[2]==0){
		optvars$Data = switch(tolower(spec@model$risk), 
				"mad" = scale(spec@modeldata$scenario, scale = FALSE),
				"ev" = scale(spec@modeldata$scenario, scale = FALSE),
				"minimax" = spec@modeldata$scenario,
				"cvar" = spec@modeldata$scenario,
				"cdar" = spec@modeldata$scenario,
				"lpm"  = spec@modeldata$scenario,
				spec@modeldata$scenario)
	} else{
		optvars$Data = spec@modeldata$scenario
	}
	optvars$benchmark  = spec@modeldata$benchmark
	#optvars$mbenchmark = spec@modeldata$mbenchmark
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	
	risk = spec@model$risk
	mn = dim(optvars$Data)
	N = mn[1]
	m = mn[2]
	
	optvars$wm = m
	optvars$N  = N
	if(tolower(risk) == "lpm"){
		if(spec@model$options$threshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$threshold = 0
		} else{
			optvars$threshold = spec@model$options$threshold
		}
	} else{
		optvars$threshold = 0
	}
	optvars$moment = spec@model$options$moment
			
	if(tolower(risk) == "cvar"){
		optvars$fm = m+1
		x0 = as.numeric(c(quantile(optvars$Data %*% rep(1/m, m), spec@model$options$alpha), rep(1/m, m)))
		optvars$LB = c(-10, spec@constraints$LB)
		optvars$UB = c(  0, spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else if(tolower(risk) == "minimax"){
		optvars$fm = m+1
		x0 = as.numeric(c(-min(optvars$Data %*% rep(1/m, m)), rep(1/m, m)))
		optvars$LB = c( 0, spec@constraints$LB)
		optvars$UB = c( 1, spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else{
		optvars$fm = m
		x0 = rep(1/m, m)
		optvars$LB = spec@constraints$LB
		optvars$UB = spec@constraints$UB
		optvars$alpha = 0.05
	}
	optvars$budget = spec@constraints$budget
	optvars$leverage = spec@constraints$leverage

	# In the minRisk problem it is possible to have no inequalities in which case
	# we must tell the solver (no target)
	optvars$ineqfun  = spec@constraints$ineqfun
	optvars$ineqgrad = spec@constraints$ineqgrad
	optvars$eqfun    = spec@constraints$eqfun
	optvars$eqgrad   = spec@constraints$eqgrad

	optvars$x0 = x0
	return(optvars)
}

.spec2optNLP = function(spec, pcontrol = list(ubounds = 1e4, 
				mbounds = 1e5, penalty = 1e4, startp = 150)){
	optvars = list()
	optvars$widx = spec@model$widx
	optvars$midx = spec@model$midx
	optvars$vidx = spec@model$vidx
	optvars$index = spec@model$indx
	if(spec@model$indx[2]==0){
		optvars$Data = switch(tolower(spec@model$risk),  
				"mad" = scale(spec@modeldata$scenario, scale = FALSE),
				"ev" = scale(spec@modeldata$scenario, scale = FALSE),
				"minimax" = spec@modeldata$scenario,
				"cvar" = spec@modeldata$scenario,
				"cdar" = spec@modeldata$scenario,
				"lpm"  = spec@modeldata$scenario,
				spec@modeldata$scenario)
	} else{
		optvars$Data = spec@modeldata$scenario
	}
	optvars$benchmark  = spec@modeldata$benchmark
	#optvars$mbenchmark = spec@modeldata$mbenchmark
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	
	risk = spec@model$risk
	mn = dim(optvars$Data)
	N = mn[1]
	m = mn[2]
	
	optvars$wm = m
	optvars$N  = N
	if(tolower(risk) == "lpm"){
		if(spec@model$options$threshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$threshold = 0
		} else{
			optvars$threshold = spec@model$options$threshold
		}
	} else{
		optvars$threshold = 0
	}
	if(tolower(risk) == "lpmupm"){
		if(spec@model$options$lthreshold == 999 || spec@model$options$uthreshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$lthreshold = 0
			optvars$uthreshold = 0
		} else{
			optvars$lthreshold = spec@model$options$lthreshold
			optvars$uthreshold = spec@model$options$uthreshold
			if(optvars$lthreshold>optvars$uthreshold) stop("\nparma: lower threshold higher than upper threshold in LPMUPM!")
		}
		optvars$lmoment = spec@model$options$lmoment
		optvars$umoment = spec@model$options$umoment
	} else{
		optvars$lthreshold = 0
		optvars$uthreshold = 0
		optvars$umoment = 1
		optvars$lmoment = 1
	}
	optvars$moment = spec@model$options$moment
	# fLB and fUB are the unconstrained upper/lower bounds for the fractional
	# problem
	if(tolower(risk) == "cvar"){
		optvars$fm = m+2
		x0 = as.numeric(c(5*quantile(optvars$Data %*% rep(1/m, m), spec@model$options$alpha), 5*rep(1/m, m), 5))
		optvars$fLB = c( -pcontrol$ubounds, rep(-pcontrol$ubounds, m), 1e-8)
		optvars$fUB = c(  pcontrol$ubounds, rep( pcontrol$ubounds, m), pcontrol$mbounds)
		#xidx = which(spec@constraints$LB<0)
		#if(length(xidx>0)) optvars$fLB[xidx] = -pcontrol$ubounds
		optvars$LB = c(spec@constraints$LB)
		optvars$UB = c(spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else if(tolower(risk) == "minimax"){
		optvars$fm = m+2
		x0 = as.numeric(c(-min(optvars$Data %*% rep(1/m, m)), rep(1/m, m), 1))
		optvars$fLB = c(-pcontrol$ubounds, rep(-pcontrol$ubounds, m), 1e-8)
		optvars$fUB = c( pcontrol$ubounds, rep( pcontrol$ubounds, m), pcontrol$mbounds)
		optvars$LB = c(spec@constraints$LB)
		optvars$UB = c(spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else{
		optvars$fm = m+1
		x0  = c(rep(1/m, m), 2)
		optvars$fLB = c(rep(-pcontrol$ubounds, m), 1e-8)
		optvars$fUB = c(rep( pcontrol$ubounds, m), pcontrol$mbounds)
		#xidx = which(spec@constraints$LB<0)
		#if(length(xidx>0)) optvars$fLB[xidx] = -10000
		optvars$LB = spec@constraints$LB
		optvars$UB = spec@constraints$UB
		optvars$alpha =  0.05
	}
	
	optvars$budget = spec@constraints$budget
	optvars$leverage = spec@constraints$leverage
	
	# In the minRisk problem it is possible to have no inequalities in which case
	# we must tell the solver (no target)
	optvars$ineqfun  = spec@constraints$ineqfun
	optvars$ineqgrad = spec@constraints$ineqgrad
	optvars$eqfun    = spec@constraints$eqfun
	optvars$eqgrad   = spec@constraints$eqgrad
	optvars$x0 = x0
	return(optvars)
}


.spec2minGNLP = function(spec, pcontrol = list(ubounds = 1e4, 
				mbounds = 1e5, penalty = 1e4, startp = 150)){
	optvars = list()
	optvars$widx = spec@model$widx
	optvars$midx = spec@model$midx
	optvars$vidx = spec@model$vidx
	optvars$index = spec@model$indx
	if(spec@model$indx[2]==0){
		optvars$Data = switch(tolower(spec@model$risk), 
				"mad" = scale(spec@modeldata$scenario, scale = FALSE),
				"ev" = scale(spec@modeldata$scenario, scale = FALSE),
				"minimax" = spec@modeldata$scenario,
				"cvar" = spec@modeldata$scenario,
				"cdar" = spec@modeldata$scenario,
				"lpm"  = spec@modeldata$scenario,
				spec@modeldata$scenario)
	} else{
		optvars$Data = spec@modeldata$scenario
	}
	optvars$benchmark  = spec@modeldata$benchmark
	#optvars$mbenchmark = spec@modeldata$mbenchmark
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	optvars$penalty  = pcontrol$penalty
	risk = spec@model$risk
	mn = dim(optvars$Data)
	N = mn[1]
	m = mn[2]
	
	optvars$wm = m
	optvars$N  = N
	if(tolower(risk) == "lpm"){
		if(spec@model$options$threshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$threshold = 0
		} else{
			optvars$threshold = spec@model$options$threshold
		}
	} else{
		optvars$threshold = 0
	}
	optvars$moment = spec@model$options$moment
	
	if(tolower(risk) == "cvar"){
		optvars$fm = m+1
		x0 = as.numeric(c(quantile(optvars$Data %*% rep(1/m, m), spec@model$options$alpha), rep(1/m, m)))
		optvars$LB = c(-10, spec@constraints$LB)
		optvars$UB = c(  0, spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else if(tolower(risk) == "minimax"){
		optvars$fm = m+1
		x0 = as.numeric(c(-min(optvars$Data %*% rep(1/m, m)), rep(1/m, m)))
		optvars$LB = c(0, spec@constraints$LB)
		optvars$UB = c(1, spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else{
		optvars$fm = m
		x0 = rep(1/m, m)
		optvars$LB = spec@constraints$LB
		optvars$UB = spec@constraints$UB
		optvars$alpha = 0.05
	}
	optvars$budget = spec@constraints$budget
	optvars$leverage = spec@constraints$leverage
	
	# In the minRisk problem it is possible to have no inequalities in which case
	# we must tell the solver (no target)
	optvars$ineqfun  = spec@constraints$ineqfun
	optvars$ineqgrad = spec@constraints$ineqgrad
	optvars$eqfun    = spec@constraints$eqfun
	optvars$eqgrad   = spec@constraints$eqgrad
	optvars$x0 = x0
	return(optvars)
}

.spec2optGNLP = function(spec, pcontrol = list(ubounds = 1e4, 
				mbounds = 1e5, penalty = 1e4, startp = 150)){
	optvars = list()
	optvars$widx = spec@model$widx
	optvars$midx = spec@model$midx
	optvars$vidx = spec@model$vidx
	optvars$index = spec@model$indx
	if(spec@model$indx[2]==0){
		optvars$Data = switch(tolower(spec@model$risk), 
				"mad" = scale(spec@modeldata$scenario, scale = FALSE),
				"ev" = scale(spec@modeldata$scenario, scale = FALSE),
				"minimax" = spec@modeldata$scenario,
				"cvar" = spec@modeldata$scenario,
				"cdar" = spec@modeldata$scenario,
				"lpm"  = spec@modeldata$scenario,
				spec@modeldata$scenario)
	} else{
		optvars$Data = spec@modeldata$scenario
	}
	optvars$benchmark  = spec@modeldata$benchmark
	#optvars$mbenchmark = spec@modeldata$mbenchmark
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	penalty = pcontrol$penalty
	
	risk = spec@model$risk
	mn = dim(optvars$Data)
	N = mn[1]
	m = mn[2]
	
	optvars$wm = m
	optvars$N  = N
	if(tolower(risk) == "lpm"){
		if(spec@model$options$threshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$threshold = 0
		} else{
			optvars$threshold = spec@model$options$threshold
		}
	} else{
		optvars$threshold = 0
	}
	if(tolower(risk) == "lpmupm"){
		if(spec@model$options$lthreshold == 999 || spec@model$options$uthreshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$lthreshold = 0
			optvars$uthreshold = 0
		} else{
			optvars$lthreshold = spec@model$options$lthreshold
			optvars$uthreshold = spec@model$options$uthreshold
			if(optvars$lthreshold>optvars$uthreshold) stop("\nparma: lower threshold higher than upper threshold in LPMUPM!")
		}
		optvars$lmoment = spec@model$options$lmoment
		optvars$umoment = spec@model$options$umoment
	} else{
		optvars$lthreshold = 0
		optvars$uthreshold = 0
		optvars$umoment = 1
		optvars$lmoment = 1
	}
	optvars$moment = spec@model$options$moment
	# fLB and fUB are the unconstrained upper/lower bounds for the fractional
	# problem
	if(tolower(risk) == "cvar"){
		optvars$fm = m+2
		x0 = as.numeric(c(5*quantile(optvars$Data %*% rep(1/m, m), spec@model$options$alpha), 5*rep(1/m, m), 5))
		optvars$fLB = c( -pcontrol$ubounds, rep(-pcontrol$ubounds, m), 1e-8)
		optvars$fUB = c(  pcontrol$ubounds, rep( pcontrol$ubounds, m), pcontrol$mbounds)
		#xidx = which(spec@constraints$LB<0)
		#if(length(xidx>0)) optvars$fLB[xidx] = -pcontrol$ubounds
		optvars$LB = c(spec@constraints$LB)
		optvars$UB = c(spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else if(tolower(risk) == "minimax"){
		optvars$fm = m+2
		x0 = as.numeric(c(-min(optvars$Data %*% rep(1/m, m)), rep(1/m, m), 1))
		optvars$fLB = c(-pcontrol$ubounds, rep(-pcontrol$ubounds, m), 1e-8)
		optvars$fUB = c( pcontrol$ubounds, rep( pcontrol$ubounds, m), pcontrol$mbounds)
		optvars$LB = c(spec@constraints$LB)
		optvars$UB = c(spec@constraints$UB)
		optvars$alpha =  spec@model$options$alpha
	} else{
		optvars$fm = m+1
		x0  = c(rep(1/m, m), 2)
		optvars$fLB = c(rep(-pcontrol$ubounds, m), 1e-8)
		optvars$fUB = c(rep( pcontrol$ubounds, m), pcontrol$mbounds)
		#xidx = which(spec@constraints$LB<0)
		#if(length(xidx>0)) optvars$fLB[xidx] = -10000
		optvars$LB = spec@constraints$LB
		optvars$UB = spec@constraints$UB
		optvars$alpha =  0.05
	}
	optvars$budget = spec@constraints$budget
	optvars$leverage = spec@constraints$leverage
	optvars$ineqfun  = spec@constraints$ineqfun
	optvars$ineqgrad = spec@constraints$ineqgrad
	optvars$eqfun    = spec@constraints$eqfun
	optvars$eqgrad   = spec@constraints$eqgrad
	optvars$x0 = x0
	return(optvars)
}

.spec2QP = function(spec){
	optvars = list()
	# optvars$widx = spec@model$widx
	# optvars$midx = spec@model$midx
	optvars$index = spec@model$indx
	optvars$S = spec@modeldata$S
	optvars$benchmarkS = spec@modeldata$benchmarkS
	#optvars$benchmarkM = spec@modeldata$benchmarkM
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	optvars$budget = spec@constraints$budget
	optvars$ineq.mat  = spec@constraints$ineq.mat
	optvars$eq.mat    = spec@constraints$eq.mat
	optvars$LB = spec@constraints$LB
	optvars$UB = spec@constraints$UB
	optvars$ineq.LB  = spec@constraints$ineq.LB
	optvars$ineq.UB  = spec@constraints$ineq.UB	
	optvars$eqB      = spec@constraints$eqB	
	return(optvars)
}

.spec2SOCP = function(spec){
	optvars = list()
	# optvars$widx = spec@model$widx
	# optvars$midx = spec@model$midx
	optvars$index = spec@model$indx
	optvars$S = spec@modeldata$S
	optvars$Q = spec@modeldata$Q
	optvars$qB = spec@modeldata$qB
	optvars$riskB = spec@modeldata$riskB
	optvars$benchmarkS = spec@modeldata$benchmarkS
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	optvars$budget = spec@constraints$budget
	optvars$leverage = spec@constraints$leverage
	optvars$ineq.mat  = spec@constraints$ineq.mat
	optvars$eq.mat    = spec@constraints$eq.mat
	optvars$LB = spec@constraints$LB
	optvars$UB = spec@constraints$UB
	optvars$ineq.LB  = spec@constraints$ineq.LB
	optvars$ineq.UB  = spec@constraints$ineq.UB	
	optvars$eqB      = spec@constraints$eqB	
	return(optvars)
}

.spec2LP = function(spec){
	optvars = list()
	optvars$index = spec@model$indx
	if(spec@model$indx[2]==0){
		optvars$Data = switch(tolower(spec@model$risk), 
				"mad" = scale(spec@modeldata$scenario, scale = FALSE),
				"ev" = scale(spec@modeldata$scenario, scale = FALSE),
				"minimax" = spec@modeldata$scenario,
				"cvar" = spec@modeldata$scenario,
				"cdar" = spec@modeldata$scenario,
				"lpm"  = spec@modeldata$scenario,
				spec@modeldata$scenario)
	} else{
		optvars$Data = spec@modeldata$scenario
	}
	optvars$probability = spec@modeldata$probability
	optvars$benchmark  = spec@modeldata$benchmark
	#optvars$mbenchmark = spec@modeldata$mbenchmark
	optvars$mu = spec@modeldata$forecast
	optvars$mutarget = spec@modeldata$target
	risk = spec@model$risk
	mn = dim(optvars$Data)
	N = mn[1]
	m = mn[2]	
	optvars$wm = m
	optvars$N  = N
	if(tolower(risk) == "lpm"){
		if(spec@model$options$threshold == 999){
			optvars$Data = scale(optvars$Data, scale = FALSE)
			optvars$threshold = 0
		} else{
			optvars$threshold = spec@model$options$threshold
		}
	} else{
		optvars$threshold = 0
	}
	optvars$moment = 1
	optvars$alpha =  spec@model$options$alpha
	optvars$LB = spec@constraints$LB
	optvars$UB = spec@constraints$UB
	optvars$budget = spec@constraints$budget
	optvars$ineq.mat = spec@constraints$ineq.mat
	optvars$ineq.LB  = spec@constraints$ineq.LB
	optvars$ineq.UB  = spec@constraints$ineq.UB	
	optvars$eq.mat   = spec@constraints$eq.mat
	optvars$eqB      = spec@constraints$eqB	
	optvars$max.pos  = spec@constraints$max.pos
	return(optvars)
}

####----------------------------------------------------------------------------

.parmasolve = function(spec, type = NULL, solver = NULL, solver.control = list(), 
		x0 = NULL, w0 = NULL, parma.control = list(ubounds = 1e4, 
				mbounds = 1e5, penalty = 1e4, eqSlack = 1e-5), ...){
	tic = Sys.time()
	# pass the problem to the correct function and do some checks
	# c("data", "benchmark", "target", "risk", "risktype", "leverage", "aux1", "aux2")
	
	if(is.null(parma.control)) parma.control = list()
	
	mm = match(names(parma.control), c("ubounds", "mbounds", "penalty","eqSlack"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, parma.control[idx[i]])
		warning(paste(c("unidentified option(s) in parma.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	if(is.null(parma.control$ubounds)) parma.control$ubounds = 1e4 else parma.control$ubounds = parma.control$ubounds[1]
	if(is.null(parma.control$mbounds)) parma.control$mbounds = 1e5 else parma.control$mbounds = parma.control$mbounds[1]
	if(is.null(parma.control$penalty)) parma.control$penalty = 1e4 else parma.control$penalty = parma.control$penalty[1]
	if(is.null(parma.control$eqSlack)) parma.control$eqSlack = 1e-5 else parma.control$eqSlack = parma.control$eqSlack[1]
	
	available.problems = toupper(c("LP", "MILP", "QP", "MIQP", "SOCP", "NLP", "MINLP", "GNLP")[which(spec@model$type==1)])
	
	if(!is.null(type[1])){
		type = toupper(type[1])
		tmp = match.arg(type, available.problems)
		type = tmp
	} else{
		type = available.problems[1]
	}
	if(spec@model$indx[5]==1){
		optvars = switch(tolower(type),
				"lp"   = .spec2LP(spec), 
				"milp" = .spec2LP(spec),
				"nlp"  = .spec2minNLP(spec),
				"qp"   = .spec2QP(spec),
				"socp" = .spec2SOCP(spec),
				"gnlp" = .spec2minGNLP(spec, parma.control))
	} else{
		optvars = switch(tolower(type),
				"lp"   = .spec2LP(spec), 
				"milp" = .spec2LP(spec),
				"nlp"  = .spec2optNLP(spec, parma.control),
				"qp"   = .spec2QP(spec),
				"socp" = .spec2SOCP(spec),
				"gnlp" = .spec2optGNLP(spec, parma.control))
	}
	uservars = spec@modeldata$uservars	
	if(!is.null(x0)){
		if(length(optvars$x0)!=length(x0)) stop("\nparma: wrong length for x0!")
		optvars$x0 = x0
	}
	if(!is.null(w0)){
		if(length(optvars$widx)!=length(w0)) stop("\nparma: wrong length for w0!")
		optvars$x0[optvars$widx] = w0
	}
	if(type=="LP" && is.null(solver)) solver="GLPK"
	sol = switch(toupper(type),
			"LP"   = lpport(optvars, solver, ...),
			"MILP" = milpport(optvars, solver, ...),
			"NLP"  = nlpport(optvars, uservars, control = solver.control, ...),
			"QP"   = qpport(optvars,   ...),
			"SOCP" = socpport(optvars, control = solver.control, eqSlack = parma.control$eqSlack, ...),
			"GNLP" = gnlpport(optvars, uservars, solver = solver, control = solver.control, ...))
	# arbitrage check
	if(type!="QP" & type!="SOCP"){
		arbitrage = .arbcheck(sol$weights, spec@modeldata$scenario, spec@model$options, spec@model$risk)
	} else{
		arbitrage = c(0, 0)
	}
	spec@model$asset.names = spec@modeldata$asset.names
	spec@model$type = type
	sol$solver = solver
	sol$arbitrage = arbitrage
	toc = Sys.time() - tic
	spec@model$elapsed = toc
	ret = new("parmaPort", 
			solution = sol, 
			model = spec@model)
	return(ret)
}
################################################################################
.parmafrontier = function(spec, n.points = 100, miny = NULL, maxy = NULL, 
		type = NULL, solver = NULL, solver.control = list(), 
		parma.control = list(ubounds = 10000, mbounds = 1e+05, penalty = 10000), 
		cluster = NULL)
{
	if(!is.null(spec@modeldata$S))
	{
		ans = m.parmafrontier(spec, n.points = n.points, type = type, 
				miny = miny, maxy = maxy, cluster = cluster)
	} else {
		ans = s.parmafrontier(spec, n.points = n.points, miny = miny, 
				maxy = maxy, type = type, solver = solver, 
				solver.control = solver.control, 
				parma.control =  parma.control, cluster = cluster)		
	}
	return(ans)
}

s.parmafrontier = function(spec, n.points = 100, miny = NULL, 
		maxy = NULL, type = NULL, solver = NULL, solver.control = list(), 
		parma.control = list(ubounds = 10000, mbounds = 1e+05, penalty = 10000), 
		cluster = NULL)
{
	targettype = parmaget(spec, "targetType")
	risktype = parmaget(spec, "riskType")
	if(risktype!="minrisk") stop("\nspec riskType must be minrisk...fix and resubmit.")
	if(targettype!="equality") stop("\nspec targetType must be equality...fix and resubmit.")
	m = NCOL(spec@modeldata$scenario)
	if(is.null(spec@modeldata$forecast)){
		f = abs(apply(spec@modeldata$scenario, 2, "mean"))
		minb = min(f)
		maxb = max(f)
	} else{
		minb = min(abs(spec@modeldata$forecast))
		maxb = max(abs(spec@modeldata$forecast))
	}
	if(is.null(miny)){
		xspec = spec
		parmaset(xspec)<-list(target=0)
		parmaset(xspec)<-list(targetType="inequality")
		solx = try(parmasolve(xspec, type = type, solver = solver, 
							solver.control = solver.control, 
							parma.control = parma.control), silent = TRUE)
		if(!inherits(solx, "try-error")) minb = parmareward(solx)
	} else{
		minb = miny
	}
	if(!is.null(maxy)){
		maxb = maxy
	}
	fs = seq(minb, maxb, length.out = n.points)
	fmat = matrix(NA, ncol = m+3, nrow = n.points)
	if(!is.null(cluster)){
		clusterEvalQ(cluster, require(parma))
		clusterExport(cluster, c("fs", "spec", "type", "solver", 
						"solver.control","parma.control"), envir = environment())
		sol = parLapply(cluster, 1:n.points, fun = function(i){
					xspec = spec
					parmaset(xspec)<-list(target=fs[i])
					tmp = parmasolve(xspec, type = type, solver = solver, 
							solver.control = solver.control, 
							parma.control = parma.control)
					return(tmp)
				})
		for(i in 1:n.points){
			fmat[i,1:m] = weights(sol[[i]])
			fmat[i,m+1] = parmarisk(sol[[i]])
			fmat[i,m+2] = parmareward(sol[[i]])
			fmat[i,m+3] = parmastatus(sol[[i]])
		}
	} else{
		for(i in 1:n.points){
			xspec = spec
			parmaset(xspec)<-list(target=fs[i])
			tmp = parmasolve(xspec, type = type, solver = solver, 
					solver.control = solver.control, 
					parma.control = parma.control)
			fmat[i,1:m] = weights(tmp)
			fmat[i,m+1] = parmarisk(tmp)
			fmat[i,m+2] = parmareward(tmp)
			fmat[i,m+3] = parmastatus(tmp)
		}
	}
	colnames(fmat) = c(spec@modeldata$asset.names, spec@model$risk, "reward", "status")
	return(fmat)
}

m.parmafrontier = function(spec, n.points = 100, type = "QP", 
		solver.control = list(abs.tol = 1e-8, rel.tol = 1e-8, Nu=2, max.iter=5250, 
				BigM.K = 4, BigM.iter = 15), miny = NULL, maxy = NULL, cluster = NULL)
{
	targettype = parmaget(spec, "targetType")
	risktype = parmaget(spec, "riskType")
	if(risktype!="minrisk") stop("\nspec riskType must be minrisk...fix and resubmit.")
	if(targettype!="equality") stop("\nspec targetType must be equality...fix and resubmit.")
	m = NCOL(spec@modeldata$S)
	if(is.null(spec@modeldata$forecast)){
		stop("\nparma: cannot have a NULL forecast vector in QP formulation.")
	} else{
		minb = min(abs(spec@modeldata$forecast))
		maxb = max(abs(spec@modeldata$forecast))
	}
	if(is.null(miny)){
		xspec = spec
		parmaset(xspec)<-list(target=0)
		parmaset(xspec)<-list(targetType="inequality")
		solx = try(parmasolve(xspec), silent = TRUE)
		if(!inherits(solx, "try-error")) minb = parmareward(solx)
	} else{
		minb = miny
	}
	if(!is.null(maxy)){
		maxb = maxy
	}
	fs = seq(minb, maxb, length.out = n.points)
	fmat = matrix(NA, ncol = m+2, nrow = n.points)
	if(!is.null(cluster)){
		clusterEvalQ(cluster, require(parma))
		clusterExport(cluster, c("fs", "spec","type","solver.control"), envir = environment())
		sol = parLapply(cluster, 1:n.points, fun = function(i){
					xspec = spec
					parmaset(xspec)<-list(target=fs[i])
					tmp = parmasolve(xspec, type = type, solver.control = solver.control)
					return(tmp)
				})
		for(i in 1:n.points){
			fmat[i,1:m] = weights(sol[[i]])
			fmat[i,m+1] = parmarisk(sol[[i]])
			fmat[i,m+2] = parmareward(sol[[i]])
		}
	} else{
		for(i in 1:n.points){
			xspec = spec
			parmaset(xspec)<-list(target=fs[i])
			tmp = parmasolve(xspec, type = type, solver.control = solver.control)
			fmat[i,1:m] = weights(tmp)
			fmat[i,m+1] = parmarisk(tmp)
			fmat[i,m+2] = parmareward(tmp)
		}
	}
	colnames(fmat) = c(spec@modeldata$asset.names, spec@model$risk, "reward")
	return(fmat)
}

.checkconsfun = function(fun, name = "ineqfun"){
	if(!is.list(fun)) stop(paste("\n",name," must be a list of functions",sep=""))
	n = length(fun)
	for(i in 1:n){
		if(!is.function(fun[[i]])) stop(paste("\n",name," list constains non functions at position: ",i,sep=""))
	}
	return(0)
}