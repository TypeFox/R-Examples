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

# Risk Functions
parma.test1 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	set.seed(20)
	w = rexp(15, 2)
	w = w/sum(w)
	a = -0.2
	b = 2
	pmatrix = matrix(0, ncol = 8, nrow = 4)
	colnames(pmatrix) = c("mad", "ev", "sd", "minimax", "cvar", "cdar", "lpm[m,t=c]", "lpm[m,t=mu]")
	rownames(pmatrix) = c("Scaling", "Location[1]", "Location[2]", "Subadditivity")
	#1. MAD and properties
	m1 = riskfun(w, Data, risk = "mad")
	# scaling property: YES
	pmatrix[1,1] = as.integer(riskfun(w, b*Data, risk = "mad")/abs(b) == m1)
	# location invariant: YES
	pmatrix[2,1] = as.integer(riskfun(w, a+Data, risk = "mad") == m1)
	# Subadditivity: Yes
	pmatrix[4,1] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "mad") + 
					riskfun(w[6:15], Data[,6:15], risk = "mad")) >= m1)
		
	#2. EV and properties
	m2 = riskfun(w[1:15], Data, risk = "ev")
	# scaling property: YES
	pmatrix[1,2] = as.integer(riskfun(w, b*Data, risk = "ev")/b^2 == m2)
	# location invariant: YES
	pmatrix[2,2] = as.integer(riskfun(w, a+Data, risk = "ev") == m2)
	# Subadditivity: NO (Variance is Superadditive)
	pmatrix[4,2] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "ev") + 
						riskfun(w[6:15], Data[,6:15], risk = "ev")) >= m2)
	
	#3. SD and properties
	m3 = sqrt(riskfun(w, Data, risk = "ev"))
	# scaling property: YES
	pmatrix[1,3] = as.integer(sqrt(riskfun(w, b*Data, risk = "ev"))/abs(b) == m3)
	# location invariant: YES
	pmatrix[2,3] = as.integer(sqrt(riskfun(w, a+Data, risk = "ev")) == m3)
	# Subadditivity: YES
	pmatrix[4,3] = as.integer( (sqrt(riskfun(w[1:5], Data[,1:5], risk = "ev")) + 
						sqrt(riskfun(w[6:15], Data[,6:15], risk = "ev"))) >= m3)
	
	#4. Minimax and properties
	m4 = riskfun(w, Data, risk = "minimax")
	# scaling property: YES*
	pmatrix[1,4] = as.integer(riskfun(w, abs(b)*Data, risk = "minimax")/abs(b) == m4)
	# location invariant: NO[1]/YES[2]*
	pmatrix[2,4] = as.integer(riskfun(w, a+Data, risk = "minimax") == m4)
	pmatrix[3,4] = as.integer(round(riskfun(w, a+Data, risk = "minimax")+a,12) == round(m4,12))
	# Subadditivity: YES
	pmatrix[4,4] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "minimax") + 
						riskfun(w[6:15], Data[,6:15], risk = "minimax")) >= m4)
	
	#5. CVaR and properties
	m5 = riskfun(w, Data, risk = "cvar", alpha = 0.05)
	# scaling property: YES*
	pmatrix[1,5] = as.integer(round(riskfun(w, abs(b)*Data, risk = "cvar", alpha = 0.05)/abs(b),12) == round(m5,12))
	# location invariant: NO[1]/YES[2]*
	pmatrix[2,5] = as.integer(round(riskfun(w, a+Data, risk = "cvar", alpha = 0.05),12) == round(m5,12))
	pmatrix[3,5] = as.integer(round(riskfun(w, a+Data, risk = "cvar", alpha = 0.05)+a,12) == round(m5,12))
	# Subadditivity: YES
	pmatrix[4,5] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "cvar", alpha = 0.05) + 
						riskfun(w[6:15], Data[,6:15], risk = "cvar", alpha = 0.05)) >= m5)
	
	#6. CDaR and properties
	m6 = riskfun(w, Data, risk = "cdar", alpha = 0.05)
	# scaling property: YES*
	pmatrix[1,6] = as.integer(round(riskfun(w, abs(b)*Data, risk = "cdar", alpha = 0.05)/abs(b),12) == round(m6,12))
	# location invariant: NO[1]/NO[2]
	pmatrix[2,6] = as.integer(round(riskfun(w, a+Data, risk = "cdar", alpha = 0.05),12) == round(m6,12))
	pmatrix[3,6] = as.integer(round(riskfun(w, a+Data, risk = "cdar", alpha = 0.05)+a,12) == round(m6,12))
	# Subadditivity: YES
	pmatrix[4,6] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "cdar", alpha = 0.05) + 
						riskfun(w[6:15], Data[,6:15], risk = "cdar", alpha = 0.05)) >= m6)
	
	#7. LPM[m, t=C] and properties
	m7 = riskfun(w, Data, risk = "lpm",  moment = 1, threshold = 0.02)
	# scaling property: YES*
	pmatrix[1,7] = as.integer(round(riskfun(w, abs(b)*Data, risk = "lpm", moment = 1, threshold = abs(b)*0.02)/abs(b),12) == round(m7,12))
	# location invariant: YES[1]/NO[2]
	pmatrix[3,7] = as.integer(round(riskfun(w, a+Data, risk = "lpm",  moment = 1, threshold = a+0.02),12) == round(m7,12))
	# Subadditivity: NO
	pmatrix[4,7] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "lpm",  moment = 1, threshold = 0.02) + 
						riskfun(w[6:15], Data[,6:15], risk = "lpm",  moment = 1, threshold = 0.02)) >= m7)
	
	#8. LPM[m, t=mean] and properties
	m8 = riskfun(w, Data, risk = "lpm",  moment = 2, threshold = 999)
	# equivalent to: riskfun(w, scale(Data, scale=FALSE), risk = "lpm",  moment = 2, threshold = 0)
	# scaling property: YES*
	pmatrix[1,8] = as.integer(round(riskfun(w, abs(b)*Data, risk = "lpm", moment = 2, threshold = 999)/abs(b),12) == round(m8,12))	
	# location invariant: YES[1]/NO[2]
	pmatrix[2,8] = as.integer(round(riskfun(w, a+Data, risk = "lpm",  moment = 2, threshold = 999),12) == round(m8,12))
	# Subadditivity: YES
	pmatrix[4,8] = as.integer( (riskfun(w[1:5], Data[,1:5], risk = "lpm",  moment = 2, threshold = 999) + 
						riskfun(w[6:15], Data[,6:15], risk = "lpm",  moment = 2, threshold = 999)) >= m8)
	
	
	zz <- file("parma_test1-1.txt", open="wt")
	sink(zz)
	pm = as.data.frame(pmatrix)
	pm[pm==1] = "TRUE"
	pm[pm==0] = "FALSE"
	print(pm)
	sink(type="message")
	sink()
	close(zz)
	
	postscript("parma_test1-1.eps")
	pp = ecdf(sort(Data %*% w))
	plot(sort(Data %*% w), pp(sort(Data %*% w)), type = "l", col = "black", ylab = "P", xlab="x",
			main = "CDF")
	abline(v = mean(Data %*% w), col = "grey", lty = 2)
	abline(h = pp(mean(Data %*% w)), col = "grey", lty=2)
	points(-m1, pp(-m1), col = "green", pch=21, bg = "green")
	points(-m3, pp(-m3), col = "orange", pch=21, bg = "orange")
	points(-m4, pp(-m4), col = "yellow", pch=21, bg = "yellow")
	points(-m5, pp(-m5), col = "tomato1", pch=21, bg = "tomato1")
	points(-m8, pp(-m8), col = "steelblue", pch=21, bg = "steelblue")
	legend("topleft", c("MAD","SD", "Min", "CVaR[5%]", "LPM[m=2,t=mu]"), 
			col = c("green", "orange", "yellow", "tomato1", "steelblue"),
			pch = rep(21, 5), 
			pt.bg = c("green", "orange", "yellow", "tomato1", "steelblue"), bty="n")
	dev.off()
	
	toc = Sys.time() - tic
	return(toc)
	
}


# MinRisk LP vs NLP vs QP
parma.test2 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	
	pspec1 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[1], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol1lpmin = parmasolve(pspec1, type = "LP")
	sol1nlpmin = parmasolve(pspec1, type = "NLP")
	
	d1 = cbind(weights(sol1lpmin), weights(sol1nlpmin))
	d1 = rbind(d1, unname(cbind(parmarisk(sol1lpmin), parmarisk(sol1nlpmin))))
	d1 = rbind(d1, cbind(parmareward(sol1lpmin), parmareward(sol1nlpmin)))
	d1 = rbind(d1, cbind(as.numeric(tictoc(sol1lpmin)), as.numeric(tictoc(sol1nlpmin))))
	rownames(d1)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d1) = c("LP", "NLP")
	
	options(width = 120)
	zz <- file("parma_test2-1.txt", open="wt")
	sink(zz)	
	show(pspec1)
	print(round(d1,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec2 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[2], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol2lpmin = parmasolve(pspec2, type = "LP")
	sol2nlpmin = parmasolve(pspec2, type = "NLP")
	
	d2 = cbind(weights(sol2lpmin), weights(sol2nlpmin))
	d2 = rbind(d2, unname(cbind(parmarisk(sol2lpmin), parmarisk(sol2nlpmin))))
	d2 = rbind(d2, cbind(parmareward(sol2lpmin), parmareward(sol2nlpmin)))
	d2 = rbind(d2, cbind(as.numeric(tictoc(sol2lpmin)), as.numeric(tictoc(sol2nlpmin))))
	rownames(d2)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d2) = c("LP", "NLP")
	
	zz <- file("parma_test2-2.txt", open="wt")
	sink(zz)	
	show(pspec2)
	print(round(d2,6))
	sink(type="message")
	sink()
	close(zz)
	
	pspec3 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			options = list(alpha = 0.05), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol3lpmin = parmasolve(pspec3, type = "LP")
	sol3nlpmin = parmasolve(pspec3, type = "NLP")
	
	d3 = cbind(weights(sol3lpmin), weights(sol3nlpmin))
	d3 = rbind(d3, unname(cbind(parmarisk(sol3lpmin), parmarisk(sol3nlpmin))))
	d3 = rbind(d3, cbind(parmareward(sol3lpmin), parmareward(sol3nlpmin)))
	d3 = rbind(d3, cbind(sol3lpmin@solution$VaR, sol3nlpmin@solution$VaR))	
	d3 = rbind(d3, cbind(as.numeric(tictoc(sol3lpmin)), as.numeric(tictoc(sol3nlpmin))))
	rownames(d3)[16:19] = c("risk", "reward", "VaR", "elapsed(secs)")
	colnames(d3) = c("LP", "NLP")
	
	zz <- file("parma_test2-3.txt", open="wt")
	sink(zz)	
	show(pspec3)
	print(round(d3,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec4 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[4], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			options = list(alpha = 0.05), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol4lpmin = parmasolve(pspec4, type = "LP")
	
	d4 = cbind(weights(sol4lpmin))
	d4 = rbind(d4, unname(cbind(parmarisk(sol4lpmin))))
	d4 = rbind(d4, cbind(parmareward(sol4lpmin)))
	d4 = rbind(d4, cbind(sol4lpmin@solution$DaR))	
	d4 = rbind(d4, cbind(as.numeric(tictoc(sol4lpmin))))
	rownames(d4)[16:19] = c("risk", "reward", "DaR", "elapsed(secs)")
	colnames(d4) = c("LP")
	
	zz <- file("parma_test2-4.txt", open="wt")
	sink(zz)	
	show(pspec4)
	print(round(d4,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec5 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	pspec5q = parmaspec(S = cov(Data), forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol5nlpmin = parmasolve(pspec5, type = "NLP")
	sol5qpmin = parmasolve(pspec5q, type = "QP")
	
	d5 = cbind(weights(sol5qpmin), weights(sol5nlpmin))
	d5 = rbind(d5, sqrt(unname(cbind(parmarisk(sol5qpmin), parmarisk(sol5nlpmin)))))
	d5 = rbind(d5, cbind(parmareward(sol5qpmin), parmareward(sol5nlpmin)))
	d5 = rbind(d5, cbind(as.numeric(tictoc(sol5qpmin)), as.numeric(tictoc(sol5nlpmin))))
	rownames(d5)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d5) = c("QP", "NLP")
	
	zz <- file("parma_test2-5.txt", open="wt")
	sink(zz)	
	show(pspec5)
	show(pspec5q)
	print(round(d5,6))
	sink(type="message")
	sink()
	close(zz)
	
	pspec6 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[6], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[1], 
			options = list(threshold = 999, moment=1), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol6lpmin = parmasolve(pspec6, type = "LP")
	sol6nlpmin = parmasolve(pspec6, type = "NLP")
	
	d6 = cbind(weights(sol6lpmin), weights(sol6nlpmin))
	d6 = rbind(d6, unname(cbind(parmarisk(sol6lpmin), parmarisk(sol6nlpmin))))
	d6 = rbind(d6, cbind(parmareward(sol6lpmin), parmareward(sol6nlpmin)))
	d6 = rbind(d6, cbind(as.numeric(tictoc(sol6lpmin)), as.numeric(tictoc(sol6nlpmin))))
	rownames(d6)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d6) = c("LP", "NLP")
	
	zz <- file("parma_test2-6.txt", open="wt")
	sink(zz)	
	show(pspec6)
	print(round(d6,6))
	sink(type="message")
	sink()
	close(zz)
	
	# Mean Squared Error LP-NLP, QP-NLP
	errmat = matrix(NA, ncol = 5, nrow = 4)
	colnames(errmat) = c("MAD", "MiniMax", "CVaR", "EV", "LPM[1]")
	rownames(errmat) = c("MSE[weights]", "MAE[weights]", "MaxE[weights]", "Err[risk]")
	errmat[1,1] = mean( (weights(sol1lpmin)-weights(sol1nlpmin))^2)
	errmat[2,1] = mean( abs(weights(sol1lpmin)-weights(sol1nlpmin)))
	errmat[3,1] = max( abs(weights(sol1lpmin)-weights(sol1nlpmin)))
	errmat[4,1] = parmarisk(sol1nlpmin)-parmarisk(sol1lpmin)
	
	errmat[1,2] = mean( (weights(sol2lpmin)-weights(sol2nlpmin))^2)
	errmat[2,2] = mean( abs(weights(sol2lpmin)-weights(sol2nlpmin)))
	errmat[3,2] = max( abs(weights(sol2lpmin)-weights(sol2nlpmin)))
	errmat[4,2] = parmarisk(sol2lpmin)-parmarisk(sol2nlpmin)
	
	errmat[1,3] = mean( (weights(sol3lpmin)-weights(sol3nlpmin))^2)
	errmat[2,3] = mean( abs(weights(sol3lpmin)-weights(sol3nlpmin)))
	errmat[3,3] = max( abs(weights(sol3lpmin)-weights(sol3nlpmin)))
	errmat[4,3] = parmarisk(sol3lpmin)-parmarisk(sol3nlpmin)
	
	errmat[1,4] = mean( (weights(sol5qpmin)-weights(sol5nlpmin))^2)
	errmat[2,4] = mean( abs(weights(sol5qpmin)-weights(sol5nlpmin)))
	errmat[3,4] = max( abs(weights(sol5qpmin)-weights(sol5nlpmin)))
	errmat[4,4] = parmarisk(sol5qpmin)-parmarisk(sol5nlpmin)
	
	errmat[1,5] = mean( (weights(sol6lpmin)-weights(sol6nlpmin))^2)
	errmat[2,5] = mean( abs(weights(sol6lpmin)-weights(sol6nlpmin)))
	errmat[3,5] = max( abs(weights(sol6lpmin)-weights(sol6nlpmin)))
	errmat[4,5] = parmarisk(sol6lpmin)-parmarisk(sol6nlpmin)
	
	zz <- file("parma_test2-7.txt", open="wt")
	sink(zz)	
	print(errmat)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
}

# OptRisk LP vs NLP vs QP
parma.test3 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	
	pspec1 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[1], 
			target =NULL, targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol1lpopt = parmasolve(pspec1, type = "LP", solver="GLPK")
	sol1nlpopt = parmasolve(pspec1, type = "NLP")
	
	d1 = cbind(weights(sol1lpopt), weights(sol1nlpopt))
	d1 = rbind(d1, unname(cbind(parmarisk(sol1lpopt), parmarisk(sol1nlpopt))))
	d1 = rbind(d1, cbind(parmareward(sol1lpopt), parmareward(sol1nlpopt)))
	d1 = rbind(d1, cbind(parmarisk(sol1lpopt)/parmareward(sol1lpopt), parmarisk(sol1nlpopt)/parmareward(sol1nlpopt)))
	d1 = rbind(d1, cbind(as.numeric(tictoc(sol1lpopt)), as.numeric(tictoc(sol1nlpopt))))
	rownames(d1)[16:19] = c("risk", "reward", "risk/reward", "elapsed(secs)")
	colnames(d1) = c("LP", "NLP")
	
	options(width = 120)
	zz <- file("parma_test3-1.txt", open="wt")
	sink(zz)	
	show(pspec1)
	print(round(d1,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec2 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[2], 
			target = NULL, targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			LB = rep(0, 15), UB = rep(0.5, 15), budget = 1)
	
	sol2lpopt = parmasolve(pspec2, type = "LP", solver="GLPK")
	sol2nlpopt = parmasolve(pspec2, type = "NLP")
	
	d2 = cbind(weights(sol2lpopt), weights(sol2nlpopt))
	d2 = rbind(d2, unname(cbind(parmarisk(sol2lpopt), parmarisk(sol2nlpopt))))
	d2 = rbind(d2, cbind(parmareward(sol2lpopt), parmareward(sol2nlpopt)))
	d2 = rbind(d2, cbind(as.numeric(tictoc(sol2lpopt)), as.numeric(tictoc(sol2nlpopt))))
	rownames(d2)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d2) = c("LP", "NLP")
	
	zz <- file("parma_test3-2.txt", open="wt")
	sink(zz)	
	show(pspec2)
	print(round(d2,6))
	sink(type="message")
	sink()
	close(zz)
	
	pspec3 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol3lpopt = parmasolve(pspec3, type = "LP", solver="GLPK")
	sol3nlpopt = parmasolve(pspec3, type = "NLP")
	
	d3 = cbind(weights(sol3lpopt), weights(sol3nlpopt))
	d3 = rbind(d3, unname(cbind(parmarisk(sol3lpopt), parmarisk(sol3nlpopt))))
	d3 = rbind(d3, cbind(parmareward(sol3lpopt), parmareward(sol3nlpopt)))
	d3 = rbind(d3, cbind(sol3lpopt@solution$VaR, sol3nlpopt@solution$VaR))	
	d3 = rbind(d3, cbind(as.numeric(tictoc(sol3lpopt)), as.numeric(tictoc(sol3nlpopt))))
	rownames(d3)[16:19] = c("risk", "reward", "VaR", "elapsed(secs)")
	colnames(d3) = c("LP", "NLP")
	
	zz <- file("parma_test3-3.txt", open="wt")
	sink(zz)	
	show(pspec3)
	print(round(d3,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec4 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[4], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol4lpopt = parmasolve(pspec4, type = "LP", solver="GLPK")
	
	d4 = cbind(weights(sol4lpopt))
	d4 = rbind(d4, unname(cbind(parmarisk(sol4lpopt))))
	d4 = rbind(d4, cbind(parmareward(sol4lpopt)))
	d4 = rbind(d4, cbind(sol4lpopt@solution$DaR))	
	d4 = rbind(d4, cbind(as.numeric(tictoc(sol4lpopt))))
	rownames(d4)[16:19] = c("risk", "reward", "DaR", "elapsed(secs)")
	colnames(d4) = c("LP")
	
	zz <- file("parma_test3-4.txt", open="wt")
	sink(zz)	
	show(pspec4)
	print(round(d4,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec5 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	pspec5q = parmaspec(S = cov(Data), forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol5nlpopt = parmasolve(pspec5, type = "NLP")
	sol5qpopt = parmasolve(pspec5q, type = "QP")
	
	d5 = cbind(weights(sol5qpopt), weights(sol5nlpopt))
	d5 = rbind(d5, sqrt(unname(cbind(parmarisk(sol5qpopt), parmarisk(sol5nlpopt)))))
	d5 = rbind(d5, cbind(parmareward(sol5qpopt), parmareward(sol5nlpopt)))
	d5 = rbind(d5, cbind(as.numeric(tictoc(sol5qpopt)), as.numeric(tictoc(sol5nlpopt))))
	rownames(d5)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d5) = c("QP", "NLP")
	
	zz <- file("parma_test2-5.txt", open="wt")
	sink(zz)	
	show(pspec5)
	show(pspec5q)
	print(round(d5,6))
	sink(type="message")
	sink()
	close(zz)
	
	pspec6 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[6], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			options = list(threshold = 999, moment=1), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol6lpopt = parmasolve(pspec6, type = "LP", solver="GLPK")
	sol6nlpopt = parmasolve(pspec6, type = "NLP")
	
	d6 = cbind(weights(sol6lpopt), weights(sol6nlpopt))
	d6 = rbind(d6, unname(cbind(parmarisk(sol6lpopt), parmarisk(sol6nlpopt))))
	d6 = rbind(d6, cbind(parmareward(sol6lpopt), parmareward(sol6nlpopt)))
	d6 = rbind(d6, cbind(as.numeric(tictoc(sol6lpopt)), as.numeric(tictoc(sol6nlpopt))))
	rownames(d6)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d6) = c("LP", "NLP")
	
	zz <- file("parma_test3-6.txt", open="wt")
	sink(zz)	
	show(pspec6)
	print(round(d6,6))
	sink(type="message")
	sink()
	close(zz)
	
	
	pspec7 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[6], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			options = list(threshold = 999, moment=2), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	
	sol7nlpopt = parmasolve(pspec7, type = "NLP")
	
	pspec7x = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[6], 
			target = mean(colMeans(Data)), targetType = "equality",
			riskType = c("minrisk", "optimal")[2], 
			options = list(threshold = 999, moment=3), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	sol7xnlpopt = parmasolve(pspec7x, type = "NLP")
	
	skewness = function(X) {apply(X, 2, FUN = function(x) mean((x-mean(x))^3)/sd(x)^3 )}
	mcorr = function(X){ apply(cor(X), 2, FUN = function(x) mean(x[x!=1])) }
	d7 = cbind(weights(sol6nlpopt), weights(sol7nlpopt), weights(sol7xnlpopt), apply(Data, 2, "sd"), skewness(Data), mcorr(Data) )
	d7 = rbind(d7, cbind(parmarisk(sol6nlpopt), parmarisk(sol7nlpopt), parmarisk(sol7xnlpopt), NA, NA, NA))
	d7 = rbind(d7, cbind(parmareward(sol6nlpopt), parmareward(sol7nlpopt), parmareward(sol7xnlpopt), NA, NA, NA))
	d7 = rbind(d7, cbind(as.numeric(tictoc(sol6nlpopt)), as.numeric(tictoc(sol7nlpopt)), as.numeric(tictoc(sol7xnlpopt)), NA, NA, NA))
	rownames(d7)[16:18] = c("risk", "reward", "elapsed(secs)")
	colnames(d7) = c("LPM[1]", "LPM[2]", "LPM[3]", "sd", "skewness", "av.corr")
	
	zz <- file("parma_test3-7.txt", open="wt")
	sink(zz)	
	show(pspec7)
	print(round(d7,6))
	sink(type="message")
	sink()
	close(zz)
	
	# Mean Squared Error LP-NLP, QP-NLP
	errmat = matrix(NA, ncol = 5, nrow = 4)
	colnames(errmat) = c("MAD", "MiniMax", "CVaR", "EV", "LPM[1]")
	rownames(errmat) = c("MSE[weights]", "MAE[weights]", "MaxE[weights]", "Err[risk]")
	errmat[1,1] = mean( (weights(sol1lpopt)-weights(sol1nlpopt))^2)
	errmat[2,1] = mean( abs(weights(sol1lpopt)-weights(sol1nlpopt)))
	errmat[3,1] = max( abs(weights(sol1lpopt)-weights(sol1nlpopt)))
	errmat[4,1] = parmarisk(sol1nlpopt)-parmarisk(sol1lpopt)
	
	errmat[1,2] = mean( (weights(sol2lpopt)-weights(sol2nlpopt))^2)
	errmat[2,2] = mean( abs(weights(sol2lpopt)-weights(sol2nlpopt)))
	errmat[3,2] = max( abs(weights(sol2lpopt)-weights(sol2nlpopt)))
	errmat[4,2] = parmarisk(sol2lpopt)-parmarisk(sol2nlpopt)
	
	errmat[1,3] = mean( (weights(sol3lpopt)-weights(sol3nlpopt))^2)
	errmat[2,3] = mean( abs(weights(sol3lpopt)-weights(sol3nlpopt)))
	errmat[3,3] = max( abs(weights(sol3lpopt)-weights(sol3nlpopt)))
	errmat[4,3] = parmarisk(sol3lpopt)-parmarisk(sol3nlpopt)
	
	errmat[1,4] = mean( (weights(sol5qpopt)-weights(sol5nlpopt))^2)
	errmat[2,4] = mean( abs(weights(sol5qpopt)-weights(sol5nlpopt)))
	errmat[3,4] = max( abs(weights(sol5qpopt)-weights(sol5nlpopt)))
	errmat[4,4] = parmarisk(sol5qpopt)-parmarisk(sol5nlpopt)
	
	errmat[1,5] = mean( (weights(sol6lpopt)-weights(sol6nlpopt))^2)
	errmat[2,5] = mean( abs(weights(sol6lpopt)-weights(sol6nlpopt)))
	errmat[3,5] = max( abs(weights(sol6lpopt)-weights(sol6nlpopt)))
	errmat[4,5] = parmarisk(sol6lpopt)-parmarisk(sol6nlpopt)
	
	zz <- file("parma_test3-8.txt", open="wt")
	sink(zz)	
	print(errmat)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time() - tic
	return(toc)
}


##### LPM for different moments with L/S and Leverage constraint #####
parma.test4 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	

	spec1 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk ="LPM", riskType = "optimal", 
			options = list(moment = 1, threshold = 999), 
			LB = rep(-0.5, 15), UB = rep(0.5, 15), leverage = 1)
	
	spec2 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk ="LPM", riskType = "optimal", 
			options = list(moment = 2, threshold = 999), 
			LB = rep(-0.5, 15), UB = rep(0.5, 15), leverage = 1)
	
	spec3 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk ="LPM", riskType = "optimal", 
			options = list(moment = 3, threshold = 999), 
			LB = rep(-0.5, 15), UB = rep(0.5, 15), leverage = 1)
	
	spec4 = parmaspec(scenario = Data, forecast = colMeans(Data), 
			risk ="LPM", riskType = "optimal", 
			options = list(moment = 4, threshold = 999), 
			LB = rep(-0.5, 15), UB = rep(0.5, 15), leverage = 1)
	# we can definitely control the degree of accuracy of the final
	# solution by adjusting the termination criteria.
	sol1nlp = parmasolve(spec1, solver.control=list(ftol_rel = 1e-12, 
					ftol_abs = 1e-12, xtol_rel = 1e-12, maxeval = 5000, 
					print_level=1))
	
	sol2nlp = parmasolve(spec2, solver.control=list(ftol_rel = 1e-12, 
					ftol_abs = 1e-12, xtol_rel = 1e-12, maxeval = 5000, 
					print_level=1))
	
	sol3nlp = parmasolve(spec3, solver.control=list(ftol_rel = 1e-12, 
					ftol_abs = 1e-12, xtol_rel = 1e-12, maxeval = 5000, 
					print_level=1))
	
	sol4nlp = parmasolve(spec4, solver.control=list(ftol_rel = 1e-12, 
					ftol_abs = 1e-12, xtol_rel = 1e-12, maxeval = 5000, 
					print_level=1))
	
	
	S1sol = round(cbind(weights(sol1nlp), weights(sol2nlp), weights(sol3nlp),  weights(sol4nlp)), 4)
	idx = unique(which(abs(S1sol)>0.001, arr.ind=TRUE)[,1])
	S1sol = S1sol[idx, ]
	S1sol = rbind(S1sol, 
			cbind(parmarisk(sol1nlp)/parmareward(sol1nlp), parmarisk(sol2nlp)/parmareward(sol2nlp), 
					parmarisk(sol3nlp)/parmareward(sol3nlp), parmarisk(sol4nlp)/parmareward(sol4nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "Risk/Reward"
	S1sol = rbind(S1sol, cbind(tictoc(sol1nlp), tictoc(sol2nlp), tictoc(sol3nlp), tictoc(sol4nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "Elapsed(sec)"
	colnames(S1sol) = c("LPM(1)", "LPM(2)", "LPM(3)", "LPM(4)")
	
	options(width = 120)
	zz <- file("parma_test4_LPM.txt", open="wt")
	sink(zz)	
	cat("\nLPM (Long-Short)")
	cat("\n--------------------------------------------\n")
	print(S1sol, digits=4)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
}

##### Benchmark Relative Optimization (MinRisk)#####
parma.test5 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	
	set.seed(12)
	wb = rexp(10, 2)
	wb = wb/sum(wb)
	idx = c(1,2,3,4,5,7,8,9,10,11)
	# we want to achieve an annualized excess of 2% over the benchmark
	# therefore :
	benchmark = Data[,idx] %*% wb
	forecast = colMeans(Data[,idx])
	activeforc = forecast - mean(benchmark)
	# Excess of benchmark target:
	target = 0.02/252
	# Require positive total weights:
	# LB = -wb (i.e. max actual deviation from benchmark leads to a zero weight)
	pspec1min = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[1], 
			riskType = c("minrisk", "optimal")[1], 
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol1lpmin = parmasolve(pspec1min, type = "LP", verbose = TRUE)
	sol1nlpmin = parmasolve(pspec1min, type = "NLP", solver.control=list(print_level=1))
	
	
	pspec2min = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[2], 
			riskType = c("minrisk", "optimal")[1], 
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol2lpmin = parmasolve(pspec2min, type = "LP", verbose = TRUE)
	sol2nlpmin = parmasolve(pspec2min, type = "NLP", solver.control=list(print_level=1))
	
	
	pspec3min = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[1],
			options=list(alpha=0.05),
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol3lpmin = parmasolve(pspec3min, type = "LP", verbose = TRUE)
	sol3nlpmin = parmasolve(pspec3min, type = "NLP", solver.control=list(print_level=1))
	
	
	pspec4min = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[4], 
			riskType = c("minrisk", "optimal")[1],
			options=list(alpha=0.05),
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol4lpmin = parmasolve(pspec4min, type = "LP", verbose = TRUE)
	
	
	idx = c(1,2,3,4,5,7,8,9,10,11)
	# we want to achieve an annualized excess of 2% over the benchmark
	# therefore :
	benchmark = Data[,6]
	forecast = colMeans(Data[,idx])
	activeforc = forecast - mean(benchmark)
	# Excess of benchmark target:
	target = 0.02/252
	
	pspec5min = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[1],
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	S = cov(Data[,idx])
	bS = cov(cbind(benchmark, Data[,idx]))[1,1:11]
	pspec5qmin = parmaspec(S = S, benchmarkS = bS, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[1],
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol5nlpmin = parmasolve(pspec5min, type = "NLP", solver.control=list(print_level=1))
	sol5qpmin = parmasolve(pspec5qmin, type = "QP")
	
	# LPM does not admit a benchmark, instead pass the excesss benchmark returns
	# and set threshold to zero
	pspec6min = parmaspec(scenario = Data[,idx] - matrix(benchmark, ncol = 10, nrow = nrow(Data)), 
			benchmark = NULL, 
			forecast = activeforc, target = target,
			targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[6], 
			riskType = c("minrisk", "optimal")[1],
			options=list(moment=1, threshold = 0), LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol6lpmin = parmasolve(pspec6min, type = "LP", verbose = TRUE)
	sol6nlpmin = parmasolve(pspec6min, type = "NLP", solver.control=list(print_level=1))
	
	M = matrix(NA, ncol = 6, nrow = 17)
	colnames(M) =  c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM[1]")
	rownames(M) = c(colnames(Data[,idx]), "Active_Risk", "Active_Reward", "w-MSE(vsNLP)", "w-MAE(vsNLP)", "w-maxErr(vsNLP)", 
			"r-Diff(vsNLP)", "t-Diff(vsNLP)")
	M[1:10, ] = cbind(weights(sol1lpmin), weights(sol2lpmin), weights(sol3lpmin), 
			weights(sol4lpmin), weights(sol5qpmin), weights(sol6lpmin))
	M[11, ] = c(parmarisk(sol1lpmin), parmarisk(sol2lpmin), parmarisk(sol3lpmin), 
			parmarisk(sol4lpmin), sqrt(parmarisk(sol5qpmin)), parmarisk(sol6lpmin))
	M[12, ] = c(parmareward(sol1lpmin), parmareward(sol2lpmin), parmareward(sol3lpmin), 
			parmareward(sol4lpmin), parmareward(sol5qpmin), parmareward(sol6lpmin))
	# LV(or QP) vs NLP
	nlpw = cbind(weights(sol1nlpmin), weights(sol2nlpmin), weights(sol3nlpmin), 
			rep(NA, 10), weights(sol5nlpmin), weights(sol6nlpmin))
	nlpr = c(parmarisk(sol1nlpmin), parmarisk(sol2nlpmin), parmarisk(sol3nlpmin), 
			NA, sqrt(parmarisk(sol5nlpmin)), parmarisk(sol6nlpmin))
	
	M[13,] = apply((M[1:10,] - nlpw)^2, 2, FUN = function(x) mean(x))
	M[14,] = apply(abs(M[1:10,] - nlpw), 2, FUN = function(x) mean(x))
	M[15,] = apply(abs(M[1:10,] - nlpw), 2, FUN = function(x) max(x))
	M[16,] = M[11,] - nlpr
	M[17,] = as.numeric(c(tictoc(sol1lpmin)-tictoc(sol1nlpmin),
			tictoc(sol2lpmin)-tictoc(sol2nlpmin),
			tictoc(sol3lpmin)-tictoc(sol3nlpmin),
			NA,
			tictoc(sol5qpmin)-tictoc(sol5nlpmin),
			tictoc(sol6lpmin)-tictoc(sol6nlpmin)))
	
	options(width = 120)
	zz <- file("parma_test5_Benchmark(MinRisk).txt", open="wt")
	sink(zz)
	cat("\nActive Weights and Measures\n")
	print(round(M,5))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
}

##### Benchmark Relative Optimization (OptRisk)#####
parma.test6 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	
	set.seed(12)
	wb = rexp(10, 2)
	wb = wb/sum(wb)
	idx = c(1,2,3,4,5,7,8,9,10,11)
	# we want to achieve an annualized excess of 2% over the benchmark
	# therefore :
	benchmark = Data[,idx] %*% wb
	forecast = colMeans(Data[,idx])
	activeforc = forecast - mean(benchmark)
	# Excess of benchmark target:
	target = 0.02/252
	# Require positive total weights:
	# LB = -wb (i.e. max actual deviation from benchmark leads to a zero weight)
	pspec1opt = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[1], 
			riskType = c("minrisk", "optimal")[2], 
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol1lpopt = parmasolve(pspec1opt, type = "LP", verbose = TRUE)
	sol1nlpopt = parmasolve(pspec1opt, type = "NLP", solver.control=list(print_level=1))
	
	
	pspec2opt = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc,
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[2], 
			riskType = c("minrisk", "optimal")[2], 
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol2lpopt = parmasolve(pspec2opt, type = "LP", verbose = TRUE)
	sol2nlpopt = parmasolve(pspec2opt, type = "NLP", solver.control=list(print_level=1,xtol_rel=1e-14, ftol_rel=1e-14))
	
	
	pspec3opt = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc,
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[2],
			options=list(alpha=0.05),
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol3lpopt = parmasolve(pspec3opt, type = "LP", verbose = TRUE)
	sol3nlpopt = parmasolve(pspec3opt, type = "NLP", solver.control=list(print_level=1))
	
	
	pspec4opt = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[4], 
			riskType = c("minrisk", "optimal")[2],
			options=list(alpha=0.05),
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol4lpopt = parmasolve(pspec4opt, type = "LP", verbose = TRUE)
	
	
	idx = c(1,2,3,4,5,7,8,9,10,11)
	# we want to achieve an annualized excess of 2% over the benchmark
	# therefore :
	benchmark = Data[,6]
	forecast = colMeans(Data[,idx])
	activeforc = forecast - mean(benchmark)
	# Excess of benchmark target:
	target = 0.02/252
	
	pspec5opt = parmaspec(scenario = Data[,idx], benchmark = benchmark, 
			forecast = activeforc, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[2],
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	S = cov(Data[,idx])
	bS = cov(cbind(benchmark, Data[,idx]))[1,1:11]
	pspec5qopt = parmaspec(S = S, benchmarkS = bS, 
			forecast = activeforc, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[2],
			LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol5nlpopt = parmasolve(pspec5opt, type = "NLP", solver.control=list(print_level=1))
	sol5qpopt = parmasolve(pspec5qopt, type = "QP")
	
	# LPM does not admit a benchmark, instead pass the excesss benchmark returns
	# and set threshold to zero
	pspec6opt = parmaspec(scenario = Data[,idx] - matrix(benchmark, ncol = 10, nrow = nrow(Data)), 
			benchmark = NULL, forecast = activeforc, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[6], 
			riskType = c("minrisk", "optimal")[2],
			options=list(moment=1, threshold = 0), LB = -wb, UB = rep(0.2, 10), budget = 0)
	
	sol6lpopt = parmasolve(pspec6opt, type = "LP", verbose = TRUE)
	sol6nlpopt = parmasolve(pspec6opt, type = "NLP", solver.control=list(print_level=1))
	
	M = matrix(NA, ncol = 6, nrow = 17)
	colnames(M) =  c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM[1]")
	rownames(M) = c(colnames(Data[,idx]), "Active_Risk", "Active_Reward", "w-MSE(vsNLP)", "w-MAE(vsNLP)", "w-maxErr(vsNLP)", 
			"r-Diff(vsNLP)", "t-Diff(vsNLP)")
	M[1:10, ] = cbind(weights(sol1lpopt), weights(sol2lpopt), weights(sol3lpopt), 
			weights(sol4lpopt), weights(sol5qpopt), weights(sol6lpopt))
	M[11, ] = c(parmarisk(sol1lpopt), parmarisk(sol2lpopt), parmarisk(sol3lpopt), 
			parmarisk(sol4lpopt), sqrt(parmarisk(sol5qpopt)), parmarisk(sol6lpopt))
	M[12, ] = c(parmareward(sol1lpopt), parmareward(sol2lpopt), parmareward(sol3lpopt), 
			parmareward(sol4lpopt), parmareward(sol5qpopt), parmareward(sol6lpopt))
	# LV(or QP) vs NLP
	nlpw = cbind(weights(sol1nlpopt), weights(sol2nlpopt), weights(sol3nlpopt), 
			rep(NA, 10), weights(sol5nlpopt), weights(sol6nlpopt))
	nlpr = c(parmarisk(sol1nlpopt), parmarisk(sol2nlpopt), parmarisk(sol3nlpopt), 
			NA, sqrt(parmarisk(sol5nlpopt)), parmarisk(sol6nlpopt))
	
	M[13,] = apply((M[1:10,] - nlpw)^2, 2, FUN = function(x) mean(x))
	M[14,] = apply(abs(M[1:10,] - nlpw), 2, FUN = function(x) mean(x))
	M[15,] = apply(abs(M[1:10,] - nlpw), 2, FUN = function(x) max(x))
	M[16,] = M[11,] - nlpr
	M[17,] = as.numeric(c(tictoc(sol1lpopt)-tictoc(sol1nlpopt),
					tictoc(sol2lpopt)-tictoc(sol2nlpopt),
					tictoc(sol3lpopt)-tictoc(sol3nlpopt),
					NA,
					tictoc(sol5qpopt)-tictoc(sol5nlpopt),
					tictoc(sol6lpopt)-tictoc(sol6nlpopt)))
	
	options(width = 120)
	zz <- file("parma_test6_Benchmark(OptRisk).txt", open="wt")
	sink(zz)
	cat("\nActive Weights and Measures\n")
	print(round(M,5))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
}

##### Custom Constraints (LP+NLP) #####
parma.test6 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	m = 15
	
	# Constraint 1:
	# The joint exposure to (1, 2, 8, 9, 10) to be less than 40%
	# The joint exposure to (3, 4, 7, 11) to be less than 40%
	# The joint exposure to (5, 6, 12, 13, 14) to be less than 40%
	s1 = s2 = s3 = rep(0, m-1)
	s1[c(1, 2, 8, 9, 10)] = 1
	s2[c(3, 4, 7, 11)] = 1
	s3[c(5, 6, 12, 13, 14)] = 1
	
	sectmatrix = cbind(rbind(t(s1), t(s2), t(s3)))
	
	# Unlike the NLP formulation, we do not have to worry about the
	# extra CVaR parameter nor the fractional scaling parameter which are
	# handled internally by the LP problem setup
	
	pspec1 = parmaspec(scenario = Data[,-6], forecast = colMeans(Data[,-6]), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 14), UB = rep(0.2, 14), budget = 1)
	
	pspec2 = parmaspec(scenario = Data[,-6], forecast = colMeans(Data[,-6]), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 14), UB = rep(0.2, 14), budget = 1, 
			ineq.mat = sectmatrix, ineq.LB = c(0,0,0), ineq.UB = c(0.4, 0.4, 0.4))
	
	sol1lp = parmasolve(pspec1, type = "LP", solver="GLPK")
	# spec 2 is LP because of the ineq.mat LP matrices
	sol2lp = parmasolve(pspec2, type = "LP", solver="GLPK")
	
	
	###############################################
	# NLP formulation:
	# User list:
	# Any constraint must take as arguments w, optvars, uservars
	# w - weights (see example)
	# optvars - program specific values (user should not tamper with this)
	# e.g. contains the VaR quantile level, LPM moment and threshold
	# uservars - list supplied by user containing any extra values required
	# to evaluate the constraints
	
	# In the nloptr solver, the inequality is in the form h(x)<=0
	
	# The actual weights are indexed by optvars$widx
	# depending on the problem, there may be other
	# parameters which are optimized and not related
	# to the weights (e.g. VaR in CVaR problem, and
	# the multiplier in all fractional risk problems).
	# -->extra parameters are indexed by: optvars$vidx
	# -->multiplier indexed by:           optvars$midx
	# For the fractional problem, ineqLB and ineqUB are multiplied by the fractional
	# multiplier: m * ineqLB <= ineq <= m * ineqUB
	uservars = list()
	uservars$sectmatrix = sectmatrix
	secfun = function(w, optvars, uservars){
		ineqc = w[optvars$midx]*c(0, 0, 0) -uservars$sectmatrix %*% w[optvars$widx]
		ineqc = c(ineqc, uservars$sectmatrix %*% w[optvars$widx] - w[optvars$midx]*c(0.4, 0.4, 0.4))
		return(ineqc)
	}
	
	# gradient
	secjac = function(w, optvars, uservars){
		# size of problem: optvars$fm
		# nrows = n.constraints = 3x2
		widx = optvars$widx
		midx = optvars$midx
		g = matrix(0, ncol = optvars$fm, nrow = 6)
		g[1:3, widx] = -uservars$sectmatrix
		g[4:6, widx] =  uservars$sectmatrix
		g[1:3, midx] = c(0,0,0)
		g[4:6, midx] = c(-0.4,-0.4,-0.4)
		return(g)
	}
	
	pspec3 = parmaspec(scenario = Data[,-6], forecast = colMeans(Data[,-6]), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 14), UB = rep(0.2, 14), budget = 1, 
			ineqfun = list(secfun), ineqgrad = list(secjac), uservars = uservars)
	
	sol3nlp = parmasolve(pspec3, type="NLP", solver.control=list(print_level=1))
	
	
	S1sol = round(cbind(weights(sol1lp), weights(sol2lp), weights(sol3nlp)), 4)
	idx = unique(which(S1sol>0.001, arr.ind=TRUE)[,1])
	S1sol = S1sol[idx, ]
	S1sol = rbind(S1sol, cbind(parmarisk(sol1lp), parmarisk(sol2lp), parmarisk(sol3nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "CvaR"
	S1sol = rbind(S1sol, cbind(sol1lp@solution$VaR, sol2lp@solution$VaR, sol3nlp@solution$VaR))
	rownames(S1sol)[dim(S1sol)[1]] = "VaR"
	S1sol = rbind(S1sol, cbind(parmareward(sol1lp), parmareward(sol2lp), parmareward(sol3nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "Reward"
	S1sol = rbind(S1sol, cbind(parmarisk(sol1lp)/parmareward(sol1lp), 
					parmarisk(sol2lp)/parmareward(sol2lp), 
					parmarisk(sol3nlp)/parmareward(sol3nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "Risk/Reward"
	xc1 = as.numeric(sectmatrix %*% weights(sol1lp))
	xc2 = as.numeric(sectmatrix %*% weights(sol2lp))
	xc3 = as.numeric(sectmatrix %*% weights(sol3nlp))
	S1sol = rbind(S1sol, cbind(xc1[1], xc2[1], xc3[1]))
	rownames(S1sol)[dim(S1sol)[1]] = "Sect_1"
	S1sol = rbind(S1sol, cbind(xc1[2], xc2[2], xc3[2]))
	rownames(S1sol)[dim(S1sol)[1]] = "Sect_2"
	S1sol = rbind(S1sol, cbind(xc1[3], xc2[3], xc3[3]))
	rownames(S1sol)[dim(S1sol)[1]] = "Sect_3"
	S1sol = rbind(S1sol, cbind(tictoc(sol1lp), tictoc(sol2lp), tictoc(sol3nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "Elapsed(sec)"
	colnames(S1sol) = c("LP (nc)", "LP (c)", "NLP (c)")
	
	zz <- file("parma_test6_userconstraints.txt", open="wt")
	sink(zz)
	cat("\nCVaR with user constraints")
	cat("\n--------------------------------------------\n")
	print(round(S1sol,4))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
	
}

##### Custom Constraints (NLP) #####
parma.test7 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	
	m = 15
	
	# User Constraint Examples
	# 1. Sector Weights
	# 2. Quadratic Constraint based on some PD user supplied matrix (purely
	# illustrative).
	# See example 6 above for NLP custom constraints setup
	
	uservars =  list()
	#--------------------------------------------------------------------------
	# Constraint 1:
	s1 = s2 = s3 = rep(0, m-1)
	s1[c(1, 2, 8, 9, 10)] = 1
	s2[c(3, 4, 7, 11)] = 1
	s3[c(5, 6, 12, 13, 14)] = 1
	sectmatrix = cbind(rbind(t(s1), t(s2), t(s3)))
	
	uservars$sectmatrix = sectmatrix
	
	secfun = function(w, optvars, uservars){
		ineqc = w[optvars$midx]*c(0, 0, 0) -uservars$sectmatrix %*% w[optvars$widx]
		ineqc = c(ineqc, uservars$sectmatrix %*% w[optvars$widx] - w[optvars$midx]*c(0.4, 0.4, 0.4))
		return(ineqc)
	}
	secjac = function(w, optvars, uservars){
		# size of problem: optvars$fm
		# nrows = n.constraints = 3x2
		widx = optvars$widx
		midx = optvars$midx
		g = matrix(0, ncol = optvars$fm, nrow = 6)
		g[1:3, widx] = -uservars$sectmatrix
		g[4:6, widx] =  uservars$sectmatrix
		g[1:3, midx] = c(0,0,0)
		g[4:6, midx] = c(-0.4,-0.4,-0.4)
		return(g)
	}
	#--------------------------------------------------------------------------
	# Constraint 2:
	# The joint variance of (1,2,8,9,10) <= 0.01^2
	S = cov(Data[,-6])[c(1, 2, 8, 9, 10), c(1, 2, 8, 9, 10)]
	uservars$Smatrix = S
	
	sqfun = function(w, optvars, uservars){
		x = w[optvars$widx]
		return( as.numeric( x[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix %*% x[c(1, 2, 8, 9, 10)]) - w[optvars$midx]*(0.01)^2 )
	}
	# --> Gradient (MUST BE PROVIDED) 
	sqjac = function(w, optvars, uservars){
		g = matrix(0, ncol = optvars$fm, nrow = 1)
		# remember that widx = index of weights in the w vector (NOT necesarily 1:m!)
		widx = optvars$widx
		x = w[widx]
		g[1, widx[c(1, 2, 8, 9, 10)]] = as.numeric(2*x[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix)
		g[1,optvars$midx] = 0.01^2
		return(g)
	}
	
	pspec = parmaspec(scenario = Data[,-6], forecast = colMeans(Data), 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 14), UB = rep(0.5, 14), budget = 1,
			ineqfun = list(secfun, sqfun), ineqgrad = list(secjac, sqjac), 
			uservars = uservars)
	
	sol1nlp = parmasolve(pspec, type="NLP", solver.control=list(print_level=1))
	
	### Same example with minrisk type (no multiplier-->no midx)
	secfun = function(w, optvars, uservars){
		ineqc = c(0, 0, 0) -uservars$sectmatrix %*% w[optvars$widx]
		ineqc = c(ineqc, uservars$sectmatrix %*% w[optvars$widx] - c(0.4, 0.4, 0.4))
		return(ineqc)
	}
	secjac = function(w, optvars, uservars){
		# size of problem: optvars$fm
		# nrows = n.constraints = 3x2
		# optvars$fm = size of w
		widx = optvars$widx
		g = matrix(0, ncol = optvars$fm, nrow = 6)
		g[1:3, widx] = -uservars$sectmatrix
		g[4:6, widx] =  uservars$sectmatrix
		return(g)
	}
	#--------------------------------------------------------------------------
	# Constraint 2:
	# The joint variance of (1,2,8,9,10) <= 0.01^2
	S = cov(Data[,-6])[c(1, 2, 8, 9, 10), c(1, 2, 8, 9, 10)]
	uservars$Smatrix = S
	# make the inequality less than or equal to solution of previous optimal
	# weights(sol1nlp)[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix %*% weights(sol1nlp)[c(1, 2, 8, 9, 10)]
	sqfun = function(w, optvars, uservars){
		x = w[optvars$widx]
		return( as.numeric( x[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix %*% x[c(1, 2, 8, 9, 10)]) - 6.465135e-06 )
	}
	# --> Gradient (MUST BE PROVIDED) 
	sqjac = function(w, optvars, uservars){
		g = matrix(0, ncol = optvars$fm, nrow = 1)
		# remember that widx = index of weights in the w vector (NOT necesarily 1:m!)
		widx = optvars$widx
		x = w[widx]
		g[1, widx[c(1, 2, 8, 9, 10)]] = as.numeric(2*x[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix)
		return(g)
	}
	# set the target equal to the reward of the optimal solution
	pspec2 = parmaspec(scenario = Data[,-6], forecast = colMeans(Data), 
			target = parmareward(sol1nlp), targetType = "equality",
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[1], 
			options = list(alpha = 0.05), 
			LB = rep(0, 14), UB = rep(0.5, 14), budget = 1,
			ineqfun = list(secfun, sqfun), ineqgrad = list(secjac, sqjac), 
			uservars = uservars)
	
	sol2nlp = parmasolve(pspec2, type="NLP", solver.control=list(print_level=1))

	
	# Examine the 2 solutions (they will be close/equal).
	# Note: we did not use an equality for the quadratic variance constraint
	# since for a convex problem we require strictly affine equality constraints.
	
	w1 = weights(sol1nlp)
	sol1C1 = as.numeric(sectmatrix %*% w1)
	sol1C2 = w1[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix %*% w1[c(1, 2, 8, 9, 10)]
	
	w2 = weights(sol2nlp)
	sol2C1 = as.numeric(sectmatrix %*% w2)
	sol2C2 = w2[c(1, 2, 8, 9, 10)] %*% uservars$Smatrix %*% w2[c(1, 2, 8, 9, 10)]
	
	r1 = parmarisk(sol1nlp)
	r2 = parmarisk(sol2nlp)
	
	rw1 = parmareward(sol1nlp)
	rw2 = parmareward(sol2nlp)
	
	
	S1sol = round(cbind(w1, w2), 4)
	idx = unique(which(S1sol>0.001, arr.ind=TRUE)[,1])
	S1sol = S1sol[idx, ]
	S1sol = rbind(S1sol, cbind(r1, r2))
	rownames(S1sol)[dim(S1sol)[1]] = "CvaR"
	S1sol = rbind(S1sol, cbind(sol1nlp@solution$VaR, sol2nlp@solution$VaR))
	rownames(S1sol)[dim(S1sol)[1]] = "VaR"
	S1sol = rbind(S1sol, cbind(rw1, rw2))
	rownames(S1sol)[dim(S1sol)[1]] = "Reward"
	S1sol = rbind(S1sol, cbind(r1/rw1, r2/rw2))
	rownames(S1sol)[dim(S1sol)[1]] = "Risk/Reward"
	S1sol = rbind(S1sol, cbind(sol1C1[1], sol2C1[1]))
	rownames(S1sol)[dim(S1sol)[1]] = "Sect_1"
	S1sol = rbind(S1sol, cbind(sol1C1[2], sol2C1[2]))
	rownames(S1sol)[dim(S1sol)[1]] = "Sect_2"
	S1sol = rbind(S1sol, cbind(sol1C1[3], sol2C1[3]))
	rownames(S1sol)[dim(S1sol)[1]] = "Sect_3"
	S1sol = rbind(S1sol, cbind(sol1C2[1]*1e5, sol2C2[1]*1e5))
	rownames(S1sol)[dim(S1sol)[1]] = "varQ_1(x1e5)"
	S1sol = rbind(S1sol, cbind(tictoc(sol1nlp), tictoc(sol2nlp)))
	rownames(S1sol)[dim(S1sol)[1]] = "Elapsed(sec)"
	colnames(S1sol) = c("NLP (opt-c)", "NLP (min-c)")

	# Risk/Reward is NOT "higher" in the min formulation since the var 
	# constraint is NOT exactly the same. The small differences in the 
	# Risk/Reward is attributable to the small differences in the lower varQ 
	# constraint in the fractional formulation.
	zz <- file("parma_test7_userconstraints.txt", open="wt")
	sink(zz)
	cat("\nCVaR with NLP user constraints")
	cat("\n--------------------------------------------\n")
	print(round(S1sol,4))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
}

##### The efficient frontier and fractional programming example #####
parma.test9 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)
	
	f = colMeans(Data)
	pspec1 = parmaspec(scenario = Data, target = 0, 
			targetType = "equality", forecast = f, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[1], 
			options = list(alpha = 0.05), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	flp = parmafrontier(pspec1, n.points = 100, type = "LP", solver="GLPK")
	r1 = flp[,"CVaR"]
	rw1 = flp[,"reward"]
	opt = which(r1/rw1 == min(r1/rw1))
	
	fnlp = parmafrontier(pspec1, n.points = 100, type = "NLP")
	r2 = fnlp[,"CVaR"]
	rw2 = fnlp[,"reward"]
	
	pspec2 = parmaspec(scenario = Data,  forecast = f, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[3], 
			riskType = c("minrisk", "optimal")[2], 
			options = list(alpha = 0.05), 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	optlp = parmasolve(pspec2, type="LP", solver="GLPK")
	
	postscript("parma_test9-1.eps", width = 12, height = 8)
	plot(r1, rw1, type = "p", col = "steelblue", xlab = "E[Risk]", ylab = "E[Return]", main = "CVaR Frontier")
	lines(r2, rw2, col = "tomato1")
	points(r1[opt], rw1[opt], col = "orange", pch = 4, lwd = 3)
	points(parmarisk(optlp), parmareward(optlp), col = "yellow", pch = 6, lwd = 3)
	legend("topleft", c("LP Frontier", "NLP Frontier", paste("Frontier Optimal (",round(r1[opt]/rw1[opt], 3),")",sep=""), 
					paste("Fractional Optimal (",round(parmarisk(optlp)/parmareward(optlp), 3),")",sep="")), 
			col = c("steelblue", "tomato1", "orange", "yellow"), pch = c(1, -1, 4, 6), 
			lwd = c(1,1,3,3), lty=c(0,1,0,0), bty = "n", cex = 0.9)
	grid()
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	dev.off()
	
	
	###########################################################################
	# One more with QP based solver (which shows that the optimal solution for
	# standard deviation NOT variance because it is superadditive rather than
	# subadditive)	
	S = cov(Data)
	pspec2 = parmaspec(S = S, target = 0, 
			targetType = "equality", forecast = f, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[1], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	fqp = parmafrontier(pspec2, n.points = 200)
	r1 =  sqrt(na.omit(fqp[,"EV"]))
	rw1 = na.omit(fqp[,"reward"])
	opt = which(r1/rw1 == min(r1/rw1))
	
	pspec2 = parmaspec(scenario = Data, target = 0, 
			targetType = "equality", forecast = f, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[1], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	fnlp = parmafrontier(pspec2, n.points = 200, type = "NLP")
	
	r2 =  sqrt(na.omit(fnlp[,"EV"]))
	rw2 = na.omit(fnlp[,"reward"])
	opt2 = which(r2/rw2 == min(r2/rw2))
	
	pspec2 = parmaspec(S = S,  forecast = f, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[2], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	optqp = parmasolve(pspec2)
	
	pspec2 = parmaspec(scenario = Data,  forecast = f, 
			risk = c("MAD", "MiniMax", "CVaR", "CDaR", "EV", "LPM")[5], 
			riskType = c("minrisk", "optimal")[2], 
			LB = rep(0, 15), UB = rep(1, 15), budget = 1)
	optnlp = parmasolve(pspec2)
	
	postscript("parma_test9-2.eps", width = 12, height = 8)
	plot(r1, rw1, type = "p", col = "steelblue", xlab = "E[Risk]", ylab = "E[Return]", main = "EV Frontier")
	lines(r2, rw2, col = "tomato1")
	points(r1[opt], rw1[opt], col = "orange", pch = 4, lwd = 3)
	points(sqrt(parmarisk(optqp)), parmareward(optqp), col = "yellow", pch = 6, lwd = 3)
	legend("topleft", c("QP Frontier", "NLP Frontier", paste("Frontier Optimal (",round(r1[opt]/rw1[opt], 3),")",sep=""), 
					paste("Fractional Optimal (",round(sqrt(parmarisk(optqp))/parmareward(optqp), 3),")",sep="")), 
			col = c("steelblue", "tomato1", "orange", "yellow"), pch = c(1, -1, 4, 6), 
			lwd = c(1,1,3,3), lty=c(0,1,0,0), bty = "n", cex = 0.9)
	grid()
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	dev.off()


	toc = Sys.time() - tic
	return(toc)

}

##### CARA based optimization #####
parma.test10 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	Data = etfdata/lag(etfdata)-1
	Data = na.omit(Data)	
	require(rmgarch)
	require(parallel)
	m = c(1,2,3,4,5,7,8,9,10, 11, 12, 13)
	
	cl = makePSOCKcluster(5)
	spec = gogarchspec(
			mean.model = list(model = c("constant", "AR", "VAR")[2], lag = 1), 
			distribution.model = c("mvnorm", "manig", "magh")[2], 
			ica = "fastica")
	# Key Note: roll is like out.sample --> we keep 200 points for out of sample
	# forecasting (tail(Data, 200)), and the function returns 200+1 forecasts
	# (the 200 from the in-sample dataset and 1 given the last datapoint).
	S1 = fmoments(Data = Data[, m], n.ahead = 1, roll = 200, spec = spec, 
			solver = "hybrid", cluster = cl, gfun = "tanh", maxiter1 = 100000, 
			rseed = 500)
	
	w = matrix(NA, ncol = 12, nrow = 200)
	port = rep(0, 200)
	actual = tail(Data[,m], 200)
	for(i in 1:200){
		# 4 Moment Approximation
		tmp = parmautility(U = "CARA", method = "moment", 
				M1 = S1@moments$forecastMu[i,], M2 =  S1@moments$forecastCov[,,i], 
				M3 = S1@moments$forecastM3[,,i], M4 = S1@moments$forecastM4[,,i], 
				RA = 25, budget = 1, LB = rep(0, 12), UB = rep(0.5, 12))
		w[i,] = weights(tmp)
	}
	port = as.numeric(apply(cbind(w, actual), 1, FUN = function(x) sum(x[1:12]*x[13:24])))
	ewport = as.numeric(apply(actual, 1, "mean"))
	ew = rep(1/12, 12)
	pexmu = pexsigma = pexpkurt = pexpskew = rep(0, 200)
	ewsigma = ewkurt = ewskew = rep(0, 200)
	for(i in 1:200){
		pexsigma[i] = sqrt(w[i,] %*% S1@moments$forecastCov[,,i] %*% w[i,])
		ewsigma[i] = sqrt(ew %*% S1@moments$forecastCov[,,i] %*% ew)

		pexpskew[i] = (w[i,] %*% S1@moments$forecastM3[,,i] %*% kronecker(w[i,], w[i,]))/pexsigma[i]^2
		ewskew[i] = (ew %*% S1@moments$forecastM3[,,i] %*% kronecker(ew, ew))/ewsigma[i]^2
		
		pexpkurt[i] = (w[i,] %*% S1@moments$forecastM4[,,i] %*% kronecker(w[i,], kronecker(w[i,],w[i,])))/pexsigma[i]^4
		ewkurt[i] = (ew %*% S1@moments$forecastM4[,,i] %*% kronecker(ew, kronecker(ew,ew)))/ewsigma[i]^4
	}
	
	postscript("parma_test10.eps", width = 12, height = 8)
	par(mfrow = c(2,2))
	plot(cumprod(1+port), type = "l", col = "steelblue", xlab = "Time", ylab = "Wealth", 
			main = "CARA Utility Optimization", ylim = c(0.85, 1.3))
	lines(cumprod(1+ewport), col = "tomato1")
	legend("topleft", c("CARA (4 Mom)", "EqualWeight"), 
			col = c("steelblue", "tomato1"), lty=c(1,1), bty = "n", cex = 0.9)
	grid()
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	
	plot(pexsigma, type = "l", col = "steelblue", xlab = "Time", ylab = "Forecast Optimal Sigma", 
			main = "Optimal Sigma", ylim = c(min(min(pexsigma), min(ewsigma)), max(max(pexsigma), max(ewsigma))))
	lines(ewsigma, col = "tomato1")
	legend("topright", c("CARA (4 Mom)", "EqualWeight"), 
			col = c("steelblue", "tomato1"), lty=c(1,1), bty = "n", cex = 0.9)
	grid()
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	
	plot(pexpskew, type = "l", col = "steelblue", xlab = "Time", ylab = "Forecast Optimal Skewness", 
			main = "Optimal Skewness", ylim = c(min(min(pexpskew), min(ewskew)), max(max(pexpskew), max(ewskew))))
	lines(ewskew, col = "tomato1")
	legend("topleft", c("CARA (4 Mom)", "EqualWeight"), 
			col = c("steelblue", "tomato1"), lty=c(1,1), bty = "n", cex = 0.9)
	grid()
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	
	plot(pexpkurt, type = "l", col = "steelblue", xlab = "Time", ylab = "Forecast Optimal Kurtosis", 
			main = "Optimal Kurtosis", ylim = c(min(min(pexpkurt), min(ewkurt)), max(max(pexpkurt), max(ewkurt))))
	lines(ewkurt, col = "tomato1")
	legend("topleft", c("CARA (4 Mom)", "EqualWeight"), 
			col = c("steelblue", "tomato1"), lty=c(1,1), bty = "n", cex = 0.9)
	grid()
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	dev.off()
	
	toc = Sys.time() - tic
	return(toc)
}
# SOCP Tests
parma.test11 = function()
{
	library(xts)
	data("etfdata")
	R = etfdata/lag(etfdata)-1
	R = na.omit(R)
	
	sectors = matrix(rep(0, 15), nrow=1)
	sectors[1,c(1,4,5,10,12,15)]= 1
	sectorsLB = 0
	sectorsUB = 0.3
	S = cov(coredata(R))
	fmu = colMeans(coredata(R))
	fmub = median(fmu)
	
	control=list(abs.tol = 1e-12, rel.tol = 1e-12, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	
	# Example 1: Minimize Risk, Subject to return
	# plus inequality constraint
	LB = rep(0, 15)
	UB = rep(0.5, 15)
	spec = parmaspec(S = S, forecast = fmu, target = fmub, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "equality",
			ineq.mat = sectors, ineq.LB = 0, ineq.UB = 0.3, budget = 1)
	
	sol.socp = parmasolve(spec, solver.control = control, type="SOCP")
	sol.qp = parmasolve(spec, type = "QP")
	
	zz <- file("parma_test11a_minrisk_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP - QP")
	cat("\n--------------------------------------------\n")
	print(round(matrix(weights(sol.socp) - weights(sol.qp), ncol=1), 5))
	sink(type="message")
	sink()
	close(zz)
	
	# frontier
	spec = parmaspec(S = S, forecast = fmu, target = NULL, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "equality",
			budget = 1)
	control=list(abs.tol = 1e-8, rel.tol = 1e-8, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	frontsocp = parmafrontier(spec, n.points = 100, miny = NULL, maxy = NULL, 
			type = "SOCP", solver.control = control)
	
	frontqp = parmafrontier(spec, n.points = 100, miny = NULL, maxy = NULL, 
			type = "QP", solver.control = control)
	# ok to ignore warnings...just indicating points outside the feasible frontier
	postscript("parma_test11a.eps", width = 12, height = 8)
	plot(frontqp[,"EV"], frontqp[,"reward"], type="l", xlab = "risk", ylab="reward",
			main = "EV Frontier\n[QP and SOCP]", cex.main=0.9)
	points(frontsocp[,"EV"], frontsocp[,"reward"], col = 2)
	legend("topleft", c("QP","SOCP"), col = 1:2, bty="n", lty=c(1,0), pch=c(-1,1))
	dev.off()
	
	# Example 2: Minimize Risk, Subject to return
	# plus inequality constraint and benchmark
	control=list(abs.tol = 1e-9, rel.tol = 1e-9, Nu=5, max.iter=1250, BigM.K = 4, BigM.iter = 15)
	LB = rep(0, 15)
	UB = rep(0.5, 15)
	set.seed(100)
	wb = rexp(15)
	wb = wb/sum(wb)
	B = apply(R, 1, function(x) sum(x*wb))
	bS = c(var(B), cov(B, R))
	
	spec = parmaspec(S = S, benchmarkS = bS, 
			forecast = fmu - mean(B), target = fmub-mean(B), LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "inequality",
			ineq.mat = sectors, ineq.LB = 0, ineq.UB = 0.5, budget = 1)
	
	sol.socp = parmasolve(spec, solver.control = control, type="SOCP")
	sol.qp = parmasolve(spec, type = "QP")
	
	zz <- file("parma_test11b_minrisk_benchmark_relative_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP - QP")
	cat("\n--------------------------------------------\n")
	print(round(matrix(weights(sol.socp) - weights(sol.qp), ncol=1), 5))
	sink(type="message")
	sink()
	close(zz)
		
	# Note: equality target constraint is difficult to solve using this SOCP solver
	# because of its quality/sophistication....
	# frontier
	spec = parmaspec(S = S, benchmarkS = bS, 
			forecast = fmu - mean(B), target = NULL, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "equality",
			ineq.mat = sectors, ineq.LB = 0, ineq.UB = 0.5, budget = 1)
	control=list(abs.tol = 1e-8, rel.tol = 1e-8, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	frontsocp = parmafrontier(spec, n.points = 100, miny = NULL, maxy = NULL, 
			type = "SOCP", solver.control = control)
	
	frontqp = parmafrontier(spec, n.points = 100, miny = NULL, maxy = NULL, 
			type = "QP", solver.control = control)
	
	postscript("parma_test11b.eps", width = 12, height = 8)
	plot(frontqp[,"EV"], frontqp[,"reward"], type="l", xlab = "risk", ylab="reward",
			main = "EV Benchmark-Relative Frontier\n[QP and SOCP]", cex.main=0.9)
	points(frontsocp[,"EV"], frontsocp[,"reward"], col = 2)
	legend("topleft", c("QP","SOCP"), col = 1:2, bty="n", lty=c(1,0), pch=c(-1,1))
	dev.off()
	
	##############################################
	# Example 3: Optimize Risk-Reward (fractional SOCP problem)
	control=list(abs.tol = 1e-9, rel.tol = 1e-9, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 25)
	LB = rep(0, 15)
	UB = rep(0.5, 15)
	spec = parmaspec(S = S, forecast = fmu, target = NULL, LB = LB, UB = UB,
			risk = "EV", riskType = "optimal", budget = 1)
	
	sol.socp = parmasolve(spec, solver.control = control, type="SOCP")
	sol.qp = parmasolve(spec, type = "QP")
	
	zz <- file("parma_test11c_fractional_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP - QP")
	cat("\n--------------------------------------------\n")
	print(round(matrix(weights(sol.socp) - weights(sol.qp), ncol=1), 5))
	sink(type="message")
	sink()
	close(zz)
	
	##############################################
	# Example 4: Optimize Benchmark Risk-Reward (fractional)
	control=list(abs.tol = 1e-8, rel.tol = 1e-8, Nu=5, max.iter=3250,
			BigM.K = 4, BigM.iter = 45)
	LB = rep(0, 15)
	UB = rep(0.5, 15)
	set.seed(22)
	wb = rexp(15)
	wb = wb/sum(wb)
	B = apply(R, 1, function(x) sum(x*wb))
	bS = c(var(B), cov(B, R))
	
	
	spec = parmaspec(S = S, forecast = fmu - mean(B), 
			benchmarkS = bS, target = NULL, LB = LB, UB = UB,
			risk = "EV", riskType = "optimal", budget = 1)
	
	sol.socp = parmasolve(spec, solver.control = control, type="SOCP")
	sol.qp = parmasolve(spec, type = "QP")
	
	
	w1 = weights(sol.socp)
	w2 = weights(sol.qp)
	
	risk1 = sqrt(parmarisk(sol.socp))/parmareward(sol.socp)
	risk2 = sqrt(parmarisk(sol.qp))/parmareward(sol.qp)
	
	zz <- file("parma_test11d_benchmark_fractional_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & QP")
	cat("\n--------------------------------------------\n")
	print(round(cbind(weights(sol.socp), weights(sol.qp)), 5))
	print(cbind(risk1, risk2))
	sink(type="message")
	sink()
	close(zz)
	
	##############################################
	# Example 5: Max-Reward, Min Risk, and optimal on Frontier
	LB = rep(0, 15)
	UB = rep(1, 15)
	spec = parmaspec(S = S, forecast = fmu, target = NULL, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "equality",budget = 1)
	control=list(abs.tol = 1e-8, rel.tol = 1e-8, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	frontsocp = parmafrontier(spec, n.points = 100, miny = NULL, maxy = NULL, 
			type = "SOCP", solver.control = control)
	
	spec1 = parmaspec(S = S, forecast = fmu, target = frontsocp[1,"reward"], LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", budget = 1)
	p1 = parmasolve(spec1, type="SOCP")
	
	spec2 = parmaspec(S = S, forecast = fmu, riskB = sqrt(tail(frontsocp[,"EV"],1)), 
			LB = LB, UB = UB, risk = "EV", riskType = "maxreward", budget = 1)
	p2 = parmasolve(spec2, type="SOCP")
	
	spec3 = parmaspec(S = S, forecast = fmu, LB = LB, UB = UB, risk = "EV", 
			riskType = "optimal", budget = 1)
	p3 = parmasolve(spec3, type="SOCP")
	
	# find optimal:
	evr = sqrt(frontsocp[,"EV"])/frontsocp[,"reward"]
	idx = which(evr == min(evr))
	
	postscript("parma_test11c.eps", width = 12, height = 8)
	plot(frontsocp[,"EV"], frontsocp[,"reward"], type="l", xlab = "risk", ylab="reward",
			main = "EV Frontier\n[SOCP]", cex.main=0.9)
	points(parmarisk(p1), parmareward(p1), col = 2)
	points(parmarisk(p2), parmareward(p2), col = 4)
	points(parmarisk(p3), parmareward(p3), col = 3)
	points(frontsocp[idx,"EV"], frontsocp[idx,"reward"], col = "steelblue", pch=4)
	legend("topleft", c("Frontier","MinRisk","Optimal (Analytical)", "Optimal (Frontier)", "MaxReward"), 
			col = c(1,2,3,"steelblue",4), 
			bty="n", lty=c(1,0,0,0,0), pch=c(-1,1,1,4,1))
	dev.off()
	
	##############################################
	# Example 6: Long/Short Optimization
	LB = rep(-1, 15)
	UB = rep( 1, 15)
	control=list(abs.tol = 1e-12, rel.tol = 1e-12, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	# minrisk
	spec1 = parmaspec(S = S, forecast = fmu, target = fmub, LB = LB, UB = UB, budget = NULL,
			risk = "EV", riskType = "minrisk", targetType = "equality", leverage = 1.5)
	
	sol.socp = parmasolve(spec1, solver.control = control, type="SOCP")
	
	
	spec2 = parmaspec(scenario = R, forecast = fmu, target = fmub, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "equality", leverage = 1.5)
	
	sol.nlp = parmasolve(spec2, type="NLP")
	
	zz <- file("parma_test11e_longshort_minrisk_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & QP")
	cat("\n--------------------------------------------\n")
	print(round(cbind(weights(sol.socp), weights(sol.nlp)),5))
	sink(type="message")
	sink()
	close(zz)
	
	# optimal
	spec1 = parmaspec(S = S, forecast = fmu, LB = LB, UB = UB, budget = NULL,
			risk = "EV", riskType = "optimal", targetType = "equality", leverage = 1.5)
	
	sol.socp = parmasolve(spec1, solver.control = control, type="SOCP")
	
	
	spec2 = parmaspec(scenario = R, forecast = fmu, LB = LB, UB = UB,
			risk = "EV", riskType = "optimal", targetType = "equality", leverage = 1.5)
	
	sol.nlp = parmasolve(spec2, type="NLP")
	
	zz <- file("parma_test11f_longshort_fractional_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & NLP")
	cat("\n--------------------------------------------\n")
	print(	round(cbind(c(weights(sol.socp),sqrt(parmarisk(sol.socp)), parmareward(sol.socp)), 
							c(weights(sol.nlp), sqrt(parmarisk(sol.nlp)), parmareward(sol.nlp))),5))
	sink(type="message")
	sink()
	close(zz)

	
	##############################################
	# Example 7: minrisk with budget AND leverage
	spec1 = parmaspec(S = S, forecast = fmu, target = fmub, LB = LB, UB = UB, budget = 1,
			risk = "EV", riskType = "minrisk", targetType = "equality", leverage = 1.5)
	
	sol.socp = parmasolve(spec1, solver.control = control, type="SOCP")
	
	# NLP does not currently accept both budget and leverage (should probably changes this)
	# but we can easily supply this as a constraint:
	eqfun = function(w, optvars, uservars)
	{
		return(sum(w[optvars$widx])-1)
	}
	eqfun.jac = function(w, optvars, uservars){
		widx = optvars$widx
		fm = optvars$fm
		j = matrix(0, nrow=1, ncol=fm)
		j[1, widx] = 1
		return(j)
	}
	
	spec2 = parmaspec(scenario = R, forecast = fmu, target = fmub, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "equality", leverage = 1.5,
			eqfun = list(eqfun), eqgrad = list(eqfun.jac))
	
	sol.nlp = parmasolve(spec2, type="NLP")
	
	zz <- file("parma_test11g_longshort_fractional_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & NLP")
	cat("\n--------------------------------------------\n")
	print(round(cbind(weights(sol.socp), weights(sol.nlp)),5))
	sink(type="message")
	sink()
	close(zz)
	
	
	# optrisk with budget and leverage
	control=list(abs.tol = 1e-8, rel.tol = 1e-8, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 35)
	spec1 = parmaspec(S = S, forecast = fmu, LB = LB, UB = UB, budget = 1,
			risk = "EV", riskType = "optimal", leverage = 1.5)
	
	sol.socp = parmasolve(spec1, solver.control = control, type="SOCP")
	
	
	eqfun = function(w, optvars, uservars)
	{
		return(sum(w[optvars$widx])-w[optvars$midx])
	}
	eqfun.jac = function(w, optvars, uservars){
		widx = optvars$widx
		fm = optvars$fm
		j = matrix(0, nrow=1, ncol=fm)
		j[1, widx] = 1
		j[1, optvars$midx] = -1
		return(j)
	}
	spec2 = parmaspec(scenario = R, forecast = fmu, target = fmub, LB = LB, UB = UB,
			risk = "EV", riskType = "optimal", targetType = "equality", leverage = 1.5,
			eqfun = list(eqfun), eqgrad = list(eqfun.jac))
	
	sol.nlp = parmasolve(spec2, type="NLP")
	
	zz <- file("parma_test11h_longshort_fractional_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & NLP")
	cat("\n--------------------------------------------\n")
	print(round(cbind(c(weights(sol.socp),sqrt(parmarisk(sol.socp))/parmareward(sol.socp)), 
							c(weights(sol.nlp), sqrt(parmarisk(sol.nlp))/parmareward(sol.nlp))),5))
	sink(type="message")
	sink()
	close(zz)
	
	
	##############################################
	# Example 8: QCQP
	LB = rep(0, 15)
	UB = rep( 1, 15)
	control=list(abs.tol = 1e-12, rel.tol = 1e-12, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	
	Q = cov(tail(R, 200))
	
	spec1 = parmaspec(S = S, Q = list(Q), qB = 0.01, 
			forecast = fmu, target = fmub, LB = LB, UB = UB, budget = 1,
			risk = "EV", riskType = "minrisk", targetType = "inequality")
	
	sol.socp = parmasolve(spec1, solver.control = control, type="SOCP")
	
	# create NLP constraint
	
	ineqfun = function(w, optvars, uservars)
	{
		sqrt(w[optvars$widx]%*%uservars$Q%*%w[optvars$widx])-0.01
	}
	
	ineq.jac = function(w, optvars, uservars)
	{
		widx = optvars$widx
		fm = optvars$fm
		j = matrix(0, nrow=1, ncol=fm)
		j[1, widx] = ( (2*w[widx]%*%uservars$Q) )/(2*sqrt(w[widx]%*%uservars$Q%*%w[widx])[1])
		return(j)
	}
	uservars = list()
	uservars$Q = Q
	spec2 = parmaspec(scenario = R, forecast = fmu, target = fmub, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "inequality", budget = 1,
			ineqfun = list(ineqfun), ineqgrad = list(ineq.jac), uservars = uservars)
	
	sol.nlp = parmasolve(spec2, type="NLP")
	
	zz <- file("parma_test11i_QCQP_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & NLP")
	cat("\n--------------------------------------------\n")
	print(round(cbind(weights(sol.socp), weights(sol.nlp)),5))
	sink(type="message")
	sink()
	close(zz)
	
	# QCQP/LongShort
	LB = rep(-1, 15)
	UB = rep( 1, 15)
	control=list(abs.tol = 1e-12, rel.tol = 1e-12, Nu=4, max.iter=1250,
			BigM.K = 4, BigM.iter = 15)
	
	Q = cov(tail(R, 200))
	
	spec1 = parmaspec(S = S, Q = list(Q), qB = 0.005, 
			forecast = fmu, target = fmub, LB = LB, UB = UB, budget = NULL,
			leverage = 1.5,
			risk = "EV", riskType = "minrisk", targetType = "inequality")
	
	sol.socp = parmasolve(spec1, solver.control = control, type="SOCP")
	
	# create NLP constraint
	
	ineqfun = function(w, optvars, uservars)
	{
		sqrt(w[optvars$widx]%*%uservars$Q%*%w[optvars$widx])-0.005
	}
	
	ineq.jac = function(w, optvars, uservars)
	{
		widx = optvars$widx
		fm = optvars$fm
		j = matrix(0, nrow=1, ncol=fm)
		j[1, widx] = ( (2*w[widx]%*%uservars$Q) )/(2*sqrt(w[widx]%*%uservars$Q%*%w[widx])[1])
		return(j)
	}
	uservars = list()
	uservars$Q = Q
	spec2 = parmaspec(scenario = R, forecast = fmu, target = fmub, LB = LB, UB = UB,
			risk = "EV", riskType = "minrisk", targetType = "inequality", leverage = 1.5,
			ineqfun = list(ineqfun), ineqgrad = list(ineq.jac), uservars = uservars)
	
	sol.nlp = parmasolve(spec2, type="NLP")
	
	zz <- file("parma_test11j_QCQP_socp.txt", open="wt")
	sink(zz)
	cat("\nSOCP & NLP")
	cat("\n--------------------------------------------\n")
	print(	round(cbind(weights(sol.socp), weights(sol.nlp)),5))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	
	
}
##### MILP Tests ##### 
parma.test13 = function()
{
	tic = Sys.time()
	library(xts)
	data("etfdata")
	S = etfdata/lag(etfdata)-1
	S = na.omit(S)
	
	
	w = matrix(NA, ncol = 15, nrow = 4)
	rr = matrix(NA, ncol = 2, nrow = 4)
	rmeasure = c("MAD", "MiniMax", "CVaR", "LPM")
	for(i in 1:4){
		spec = parmaspec(scenario = S, forecast = colMeans(S), 
				target = median(colMeans(S)),
			targetType = "inequality", risk = rmeasure[i], riskType = "minrisk", 
			options = list(moment = 1, threshold = 999, alpha = 0.05),
			LB = rep(0, 15), UB = rep(0.5, 15), budget = 1, max.pos = 5)
		sol = parmasolve(spec, solver = "MILP", verbose=TRUE)
		w[i,] = weights(sol)
		rr[i,1] = parmarisk(sol)
		rr[i,2] = parmareward(sol)
	}
	# eliminate all zeros
	colnames(w) = colnames(S)
	rownames(w) = rmeasure
	idx = which(apply(w, 2, FUN = function(x) sum(x)) == 0)
	wnew = w[,-idx]
	wnew = rbind(wnew, rep(0, ncol(wnew)))
	
	postscript("parma_test13.eps", width = 12, height = 8)
	barplot(t(wnew), legend.text = TRUE, space=  0, col = c("red", "coral", "steelblue", "cadetblue", "aliceblue"),
			args.legend = list(bty = "n"), main = "MILP Portfolio (#5/15)")
	mtext(paste("parma package", sep = ""), side = 4, adj = 0, padj = -0.6, col = "grey50", cex = 0.7)
	dev.off()
	
	
	w2 = w1 = matrix(NA, ncol = 15, nrow = 4)
	rr2 = rr1 = matrix(NA, ncol = 2, nrow = 4)
	rmeasure = c("MAD", "MiniMax", "CVaR", "LPM")
	ctrl = cmaes.control()
	ctrl$options$StopOnWarnings = FALSE
	ctrl$cma$active = 1
	ctrl$options$TolFun = 1e-12
	ctrl$options$StopFitness = 0
	ctrl$options$StopOnStagnation = FALSE
	ctrl$options$DispModulo=100
	ctrl$options$Restarts = 2
	ctrl$options$MaxIter = 3000
	ctrl$options$PopSize = 300
	ctrl$options$EvalParallel = TRUE
	for(i in 1:4){
		spec = parmaspec(scenario = S, forecast = colMeans(S), 
				target = median(colMeans(S)),
				targetType = "inequality", risk = rmeasure[i], riskType = "minrisk", 
				options = list(moment = 1, threshold = 999, alpha = 0.05),
				LB = rep(0, 15), UB = rep(0.5, 15), budget = 1,
				asset.names = colnames(S))
		sol1 = parmasolve(spec, type = "LP", solver="GLPK")
		sol2 = parmasolve(spec, type = "GNLP", solver = "cmaes", solver.control = ctrl)
	
		
		w1[i,] = weights(sol1)
		rr1[i,1] = parmarisk(sol1)
		rr1[i,2] = parmareward(sol1)
		
		w2[i,] = weights(sol2)
		rr2[i,1] = parmarisk(sol2)
		rr2[i,2] = parmareward(sol2)
	}
	
	res = data.frame(MAD_LP = w1[1,], MAD_GNLP = w2[1,], MiniMax_LP = w1[2,], MiniMax_GNLP = w2[2,],
			CVaR_LP = w1[3,], CVaR_GNLP = w2[3,], LPM_LP=w1[4,], LPM_GNLP=w2[4,])
	rrd = data.frame(MAD_LP = rr1[1,], MAD_GNLP = rr2[1,], MiniMax_LP = rr1[2,], MiniMax_GNLP = rr2[2,],
			CVaR_LP = rr1[3,], CVaR_GNLP = rr2[3,], LPM_LP=rr1[4,], LPM_GNLP=rr2[4,])
	rownames(rrd) = c("Risk", "Return")
	res = rbind(res, rrd)

	# Poor Performance of GNLP
	#options(width = 120)
	zz <- file("parma_test13_MILPvGNLP.txt", open="wt")
	print(round(res,4), digits=4)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time() - tic
	return(toc)
}