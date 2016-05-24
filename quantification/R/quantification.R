library(car)


prep <- function(y.series, forecast.horizon, first.period, last.period, s.up, s.same, s.down) {
	if (forecast.horizon <= 0) stop(gettextf("\nForecast horizon must be greater than zero.\n"))
	if (length(s.up) != length(s.down) | length(s.same) != length(s.down)) stop(gettextf("\nThe survey responses s.up, s.same and s.down need to be of the same length. Instead they are of length: s.up = %i, s.same = %i and s.down = %i.\n",length(s.up),length(s.same), length(s.down)))
	if (length(y.series) != length(s.up)) stop(gettextf("\nThe current inflation series needs to be of the same length as the survey responses. Instead the survey results are of length %i and the series of current inflation is of length %i.\n",length(s.up),length(y.series)))
	if (first.period >= last.period) stop(gettextf("\nfirst.period (actual value: %i) must be smaller than last.period (actual value: %i).\n",first.period,last.period))
	if (last.period >= length(s.up)) stop(gettextf("\nlast.period must point to a period in the time series provided. However, last.period is %i while the time series are only of length %i.\n",last.period,length(s.up)))
	if ((last.period + forecast.horizon) > length(s.up)) stop(gettextf("\nThe period to which the last forecast relates must be within the time series provided. However, the last forecast period would be %i while the last period in the time series has index %i.\n",(last.period+forecast.horizon),length(s.up)))
	
	z <- s.up
	is.na(z[1:(first.period - 1)]) <- TRUE
	is.na(z[(last.period + 1):length(z)]) <- TRUE
	assign("survey.up", value=z, envir=parent.frame())
	z <- s.same
	is.na(z[1:(first.period - 1)]) <- TRUE
	is.na(z[(last.period + 1):length(z)]) <- TRUE
	assign("survey.same", value=z, envir=parent.frame())
	z <- s.down
	is.na(z[1:(first.period - 1)]) <- TRUE
	is.na(z[(last.period + 1):length(z)]) <- TRUE
	assign("survey.down", value=z, envir=parent.frame())
	
	assign("y.e.mean.abs", c(rep(NA, length(s.up))), envir=parent.frame())
	assign("y.e.mean.perc", c(rep(NA, length(s.up))), envir=parent.frame())
	assign("delta.y.e.mean.abs", c(rep(NA, length(s.up))), envir=parent.frame())
	assign("delta.y.e.mean.perc", c(rep(NA, length(s.up))), envir=parent.frame())
	assign("delta.y.e.sd.abs", c(rep(NA, length(s.up))), envir=parent.frame())
	assign("delta.y.e.sd.perc", c(rep(NA, length(s.up))), envir=parent.frame())
}


ra <- function(y.series, survey.up, survey.same, survey.down, forecast.horizon, first.period = 1, last.period = (length(survey.up) - forecast.horizon), distrib.type = "normal", distrib.mean = 0, distrib.sd = 1, distrib.log.location = 0, distrib.log.scale = 1, distrib.t.df = (first.period - last.period), growth.limit = NA, symmetry.error = "white", suppress.warnings=FALSE) {
	if (distrib.type != "normal" && distrib.type != "logistic" && distrib.type != "t") stop(gettextf("\n\"%s\" is not a valid distribution type. distrib.type must be either \"normal\" or \"logistic\"\" or \"t\".\n", distrib.type))
	if (symmetry.error != "white" && symmetry.error != "small.sample") stop(gettextf("\n%s is no valid value for symmetry.error. symmetry.error must either be \"white\" or \"small.sample\".\n", symmetry.error))
	prep(y.series, forecast.horizon, first.period, last.period, survey.up, survey.same, survey.down)	
	distr.func <- c(qnorm,qlogis,qt)
	par1 <- distrib.mean
	par2 <- distrib.sd
	if (distrib.type == "normal") i <- 1
	if (distrib.type == "logistic") {
		i <- 2
		par1 <- distrib.log.location
		par2 <- distrib.log.scale
	}
	if (distrib.type == "t") {
		i <- 3
		par1 <- distrib.t.df
		par2 <- 0
	}
	a <- distr.func[[i]](survey.down / (survey.up + survey.same + survey.down), par1, par2)
	b <- distr.func[[i]]((survey.down + survey.same) / (survey.up + survey.same + survey.down), par1, par2)
	
	f <- a / (a - b)
	r <- (-1) * b / (a - b)
	
	reg.abs <- lm((y.series[(first.period+forecast.horizon):(last.period+forecast.horizon)] - y.series[first.period:last.period]) ~ f[first.period:last.period] + r[first.period:last.period] - 1)
	
	i <- (first.period + forecast.horizon):(last.period + forecast.horizon)
	y.growth <- (y.series[(first.period + forecast.horizon):(last.period + forecast.horizon)] / y.series[first.period:last.period] - rep(1,(last.period - first.period + 1)))	
	i.growth <- i[abs(y.growth) > abs(growth.limit / 100)]
	if((!is.na(growth.limit)) & (length(i.growth)!=0)) {
		med <- median(y.growth, na.rm=TRUE)
		cat("\nGrowth for observations", i.growth, "exceeds limit of", growth.limit, "percent and has thus been modified to equal median growth of", med, ".\n")
		y.growth[y.growth > growth.limit] <- med
	}
		
	reg.perc <- lm(y.growth ~ f[first.period:last.period] + r[first.period:last.period] - 1)
		
	upper.limen.abs <- reg.abs$coeff[[1]]
	lower.limen.abs <- reg.abs$coeff[[2]]
	upper.limen.perc <- reg.perc$coeff[[1]]
	lower.limen.perc <- reg.perc$coeff[[2]]
	
	cat("\nEstimated upper / lower indifference limens (assuming expectations on percentage change) are: ", upper.limen.abs, " / ", lower.limen.abs, "\nEstimated upper / lower indifference limens (assuming expectations on percentage change) are: ", upper.limen.perc, " / ", lower.limen.perc, "\n")	
	
	if ((upper.limen.abs * lower.limen.abs > 0) && (!suppress.warnings)) {
		warning(gettextf("\n Estimated upper limen (%f) and lower limen (%f) are both positive or both negative in case of expectations on absolute change.\n", upper.limen.abs, lower.limen.abs))
	}
	if ((upper.limen.abs < lower.limen.abs) && (!suppress.warnings)) {
		warning(gettextf("\n Estimated upper limen (%f) is lower than estimated lower limen (%f) in case of expectations on absolute change.\n", upper.limen.abs, lower.limen.abs))
	}
	if ((upper.limen.perc * lower.limen.perc > 0) && (!suppress.warnings)) {
		warning(gettextf("\n Estimated upper limen (%f) and lower limen (%f) are both positive or both negative in case of expectations on percentage change.\n", upper.limen.perc, lower.limen.perc))
	}
	if ((upper.limen.perc < lower.limen.perc) && (!suppress.warnings)) {
		warning(gettextf("\n Estimated upper limen (%f) is lower than estimated lower limen (%f) in case of expectations on percentage change.\n", upper.limen.perc, lower.limen.perc))
	}
	
	if (symmetry.error == "white") { adj <- c("hc0") }
	else { adj <- c("hc3") }
	
	hyp.abs <- linearHypothesis(reg.abs, hypothesis.matrix = c(1,-1), test = "F", white.adjust = adj)
	hyp.perc <- linearHypothesis(reg.perc, hypothesis.matrix = c(1,-1), test = "F", white.adjust = adj)
	
	z <- c()
	is.na(z) <- TRUE
	i <- c(first.period:last.period)
	
	delta.y.e.mean.abs <- c(rep(z, first.period + forecast.horizon - 1), unname(reg.abs$fitted.values), rep(z, length(y.series) - last.period - forecast.horizon))
	delta.y.e.mean.perc <- c(rep(z, first.period + forecast.horizon - 1), unname(reg.perc$fitted.values), rep(z, length(y.series) - last.period - forecast.horizon))
	
	y.e.mean.abs[i + forecast.horizon] <- y.series[i] + delta.y.e.mean.abs[i + forecast.horizon]
	y.e.mean.perc[i + forecast.horizon] <- y.series[i] * (1 + delta.y.e.mean.perc[i + forecast.horizon])
	nob <- last.period - first.period + 1
	mae.abs <- sum(abs(y.e.mean.abs - y.series), na.rm=TRUE) / nob
	rmse.abs <- sqrt(sum((y.e.mean.abs - y.series)^2,na.rm=TRUE)) / nob
	mae.perc <- sum(abs(y.e.mean.perc - y.series), na.rm=TRUE) / nob
	rmse.perc <- sqrt(sum((y.e.mean.perc - y.series)^2,na.rm=TRUE)) / nob
	ra.result <- list(delta.y.e.mean.abs = c(delta.y.e.mean.abs), 
										delta.y.e.mean.perc = c(delta.y.e.mean.perc), 
										y.e.mean.abs = c(y.e.mean.abs), 
										y.e.mean.perc = c(y.e.mean.perc), 
										upper.limen.abs = c(upper.limen.abs),
										lower.limen.abs = c(lower.limen.abs),
										upper.limen.perc = c(upper.limen.perc),
										lower.limen.perc = c(lower.limen.perc),
										nob=c(nob),
										mae.abs = c(mae.abs),
										rmse.abs = c(rmse.abs),
										mae.perc = c(mae.perc),
										rmse.perc = c(rmse.perc),
										symmetry.abs = c(unname(hyp.abs[2,4])),
										symmetry.perc = c(unname(hyp.perc[2,4])))
	ra <- ra.result
}


cp <- function(y.series, survey.up, survey.same, survey.down, forecast.horizon, first.period = 1, last.period = (length(survey.up) - forecast.horizon), limen.type = "carlson.parkin", const.limen = 0, user.symm.limen = 0, user.upper.limen = 0, user.lower.limen = 0, correct.zero = TRUE, correct.by = 0.01, growth.limit = NA, distrib.type = "normal", distrib.mean = 0, distrib.sd = 1, distrib.log.location = 0, distrib.log.scale = 1, distrib.t.df = (last.period - first.period), suppress.warnings=FALSE) {
	prep(y.series, forecast.horizon, first.period, last.period,survey.up, survey.same, survey.down)
	if (distrib.type != "normal" && distrib.type != "logistic" && distrib.type != "t") stop(gettextf("\n\"%s\" is not a valid distribution type. distrib.type must be either \"normal\" or \"logistic\"\" or \"t\".\n", distrib.type))
	if (limen.type != "carlson.parkin" && limen.type != "constant" && limen.type != "weber.fechner" && limen.type != "symm.series" && limen.type != "asymm.series") stop(gettextf("\n\"%s\" is not a valid limen type. limen.type must either be \"carlson.parkin\" or \"constant\" or \"weber.fechner\" or \"symm.series\" or \"asymm.series\".\n", limen.type))
	if (limen.type =="constant" && const.limen == 0) stop(gettextf("\nconst.limen must be a positive value when using \"constant\" as limen.type.\n"))
	if (limen.type == "symm.series") {
		if (length(survey.up) != length(user.symm.limen)) stop(gettextf("\nuser.symm.limen must be of the same length as the survey responses. Instead user.symm.limen is of length %i and the survey responses are of length %i.\n", length(user.symm.limen), length(survey.up)))                       
	}
	if (limen.type == "asymm.series") {
		if ((length(survey.up) != length(user.upper.limen)) || (length(survey.up) != length(user.lower.limen))) stop(gettextf("\nuser.uper.limen must and user.lower.limen be of the same length as the survey responses. Instead user.upper.limen is of length %i, user.lower.limen is of length %i and the survey responses are of length %i.\n", length(user.upper.limen), length(user.lower.limen), length(survey.up)))                       
	}
	
	for(i in first.period:last.period) {
		if (!is.na(survey.up[i]) && survey.up[i] == 0) {
			if (correct.zero == TRUE) {
				survey.up[i] <- correct.by
				if (!suppress.warnings) warning(gettextf("\nsurvey.up for oberservation %i equals 0. This leads to counter-intuitive results for expected inflation. The observation has been automatically corrected to %f. Use option correct.zero=FALSE to switch off automatic correction.\n", i, correct.by))
			}
			else {
				if (!suppress.warnings) warning(gettextf("\nsurvey.up for oberservation %i equals 0. This leads to counter-intuitive results for expected inflation. Use option correct.zero=TRUE to correct this problem automatically.\n", i))
			}
		}
		if (!is.na(survey.down[i]) && survey.down[i] == 0) {
			if (correct.zero == TRUE) {
				survey.down[i] <- correct.by
				if (!suppress.warnings) warning(gettextf("\nsurvey.down for oberservation %i equals 0. This leads to counter-intuitive results for expected inflation. The observation has been automatically corrected to %f. Use option correct.zero=FALSE to switch off automatic correction.\n", i, correct.by))
			}
			else {
				if (!suppress.warnings) warning(gettextf("\nsurvey.down for oberservation %i equals 0. This leads to counter-intuitive results for expected inflation. Use option correct.zero=TRUE to correct this problem automatically.\n", i))
			}
		}
	}
	
	distr.func <- c(qnorm, qlogis, qt)
	par1 <- distrib.mean
	par2 <- distrib.sd
	if (distrib.type == "normal") i <- 1
	if (distrib.type == "logistic") {
		i <- 2
		par1 <- distrib.log.location
		par2 <- distrib.log.scale
	}
	if (distrib.type == "t") {
		i <- 3
		par1 <- distrib.t.df
		par2 <- 0
	}
	a <- distr.func[[i]](survey.down / (survey.up + survey.same + survey.down), par1, par2)
	b <- distr.func[[i]]((survey.down + survey.same) / (survey.up + survey.same + survey.down), par1, par2)
	
	for(i in first.period:last.period) {
		if (((!is.na(a[i]) && a[i] == 0) | (!is.na(b[i]) && b[i] == 0)) && (!suppress.warnings)) warning(gettextf("\nsurvey.up for observation %i is %f, survey.down is %f. Given the %s distribution, this leads expected inflation to be independent from the share of respondents expecting inflation to rise/fall.\n", i, survey.up[i], survey.down[i], distrib.type))
	}
	
	
	if (limen.type == "constant") limen.abs <- limen.perc <- c(rep(const.limen, length(survey.up)))
	if (limen.type == "carlson.parkin") {
		limen.abs <- (sum(y.series[(first.period+forecast.horizon):(last.period+forecast.horizon)], na.rm = TRUE) - sum(y.series[first.period:last.period], na.rm = TRUE)) / sum((a + b) / (a - b), na.rm = TRUE)
		limen.abs <- c(rep(limen.abs, length(survey.up)))
		
		i <- (first.period + forecast.horizon):(last.period + forecast.horizon)
		y.growth <- (y.series[(first.period + forecast.horizon):(last.period + forecast.horizon)] / y.series[first.period:last.period] - rep(1,(last.period - first.period + 1)))	
		if(!is.na(growth.limit)) {
			i.growth <- i[abs(y.growth) > abs(growth.limit / 100)]
			if(length(i.growth)!=0) {
				med <- median(y.growth, na.rm=TRUE)
				cat("\nGrowth for observations", i.growth, "exceeds limit of", growth.limit, "percent and has thus been modified to equal median growth of", med, ".\n")
				y.growth[y.growth > growth.limit] <- med			
			}
		}
		
		limen.perc <- sum(y.growth, na.rm = TRUE)/sum((a+b)/(a-b),na.rm=TRUE)
		limen.perc <- c(rep(limen.perc, length(survey.up)))		
		cat("\nEstimated Carlson-Parkin indifference limen for (assuming expectations on absolute change) is: ", limen.abs[first.period], "\nEstimated Carlson-Parkin indifference limen (assuming expectations on relative change) is: ", limen.perc[first.period], "\n")
	}
	if (limen.type == "weber.fechner") {
		gamma_abs <- (sum(y.series[(first.period + forecast.horizon):(last.period + forecast.horizon)], na.rm = TRUE) - sum(y.series[first.period:last.period], na.rm = TRUE)) / sum(y.series * (a + b) / (a - b), na.rm = TRUE)
		
		i <- (first.period + forecast.horizon):(last.period + forecast.horizon)
		y.growth <- (y.series[(first.period + forecast.horizon):(last.period + forecast.horizon)] / y.series[first.period:last.period] - rep(1,(last.period - first.period + 1)))	
		if(!is.na(growth.limit)) {
			i.growth <- i[abs(y.growth) > abs(growth.limit / 100)]
			if(length(i.growth)!=0) {
				med <- median(y.growth, na.rm=TRUE)
				cat("\nGrowth for observations", i.growth, "exceeds limit of", growth.limit, "percent and has thus been modified to equal median growth of", med, ".\n")
				y.growth[y.growth > growth.limit] <- med			
			}
		}
		gamma_perc <- (sum(y.growth, na.rm = TRUE))/sum(y.series * (a+b)/(a-b),na.rm=TRUE)
		
		limen.abs <- y.series * gamma_abs
		limen.perc <- y.series * gamma_perc
	}
	if (limen.type == "symm.series") limen.abs <- limen.perc <- user.symm.limen
	
	i <- c(first.period:last.period)		
	
	if (limen.type == "asymm.series") {
		delta.y.e.mean.abs[i + forecast.horizon] <- (user.upper.limen[i] * a[i] - user.lower.limen[i] * b[i]) / (a[i] - b[i])
		delta.y.e.mean.perc[i + forecast.horizon] <- (user.upper.limen[i] * a[i] - user.lower.limen[i] * b[i]) / (a[i] - b[i])
		limen.abs <- NA
		limen.perc <- NA
		delta.y.e.sd.abs[i + forecast.horizon] <- NA
		delta.y.e.sd.perc[i + forecast.horizon] <- NA
	}
	else {
		delta.y.e.mean.abs[i + forecast.horizon] <- limen.abs[i] * (a[i] + b[i]) / (a[i] - b[i])
		delta.y.e.mean.perc[i + forecast.horizon] <- limen.perc[i] * (a[i] + b[i]) / (a[i] - b[i])	
		delta.y.e.sd.abs[i + forecast.horizon] <- limen.abs[i] * (-2) / (a[i] - b[i])
		delta.y.e.sd.perc[i + forecast.horizon] <- limen.perc[i] * (-2) / (a[i] - b[i])	
	}
	
	y.e.mean.abs[i + forecast.horizon] <- y.series[i] + delta.y.e.mean.abs[i + forecast.horizon]
	y.e.mean.perc[i + forecast.horizon] <- y.series[i] * (1 + delta.y.e.mean.perc[i + forecast.horizon])
	nob <- last.period - first.period + 1
	mae.abs <- sum(abs(y.e.mean.abs - y.series), na.rm=TRUE) / nob
	rmse.abs <- sqrt(sum((y.e.mean.abs - y.series)^2,na.rm=TRUE)) / nob
	mae.perc <- sum(abs(y.e.mean.perc - y.series), na.rm=TRUE) / nob
	rmse.perc <- sqrt(sum((y.e.mean.perc - y.series)^2,na.rm=TRUE)) / nob
	cp.results <- list(y.e.mean.abs = c(y.e.mean.abs), 
										 y.e.mean.perc = c(y.e.mean.perc), 
										 delta.y.e.mean.abs = c(delta.y.e.mean.abs), 
										 delta.y.e.mean.perc = c(delta.y.e.mean.perc), 
										 delta.y.e.sd.abs = c(delta.y.e.sd.abs), 
										 delta.y.e.sd.perc = c(delta.y.e.sd.perc), 
										 limen.abs = c(limen.abs), 
										 limen.perc = c(limen.perc),
										 nob = c(nob),
										 mae.abs = c(mae.abs),
										 rmse.abs = c(rmse.abs),
										 mae.perc = c(mae.perc),
										 rmse.perc = c(rmse.perc))
	
	return(cp.results)
}


bal <- function(y.series, survey.up, survey.same, survey.down, forecast.horizon, first.period = 1, last.period = (length(survey.up) - forecast.horizon), growth.limit = NA, suppress.warnings=FALSE) {
	prep(y.series, forecast.horizon, first.period, last.period, survey.up, survey.same, survey.down)
	
	theta.abs <- sum(y.series[(first.period + forecast.horizon):(last.period + forecast.horizon)] - y.series[first.period:last.period], na.rm = TRUE) / sum((survey.up[first.period:last.period] - survey.down[first.period:last.period]) / (survey.up[first.period:last.period] + survey.same[first.period:last.period] + survey.down[first.period:last.period]),na.rm=TRUE)
	
	i <- (first.period + forecast.horizon):(last.period + forecast.horizon)
	y.growth <- (y.series[(first.period + forecast.horizon):(last.period + forecast.horizon)] / y.series[first.period:last.period] - rep(1,(last.period - first.period + 1)))	
	i.growth <- i[abs(y.growth) > abs(growth.limit / 100)]
	if((!is.na(growth.limit)) & (length(i.growth)!=0)) {
		med <- median(y.growth, na.rm=TRUE)
		cat("\nGrowth for observations", i.growth, "exceeds limit of", growth.limit, "percent and has thus been modified to equal median growth of", med, ".\n")
		y.growth[y.growth > growth.limit] <- med
	}
	
	theta.perc <- sum(y.growth, na.rm = TRUE) / sum((survey.up[first.period:last.period] - survey.down[first.period:last.period]) / (survey.up[first.period:last.period] + survey.same[first.period:last.period] + survey.down[first.period:last.period]), na.rm = TRUE)	
	cat("\nScaling factor for balance statistic (assuming expectations on absolute change) is: ", theta.abs, "\nScaling factor for balance statistic (assuming expectations on percentage change) is: ", theta.perc, "\n")
	
	i<-c(first.period:last.period)
	
	delta.y.e.mean.abs[i + forecast.horizon] <- theta.abs * ((survey.up[i] - survey.down[i]) / (survey.up[i] + survey.same[i] + survey.down[i]))
	delta.y.e.sd.abs[i + forecast.horizon] <- (theta.abs^2) * ((survey.up[i] + survey.down[i]) / (survey.up[i] + survey.same[i] + survey.down[i]) - ((survey.up[i] - survey.down[i]) / (survey.up[i] + survey.same[i] + survey.down[i]))^2)
	delta.y.e.mean.perc[i + forecast.horizon] <- theta.perc * ((survey.up[i] - survey.down[i]) / (survey.up[i] + survey.same[i] + survey.down[i]))
	delta.y.e.sd.perc[i + forecast.horizon] <- (theta.perc^2) * ((survey.up[i] + survey.down[i]) / (survey.up[i] + survey.same[i] + survey.down[i]) - ((survey.up[i] - survey.down[i]) / (survey.up[i] + survey.same[i] + survey.down[i]))^2)
	
	y.e.mean.abs[i + forecast.horizon] <- y.series[i] + delta.y.e.mean.abs[i + forecast.horizon]
	y.e.mean.perc[i + forecast.horizon] <- y.series[i] * (1 + delta.y.e.mean.perc[i + forecast.horizon])
	nob <- last.period - first.period + 1
	mae.abs <- sum(abs(y.e.mean.abs - y.series), na.rm=TRUE) / nob
	rmse.abs <- sqrt(sum((y.e.mean.abs - y.series)^2,na.rm=TRUE)) / nob
	mae.perc <- sum(abs(y.e.mean.perc - y.series), na.rm=TRUE) / nob
	rmse.perc <- sqrt(sum((y.e.mean.perc - y.series)^2,na.rm=TRUE)) / nob
	
	bal.results <- list(y.e.mean.abs = c(y.e.mean.abs), 
											y.e.mean.perc = c(y.e.mean.perc), 
											delta.y.e.mean.abs = c(delta.y.e.mean.abs), 
											delta.y.e.mean.perc = c(delta.y.e.mean.perc), 
											delta.y.e.sd.abs = c(delta.y.e.sd.abs), 
											delta.y.e.sd.perc = c(delta.y.e.sd.perc), 
											theta.abs = c(theta.abs), 
											theta.perc = c(theta.perc), 
											nob = c(nob),
											mae.abs = c(mae.abs),
											rmse.abs = c(rmse.abs),
											mae.perc = c(mae.perc),
											rmse.perc = c(rmse.perc))
	return(bal.results)
}