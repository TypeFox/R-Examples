# pgff.R 
# Copyright (C) 2014 by Matthew Clegg

# Routines implementing the weighted symmetric estimator rho of
# Pantula, Gonzales-Farias and Fuller

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# In this module, we implement two unit root tests for autoregressive
# processes.  The first test is for the process:
#
#   X[t+1] = a_0 + (1-a_1) * X[t] + eps[t+1]
#
# where eps is an IID (normal?) random variable with mean 0.
# The first test implemented for this process is the weighted symmetric
# estimator rho_ws described in the following paper:
#
#  Pantula, S.G., Gonzalez-Farias, G., and Fuller, W.A., (1994)
#  "A Comparison of Unit Root Test Criteria,"
#  Journal of Business & Economic Statistics, 12(4), 449-459. 
#
# The second test is for the process:
#
#   X[t+1] = b_0 + b_1 * t + r[t+1]
#   r[t+1] = (1 - a_1) r[t] + eps[t+1]
#
# If 0 < a_1 < 1, then this is a trend stationary autoregressive
# process.  
#
# The second is computed by de-trending the data and calculating
# the rho_ws statistic on the de-trended data.  This changes the
# distribution of rho_ws, so a modified table of quantile values
# is used in this case.

pgff_rho_ws <- function (Y, detrend=FALSE) {
	# Calculates the weighted symmetric estimator rho_ws described in
	#   Pantula, S.G., Gonzalez-Farias, G., & Fuller, W.A. (1994).  A Comparison
	#   of Unit-Root Test Criteria.  Journal of Business & Economic Statistics,
	#   12(4), 449-459.
	
	if (detrend) {
		y <- detrend(Y)
	} else {
		y <- coredata(Y) - mean(coredata(Y))
	}
	n <- length(y)
	rho_ws_num <- sum(y[1:(n-1)]*y[2:n])
	rho_ws_den <- sum(y[2:(n-1)]^2) + (1/n)*sum(y^2)
	rho_ws <- rho_ws_num/rho_ws_den
	rho_ws
}

pgff_quantiles <- function(sample_size=100, nrep=40000, 
	q=c(0.001, 0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999),
	sd=1, detrend=FALSE) {
	# Calculates quantiles of the pgff_rho_ws function under the assumptions
	# a_0 = 0 and a_1 = 1.
	qvals <- replicate(nrep, pgff_rho_ws(rar1(sample_size, sd=sd), detrend=detrend))
	quantile(qvals, q)
}

pgff_quantile_table <- function(nrep=40000,
	q=c(0.001, 0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999),
	n=c(25, 50, 100, 250, 500, 750, 1000, 1250), sd=1, detrend=FALSE) {
	# Calculates a table of quantile values by sample size of the pgff_rho_ws function
	# under the assumption a_0=0 and a_1=1.
	df <- do.call("cbind", lapply(n, function(nv) c(nv, pgff_quantiles(nv, nrep, q, sd, detrend=detrend))))
	df <- as.data.frame(df)
	colnames(df) <- n
	df <- cbind(data.frame(quantile=c(NA,q)), df) 
	df
}

# The following table was generated using the call
#   pgff_qtab <- pgff_quantile_table()
pgff_qtab <- structure(list(quantile = c(NA, 0.001, 0.01, 0.025, 0.05, 0.1, 
0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999), `25` = c(25, 0.082689115882454, 
0.322200652341753, 0.4300782608012, 0.520399979711841, 0.612185306663716, 
0.712116046191576, 0.861136943529677, 0.959228660745895, 0.99419351964997, 
1.02068910780202, 1.04052317148161, 1.06346297844783, 1.11942982973306
), `50` = c(50, 0.500370428446805, 0.641620411151757, 0.703946341631602, 
0.750111287497596, 0.801332256205306, 0.85574636308022, 0.932882292397377, 
0.981302930412182, 0.998977952350823, 1.01139764049252, 1.02155021991231, 
1.03356784684471, 1.06177323709819), `100` = c(100, 0.717214002190848, 
0.81055602604396, 0.844201902021967, 0.871747810852699, 0.899308241082093, 
0.927336161285869, 0.966794893172991, 0.990822467384842, 0.999424524925567, 
1.00558189778188, 1.01052722350505, 1.01620574215967, 1.02952504235636
), `250` = c(250, 0.888979685474781, 0.922317230379858, 0.937387758576667, 
0.94831776910876, 0.959479998854327, 0.97085795539352, 0.986617201048855, 
0.996455921812411, 0.999931382015853, 1.00236760431612, 1.00435911200533, 
1.00672550596963, 1.01241398987593), `500` = c(500, 0.943126492866217, 
0.961041170080832, 0.968420082615165, 0.974012035193666, 0.979635473359819, 
0.985253794904594, 0.993280675653231, 0.998180932027841, 0.999993774119005, 
1.00122575726077, 1.0021963591016, 1.00338939033138, 1.00625423618206
), `750` = c(750, 0.961974674865205, 0.974267324169956, 0.979092153756704, 
0.982775020713916, 0.986379212586772, 0.990248040401925, 0.995582829639505, 
0.998816234722642, 0.999979660385506, 1.00078870332178, 1.00145695793616, 
1.00219229591648, 1.00381470246242), `1000` = c(1000, 0.972875236445692, 
0.980569545097477, 0.984132620375334, 0.987022757349866, 0.989853414242705, 
0.992739020002136, 0.996690196531679, 0.999121091358807, 0.999991427658751, 
1.00058732719492, 1.00108156273596, 1.00167295826547, 1.00292182981374
), `1250` = c(1250, 0.977154244513031, 0.984256432856506, 0.987346088912999, 
0.989595670735475, 0.99183460888197, 0.994111909177945, 0.997329831215405, 
0.999309412946087, 1.0000215085751, 1.00050827317292, 1.00088755279726, 
1.00133488753129, 1.00227779580932)), .Names = c("quantile", 
"25", "50", "100", "250", "500", "750", "1000", "1250"), row.names = c("", 
"0.1%", "1%", "2.5%", "5%", "10%", "20%", "50%", "80%", "90%", 
"95%", "97.5%", "99%", "99.9%"), class = "data.frame")

# The following table was generated using the call
#   pgff_qtab_detrended <- pgff_quantile_table(detrend=TRUE)
pgff_qtab_detrended <- structure(list(quantile = c(NA, 0.001, 0.01, 0.025, 0.05, 0.1, 
0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999), `25` = c(25, -0.127174036488274, 
0.0963632180904445, 0.199723512545942, 0.285151788768107, 0.379986406203481, 
0.491037590379484, 0.669110375477428, 0.803565647828531, 0.85796248301936, 
0.895573185059302, 0.926054182619591, 0.957995846784375, 1.0116510315628
), `50` = c(50, 0.341151398814365, 0.491653061057836, 0.559774924895833, 
0.612113520598313, 0.673411576770405, 0.736614413944516, 0.836218479416276, 
0.90700639337356, 0.935611157177896, 0.955915369012772, 0.972027218835957, 
0.987955101494032, 1.0177801775292), `100` = c(100, 0.641755151597045, 
0.731920369425943, 0.769757313288344, 0.800130335752481, 0.831968270068355, 
0.86621191405031, 0.917918249840849, 0.9551478625715, 0.969790008631989, 
0.98026755726166, 0.987957576948052, 0.9965566844378, 1.01260266853289
), `250` = c(250, 0.851801395559056, 0.889638120708458, 0.906054364034369, 
0.918590995122618, 0.932499115013238, 0.946521814382138, 0.967739742911042, 
0.982616373112442, 0.988414483083285, 0.992501362446245, 0.995650507759659, 
0.99915911067781, 1.00633294330335), `500` = c(500, 0.923827795113938, 
0.943615592890087, 0.952490504980951, 0.958907929196852, 0.965801628553788, 
0.973004761929572, 0.983790657528897, 0.991320215298679, 0.99432870851561, 
0.996424242875796, 0.998099845894082, 0.999718108334274, 1.00305893192423
), `750` = c(750, 0.94825600189366, 0.962644444384074, 0.968422325469654, 
0.972616701516689, 0.977109557989421, 0.981926099456005, 0.989161853807431, 
0.994204541796216, 0.99614917116932, 0.997569251253911, 0.998657304864418, 
0.999832355745808, 1.00203540097261), `1000` = c(1000, 0.961538036906064, 
0.971499221188825, 0.975933244659408, 0.979238453104495, 0.982647750097736, 
0.986321201151243, 0.991850949983687, 0.995643274064945, 0.997113729358009, 
0.998175637948692, 0.998994250024713, 0.999835238481858, 1.0016388797613
), `1250` = c(1250, 0.969126331856519, 0.977368467058588, 0.980936695292721, 
0.983513167385501, 0.986372929162016, 0.989241192005353, 0.993571627177825, 
0.996571083299682, 0.997766694951959, 0.998567777832407, 0.999212687519612, 
0.999869964624912, 1.0012823469108)), .Names = c("quantile", 
"25", "50", "100", "250", "500", "750", "1000", "1250"), row.names = c("", 
"0.1%", "1%", "2.5%", "5%", "10%", "20%", "50%", "80%", "90%", 
"95%", "97.5%", "99%", "99.9%"), class = "data.frame")

pgff.test <- function(Y, detrend=FALSE) {
	# Tests for a unit root of an AR(1) process using the method of
	# Pantula, Gonzales-Farias and Fuller.

	if (!exists("pgff_qtab")) {
		stop("Could not find quantile table pgff_qtab")
	}
	
    DNAME <- deparse(substitute(Y))
	STAT <- pgff_rho_ws (Y, detrend=detrend)
	
	if (detrend) {
		PVAL <- quantile_table_interpolate(pgff_qtab_detrended, length(Y), STAT)
	} else {
		PVAL <- quantile_table_interpolate(pgff_qtab, length(Y), STAT)
	}
    METHOD <- "Pantula, Gonzales-Farias and Fuller Unit Root Test"
    if (detrend) METHOD <- sprintf("%s (detrended)", METHOD)
    names(STAT) <- "rho_ws"
    alternative <- "stationary"
    structure(list(statistic = STAT, alternative = alternative, 
        p.value = PVAL, method = METHOD, data.name = DNAME), 
        class = "htest")
    
}

pgff_power <- function (a0=0, a1=0.95, trend=0, n=250, nrep=10000, p.value=0.05, detrend=FALSE) {
	# Uses simulation to estimate the power of the Pantulas, Gonzalez-Farias & Fuller
	# test with k lags
	ur_power (function(X) pgff.test(X, detrend=detrend), a0=a0, a1=a1, trend=trend,
      n=n, nrep=nrep, p.value=p.value)
}

pgff_power_table <- function (nrep=1000, p.value=0.05,
	a1=c(0.995, 0.99, 0.98, 0.97, 0.96, 0.95),
    trend=0,
	n=c(100, 250, 500, 750, 1000, 1250),
	detrend=FALSE) {
	# Constructs a table of power estimates of the Pantulas, Gonzalez-Farias & Fuller
	# test with k lags
	ur_power_table(function(X) pgff.test(X, detrend=detrend), nrep=nrep, p.value=p.value, 
        a1=a1, n=n, trend=trend)
}
