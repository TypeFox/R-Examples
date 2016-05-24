#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Not all test appear here yet.
.testattr = function(test = 1, subtest = "a"){
	
	if(is.na(match(test, 1:9))) stop("\nInvalid test number (valid range is 1:9)")
	alp = c("title", "a", "b", "c", "d", "e", "f", "g", "h")
	zt = match(subtest, alp)
	if(is.na(zt)) stop("\nInvalid subtest number (valid range is a:h)")
	
	testdesc = vector(mode = "list", length = 8)
	
	testdesc[[1]][1] = "ARFIMA Model Tests"
	testdesc[[1]][2] = "Simulated Parameter Density"
	testdesc[[1]][3] = "Fitting and Filtering"
	testdesc[[1]][4] = "Unconditional Forecasting "
	testdesc[[1]][5] = "Conditional Roll Forecast"
	testdesc[[1]][6] = "Multi- Methods"
	testdesc[[1]][7] = "Roll Method"
	testdesc[[1]][8] = "Simulation"
	testdesc[[1]][9] = "Small Benchmark"

	testdesc[[2]][1] = "Fixed and Starting Parameter Tests"
	testdesc[[2]][2] = "fGARCH[ALLGARCH]"
	testdesc[[2]][3] = "eGARCH"
	testdesc[[2]][4] = "apARCH"
	testdesc[[2]][5] = "gjrGARCH"
	testdesc[[2]][6] = "sGARCH"
	
	testdesc[[3]][1] = "Fitting and Filtering Tests"
	testdesc[[3]][2] = "fGARCH"
	testdesc[[3]][3] = "sGARCH"
	testdesc[[3]][4] = "apARCH"
	testdesc[[3]][5] = "sGARCH, iGARCH & apARCH (Filtering)"
	
	testdesc[[4]][1] = "Forecast Tests"
	testdesc[[4]][2] = "sGARCH"
	testdesc[[4]][3] = "ARFIMA-ARCHM-sGARCH"
	testdesc[[4]][4] = "ARFIMA-ARCHM-gjrGARCH"
	testdesc[[4]][5] = "Forecast Performance Measures"

	testdesc[[5]][1] = "Simulation Tests"
	testdesc[[5]][2] = "sGARCH"
	testdesc[[5]][3] = "eGARCH"
	testdesc[[5]][4] = "ARFIMA-sGARCH"

	testdesc[[6]][1] = "Path Simulation Tests"
	testdesc[[6]][2] = "sGARCH"
	testdesc[[6]][3] = "gjrGARCH"
	testdesc[[6]][4] = "eGARCH"
	testdesc[[6]][5] = "apARCH"
	testdesc[[6]][6] = "fGARCH/NAGARCH)"
	testdesc[[6]][7] = "fGARCH/NGARCH"
	testdesc[[6]][8] = "fGARCH/AVGARCH"
	
	testdesc[[7]][1] = "Roll Tests"
	testdesc[[7]][2] = "Roll Test (apARCH)"
	testdesc[[7]][3] = "Roll and Resume Test"

	testdesc[[8]][1] = "Parameter Distribution Tests"
	testdesc[[8]][2] = "fGARCH/NAGARCH (Non recursive)"
	testdesc[[8]][3] = "fGARCH/ALLGARCH (recursive)"
	testdesc[[8]][4] = "ARFIMA-GARCH Simulation/Fitting Checks"
	
	if(is.na(testdesc[[test]][zt])) stop(paste("\nInvalid subtest number (valid range is a:", alp[length(testdesc[[test]])], ")", sep=""))
	ans = paste(testdesc[[test]][1],": ", testdesc[[test]][zt], sep = "")
	return(ans)
}

rugarch.runtests = function(test = 1, subtest = "a", wdir = getwd(), cluster = NULL)
{
	tmp = try(setwd(wdir), silent = TRUE)
	if(inherits(tmp, "try-error")) stop("\nInvalid wdir...")
	tat = .testattr(test, subtest)
	cat("\n")
	cat(tat)
	cat("\n")
	tstname = paste("rugarch.test", test, subtest,"(cluster)", sep="")
	ans = eval(parse(text=tstname))
	return(ans)
}