user.prompt <- function (x) {
 	ANSWER <- readline(paste("\nType 'y' to continue to", x,"or 'n' to quit: ",sep=" "))  
	if (substr(ANSWER, 1, 1) == "n")
         {stop("Function terminated by user.")}
	}

Rehm <- function()
{
    message("Replication of Rehm (2011)'...")
    
    cat("\n")
    cat("Loading data...\n")
    data(Rehm)
    cat("> data(Rehm)\n")	
    	
	## Specification 1
	cat("> out1 <- panelAR(NURR ~ gini, data=Rehm, panelVar='ccode', timeVar='year', autoCorr='ar1', panelCorrMethod='pcse', rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out1 <- panelAR(NURR ~ gini, data=Rehm, panelVar='ccode', timeVar='year', autoCorr='ar1', panelCorrMethod='pcse', rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out1)")
	print(summary(out1))
	
	user.prompt("loop over specifications 2-8")
	
	## Loop over specifications 2-8
	cat("> controls <- c('mean_ur', 'selfemp', 'cum_right', 'tradeunion', 'deficit', 'tradeopen', 'gdp_growth')\n")
	controls <- c('mean_ur', 'selfemp', 'cum_right', 'tradeunion', 'deficit', 'tradeopen', 'gdp_growth')
	cat("> for(c in 1:length(controls)){\n> eval(parse(text=paste(\"out <- panelAR(NURR ~ gini +\", controls[c] ,\", data=Rehm, panelVar='ccode',timeVar='year', autoCorr='ar1', panelCorrMethod='pcse', rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\",sep=\"\")))\n> print(summary(out))\n> }\n")
	for(c in 1:length(controls)){
		eval(parse(text=paste("out <- panelAR(NURR ~ gini +", controls[c] ,", data=Rehm, panelVar='ccode', timeVar='year', autoCorr='ar1', panelCorrMethod='pcse', rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)",sep="")))
		print(summary(out))
	}
	
	user.prompt("specification 9")
	
	## Specification 9
	cat("> out9 <- panelAR(NURR ~ gini + mean_ur + selfemp + cum_right + tradeunion + deficit + tradeopen + gdp_growth, data=Rehm, panelVar='ccode', timeVar='year', autoCorr='ar1', panelCorrMethod='pcse', rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out9 <- panelAR(NURR ~ gini+mean_ur+selfemp+cum_right+tradeunion+deficit+tradeopen+gdp_growth, data=Rehm, panelVar='ccode', timeVar='year', autoCorr='ar1', panelCorrMethod='pcse', rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out9)")
	print(summary(out9))
}

Rehm()