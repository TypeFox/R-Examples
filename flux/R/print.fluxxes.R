"print.fluxxes" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
	## prepare
	xt <- x$flux.table
	ab <- which(names(xt)=="all")
	nms.tab <- xt[,1:(ab-1)]
	co <- grep("CO2", names(xt))
	ch <- grep("CH4", names(xt))
	no <- grep("N2O", names(xt))
	htd <- xt[,-c(1:ab,co,ch,no)]
	out <- nms.tab
	## output
	cat("GHG flux rates and quality flags", "\n")
	if(length(co)!=0){
		CO2 <- xt[,co]
		CO2.tab <- data.frame(CO2.flags = paste(CO2[,4], CO2[,2], CO2[,3], ".", CO2[,10], sep=""), CO2.flux = round(CO2[,7], 3), CO2.us = format(paste(CO2[,6], CO2[,1], sep=" ")))
		out <- data.frame(out, CO2.tab)
	}
	if(length(ch)!=0){
		CH4 <- xt[,ch]
		CH4.tab <- data.frame(CH4.flags = paste(CH4[,4], CH4[,2], CH4[,3], ".", CH4[,10], sep=""), CH4.flux = round(CH4[,7], 3), CH4.us = format(paste(CH4[,6], CH4[,1], sep=" ")))
		out <- data.frame(out, CH4.tab)
	}
	if(length(no)!=0){
		N2O <- xt[,no]
		N2O.tab <- data.frame(N2O.flags = paste(N2O[,4], N2O[,2], N2O[,3], ".", N2O[,10], sep=""), N2O.flux = round(N2O[,7], 3), N2O.us = format(paste(N2O[,6], N2O[,1], sep=" ")))
		out <- data.frame(out, N2O.tab)
	}
	#htd <- apply(htd , 2, round, 3)
	out <- data.frame(out, htd)
	rownames(out) <- c(1:nrow(out))
	cat("\n", "flag meanings: 'nrmse r2 range . nomba'", "\n\n")
	## print to console
	print(out)
	cat("\n")
	invisible(x)
}