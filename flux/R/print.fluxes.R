"print.fluxes" <-
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
	## rl stuff
	if(is.null(x$range.lim)){
		rl.statement <- "is reported per measurement in the table"
		rl.out <- list(rl.statement, rl.statement, rl.statement)
	}
	else{
		rla <- x$range.lim
		rl.out <- x$range.lim
		for(i in c(1:length(rla))){
			rl <- rla[[i]]
			if(length(rl)!=1){
				mean.rl <- round(mean(rl))
				range.rl <- round(range(rl))
				rl.statement <- paste(mean.rl, "on average, ranging from", range.rl[1], "to", range.rl[2])
			} else {rl.statement <- paste(rl, "(global)")}
			rl.out[[i]] <- rl.statement
		}
	}
	## output
	cat("GHG flux rates and quality flags", "\n")
	if(length(co)!=0){
		cat("CO2 range limit:", rl.out[[1]], "\n")
		CO2 <- xt[,co]
		CO2.tab <- data.frame(CO2.flags = paste(CO2[,8], CO2[,6], CO2[,7], ".", CO2[,9], CO2[,10], sep=""), CO2.flux = round(CO2[,3], 3), CO2.us = format(paste(CO2[,1], CO2[,2], sep=" ")))
		out <- data.frame(out, CO2.tab)
	}
	if(length(ch)!=0){
		cat("CH4 range limit:", rl.out[[2]], "\n")
		CH4 <- xt[,ch]	
		CH4.tab <- data.frame(CH4.flags = paste(CH4[,8], CH4[,6], CH4[,7], ".", CH4[,9], sep=""), CH4.flux = round(CH4[,3], 3), CH4.us = format(paste(CH4[,1], CH4[,2], sep=" ")))
		out <- data.frame(out, CH4.tab)
	}		
	if(length(no)!=0){
		cat("N2O range limit:", rl.out[[3]], "\n")
		N2O <- xt[,no]	
		N2O.tab <- data.frame(N2O.flags = paste(N2O[,8], N2O[,6], N2O[,7], ".", N2O[,9], sep=""), N2O.flux = round(N2O[,3], 3), N2O.us = format(paste(N2O[,1], N2O[,2], sep=" ")))
		out <- data.frame(out, N2O.tab)
	}
	#htd <- apply(htd , 2, round, 3)
	out <- data.frame(out, htd)
	rownames(out) <- c(1:nrow(out))
	cat("\n", "flag meanings: 'nrmse r2 range . nomba leak'", "\n\n")
	## print to console
	print(out)
	cat("\n")
	invisible(x)
}