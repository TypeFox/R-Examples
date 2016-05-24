"print.fluss" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
	## prepare
	tabl <- x$flux.table
	ab <- which(names(tabl)=="all")
	new.table <- tabl[,1:(ab-1)]
	new.table$unit <- tabl$unit
	rownames(new.table) <- c(1:nrow(new.table))
	flags <- tabl[,(ab+6):(ab+10)]
	flags[,1:3] <- apply(flags[,1:3], 2, function(x) as.logical(x)*1)
	flags <- data.frame(flags)
	if(is.logical(flags[,5])){flags[,5] <- tabl[,5]*1}
	num.vals <- apply(tabl[,(ab+3):(ab+5)], 2, function(x) round(as.numeric(x), 3))
	symp <- symnum(tabl$pv, corr=FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," "))
	## rl stuff
	rl <- x$range.lim
	if(length(rl)!=1){
		mean.rl <- round(mean(rl))
		range.rl <- round(range(rl))
		rl.statement <- paste(mean.rl, "on average, ranging from", range.rl[1], "to", range.rl[2])
	} else {rl.statement <- paste(rl, "(global)")}
	## print to console
	cat(as.character(tabl$ghg[1]), "flux rates", "\n")
	cat("Range limit:", rl.statement, "\n\n")
	print(cbind(new.table, num.vals, sig=as.vector(symp), flags, x$flux.table[,-c(1:15)]))
	cat("\n\n")
	invisible(x)
}