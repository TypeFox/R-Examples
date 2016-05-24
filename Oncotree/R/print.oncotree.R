"print.oncotree" <-
function(x,...){
	prnt <- x$parent
	nevents <- length(prnt$child)-1
	cat("Oncogenetic tree from", nevents, "events\n")
	cat("Parent function:\n")
	root <- which(prnt$child=="Root")
	pfun <- cbind(prnt$child[-root], prnt$parent[-root])
	colnames(pfun) <- c("Event", "Parent")
	for (i in 1:nrow(pfun))
		cat("\t",paste(pfun[i,],collapse=" <- "),"\n")
	if (!is.null(x$eps)){
		cat("Estimated error rates: epos=",x$eps['epos'], ", eneg=", x$eps['eneg'], "\n")
	}
    invisible(x)
}

