.fxnDeps <- function(func, package = "BAMMtools", verbose=TRUE, thorough=TRUE) {
	pstr <- paste("package:",package,sep="");
	fl <- lsf.str(pstr);
	if (!func%in%fl)
		stop(sprintf("Function %s is not in package %s",func,package));
	calledby <- sapply(fl, function(x) grep(func, body(x)));
	calls <- sapply(fl, grep, body(func));
	calledby <- names(calledby[sapply(calledby, length) > 0]);
	calls <- names(calls[sapply(calls, length) > 0]);
	if (verbose) {
		x <- sprintf("%s calls %s", func, calls);
		y <- sprintf("%s is called by %s", func, calledby);
		cat(y,sep="\n");
		cat("----\n");
		cat(x,sep="\n");
	}
	if (thorough) {
		cat("----\n");
		cat("running examples\n");
		e <- .Options[["warn"]];
		options(warn = 3);
		for (i in calls) {
			cat(i,"callee\n",sep=": ");
			x <- try(eval(call("example",topic=i,package=package,echo=FALSE,verbose=FALSE,ask=FALSE)),silent=TRUE);
			if (inherits(x,"try-error")) {
				cat(sprintf("!***! callee %s has no example !***!",i),"\n");
			}
		}
		for (i in calledby) {
			cat(i,"caller\n",sep = ": ");
			x <- try(eval(call("example",topic=i,package=package,echo=FALSE,verbose=FALSE,ask=FALSE)),silent=TRUE);
			if (inherits(x,"try-error")) {
				cat(sprintf("!***! caller %s has no example !***!",i),"\n");
			}
		}
	}
	options(warn = e);
	invisible(list(calledby = calledby, calls = calls));
}
