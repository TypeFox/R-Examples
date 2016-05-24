##' @export
printsmcure <-
function(x,Var=TRUE, ...)
{
	if(is.null(Var)) Var=TRUE
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	cat("\nCure probability model:\n")
	if (Var) {
	b <- array(x$b,c(length(x$b),4))
	rownames(b) <- x$bnm
     	colnames(b) <- c("Estimate","Std.Error","Z value","Pr(>|Z|)")
	b[,2] <- x$b_sd
	b[,3] <- x$b_zvalue
	b[,4] <- x$b_pvalue}

 	if (!Var) {
	b <- array(x$b,c(length(x$b),1))
	rownames(b) <- x$bnm
     	colnames(b) <- "Estimate"
	}
	print(b)
	cat("\n")

	cat("\nFailure time distribution model:\n")
	if (Var) {	
	beta <- array(x$beta,c(length(x$beta),4))
	rownames(beta) <- x$betanm
	colnames(beta) <- c("Estimate","Std.Error","Z value","Pr(>|Z|)")
	beta[,2] <- x$beta_sd
	beta[,3] <- x$beta_zvalue
	beta[,4] <- x$beta_pvalue}
	if (!Var) {
	beta <- array(x$beta,c(length(x$beta),1))
	rownames(beta) <- x$betanm
     	colnames(beta) <- "Estimate"
	}
      print(beta)
	invisible(x)
}

