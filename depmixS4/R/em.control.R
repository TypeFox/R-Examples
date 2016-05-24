em.control <- function(maxit=500,tol=1e-8,crit=c("relative","absolute"),random.start=TRUE,classification=c("soft","hard")) {
	crit <- match.arg(crit)
	classification <- match.arg(classification)
	return(list(maxit=maxit,tol=tol,crit=crit,random.start=random.start,classification=classification))
}
