##' PCR6 rule
##'
##' PCR6 combination rule 
##'
##' @export
##' @param MassIn Matrix with \eqn{2^n} rows and \eqn{nb} columns. Parameter \eqn{n} is the number of classes (or the length of discernment frame) and \eqn{nb} is the number of experts.
##' @param TabConflict The conflict table, which can be got using the function \eqn{ConflictTable}
##' @return Two parts:
##'
##' \item{Mass}{matrix with \eqn{2^n} rows and  one column, the combined mass}
##' \item{conf}{a number, total conflict}
##' @examples
##' ## The conflict table for two experts in a discernment frame with three elements
##' TabConflict=ConflictTable(2^3,2) 
##' m1=c(0,0.4, 0.1, 0.2, 0.2, 0, 0, 0.1);
##' m2=c(0,0.2, 0.3, 0.1, 0.1, 0, 0.2, 0.1);
##' PCR6(cbind(m1,m2),TabConflict)
##'
##' @seealso \code{\link{ConflictTable}}, \code{\link{decisionDST}} 
##'
PCR6 <- function(MassIn,TabConflict){

	if (missing(TabConflict)){     # check the number of input arguments
		n=nrow(MassIn);
		nbexperts=ncol(MassIn)
		TabConflict = ConflictTable(n,nbexperts);
	}

	Conflict = ConflictPCR6(MassIn,TabConflict);

	conf=sum(Conflict);

	n=nrow(MassIn);
	m=ncol(MassIn);
	q=matrix(1,1,n);
	for (i in 1:m){
		qj=mtoq(MassIn[,i]);
		q=q*qj;
	}

	Mass=t(qtom(q))+Conflict;

	if (round(10000*conf)!=round(10000*Mass[1])){
		stop('ACCIDENT in PCR6: the calcul of conflict if false\n')
    }

	Mass[1]=0;
	re=list(Mass,conf);
	names(re)=c("Mass","conf")
    return(re)	  
} 
