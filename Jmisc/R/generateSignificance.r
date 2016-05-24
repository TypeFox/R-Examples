#' Generate t-statistics, p-value and significance from estimates and its sd. Estimates and its SD is the first and second column respectively 
#' @name generateSignificance
#' @aliases generateSignificance
#' @title Generate t-statistics, p-value and significance
#' @param x A matrix or data.frame
#' @param row_names names of row
#' @return a data.frame
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples
#' n<-1000
#' x_data<-cbind(rnorm(n,mean=0),rnorm(n,mean=1))
#' x_estimates<-cbind(apply(x_data,2,mean),apply(x_data,2,sd)/sqrt(n))
#' generateSignificance(x_estimates)
#' generateSignificance(x_estimates,row_names=c("mean0","mean1") )
generateSignificance<-function(x,row_names){
	row.names(x) = NULL 
	x<-data.frame(x)
	x$t <- x[,1] / x[,2]
	x$p_value <- 2 * (1 - pnorm( abs(x$t) ) )
	x$significance <-
		ifelse(x$p_value<0.01,'***',
			ifelse(x$p_value<0.05,'**',
				ifelse(x$p_value<0.10,'*','')
			)
		)

	names(x)<-c("estimates","sd","t","p value","significance")
	if (!missing(row_names))
		rownames(x)<-row_names
	x
}

