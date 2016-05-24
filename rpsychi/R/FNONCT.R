FNONCT <- function(x,df1,df2,prob, interval=c(0,10000), my.tol=0.000001){
		temp <- function(ncp) pf(x,df1,df2,ncp) - prob
		return(uniroot(temp, interval, tol = my.tol)$root)
}