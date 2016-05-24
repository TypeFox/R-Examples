###Based on the program written for scipy (Python) by
###Ondrej Libiger and Matt Zapala
###
###Based on some of the methods presented in:
###McArdle, B.H. and Anderson, M.J. (2001).
###Fitting multivariate models to community data:
###a comment on distance-based redundancy analysis. Ecology, 82, 290-297.


dissreg <- function(formula, data, R=1000, gower=FALSE,
					squared=TRUE, permutation="dissmatrix") {

	warning("dissreg function is deprecated. It has been renamed dissmfac.")
	return(dissmfac(formula, data, R, gower, squared=squared, permutation))
}

dissmfac <- function(formula, data, R=1000, gower=FALSE,
					squared=TRUE, permutation="dissmatrix") {
		if(permutation!="dissmatrix") {
			warning("dissmfac permutation method is now always set to dissmatrix.")
		}
		warning("dissmfac function is deprecated. It internaly use dissmfacw with different default values and was keeped for backward compatibility.")
		return(dissmfacw(formula=formula, data=data, R=R, gower=gower, squared=squared))
}

