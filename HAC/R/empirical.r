# empirical.r ############################################################################################################
# FUNCTION:            DESCRIPTION:
#  emp.copula				   Returns the values of the emprical copula for a given grid u and sample x. The estimate of the copula is based in the sample x, but is evaluated at the gird u.
#  emp.copula.self     Returns the values of the emprical copula for a given sample x. The estimate of the copula is based in the sample x and is evaluated x.
#  .emp   					   Computes the values of the empirical copula. (Internal function)
##########################################################################################################################

emp.copula = function(u, x, proc = "M", sort = "none", margins = NULL, na.rm = FALSE, ...){
	
	 n = dim(x)[1]
	 d = dim(x)[2]
	nn = dim(u)[1]
	dd = dim(u)[2]
	
	x = .margins(x, margins)
	
	if(na.rm){
		x = na.omit(x, ...)
		u = na.omit(u, ...)
    }
	
if(sort == "none"){
	.emp(u = u, x = x, proc = proc, n = n, nn = nn, d = d)}
else
if(sort == "asc"){
	sort(.emp(u = u, x = x, proc = proc, n = n, nn = nn, d = d))}
else
if(sort == "desc"){
	sort(.emp(u = u, x = x, proc = proc, n = n, nn = nn, d = d), decreasing = TRUE)}
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

emp.copula.self = function(x, proc = "M", sort = "none", margins = NULL, na.rm = FALSE, ...){switch(proc, M = emp.copula(x, x, "M", sort = sort, margins = margins, na.rm = na.rm, ...), A = emp.copula(x, x, "A", sort = sort, margins = margins, na.rm = na.rm, ...))}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.emp = function(u, x, proc, n, nn, d){
	if(proc == "M"){
		Compare = matrix(t(matrix(rep(u, n), ncol = n * d)), ncol = d, byrow = TRUE)
		Values = matrix(rep(t(x), nn), ncol = d, byrow = TRUE)		
		1/n*colSums(matrix((rowSums(Values <= Compare) == d), ncol = nn)) 
	}else{	
        loops = function(w){apply(x, 1, FUN = function(r){prod(r <= u[w,])})}
        1/n*colSums(sapply(1:nn, FUN = loops))
	}
}
