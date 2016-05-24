profileci <-
function(fmodel, y, zmin = -5, zmax = 5, lower = TRUE, upper = TRUE, level = 0.95, ...) 
{
	tmp <- postmode(fmodel, y, zmin, zmax, prior = function(z) return(1), ...)
	f <- function(z, ...) {
		fmodel(z, prior = function(z) return(1), ...)$post - tmp$post + qchisq(level, 1)/2
	}
	if (lower) {
		tmp[c("lower","f.lower")] <- uniroot(f, lower = zmin, upper = tmp$zeta, y = y, ...)[c("root","f.root")]
		tmp$f.lower <- tmp$f.lower + tmp$post - qchisq(level, 1)/2
	}
	else {
		tmp[c("lower","f.lower")] <- NA	
	}
	if (upper) {
		tmp[c("upper","f.upper")] <- uniroot(f, lower = tmp$zeta, upper = zmax, y = y, ...)[c("root","f.root")]
	  tmp$f.upper <- tmp$f.upper + tmp$post - qchisq(level, 1)/2
	}
	else {
		tmp[c("upper","f.upper")] <- NA
	}
	return(tmp)
}