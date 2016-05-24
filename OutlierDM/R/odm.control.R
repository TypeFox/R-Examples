odm.control <- function(pair.cre = 1, dist.mthd = "median", Lower = .25, Upper = .75, trans = "log2", 
						centering = TRUE, projection.type = "PCA", lbda = 1, nonlin.method = "L-BFGS-B", 
						nonlin.SS = "AsymOff", nonlin.Frank = c(2, -8, 0, 1), ncl = 2, alpha = 0.05){
## Control parameters for OutlierD
	if(pair.cre < 1) stop("pair.cre must be over 1.")
    return(list(pair.cre = pair.cre, dist.mthd = dist.mthd, Lower = Lower, Upper = Upper,
				trans = trans, centering = centering, projection.type = projection.type, lbda = lbda,
				nonlin.method = nonlin.method, nonlin.SS = nonlin.SS, nonlin.Frank = nonlin.Frank, ncl = ncl,
				alpha = alpha ))
}
