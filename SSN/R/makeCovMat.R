makeCovMat <-
function(theta, dist.hydro, a.mat, b.mat, w.matrix = NULL,
    net.zero, x.row, y.row, x.col, y.col, useTailDownWeight,
    CorModels, use.nugget, use.anisotropy, REs)
{

	nRow <- length(x.row)
	nCol <- length(x.col)

	if(is.null(net.zero)) net.zero <- matrix(1, nrow = nRow, ncol = nCol)
	V <- matrix(0, nrow = nRow, ncol = nCol)

 	# create covariance matrix component for tailup models
	npar.sofar <- 0
	if(length(grep("tailup",CorModels)) > 0){
		if(length(grep("tailup",CorModels)) > 1)
			stop("Cannot have more than 1 tailup model")
		funname <- tolower(paste(substr(unlist(strsplit(CorModels,".", 
			fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
			fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
			"tailup"] - 1], 1, 3),".tailup", sep = ""))
		tailupmod <- call(funname, dist.hydro=dist.hydro, 
			weight = w.matrix, parsil = theta[npar.sofar + 1], 
			range1 = theta[npar.sofar + 2])
		V <- V + eval(tailupmod)*net.zero
		npar.sofar <- npar.sofar + 2
	}
	# create covariance matrix component for taildown models
	if(length(grep("taildown",CorModels)) > 0){
		if(length(grep("taildown",CorModels)) > 1)
			stop("Cannot have more than 1 taildown model")
		funname <- tolower(paste(substr(unlist(strsplit(CorModels,".", 
			fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
			fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
			"taildown"] - 1], 1, 3),".taildown", sep = ""))
		taildnmod <- call(funname, dist.hydro=dist.hydro, 
			a.mat = a.mat, b.mat = b.mat, parsil = theta[npar.sofar + 1], 
			useTailDownWeight = useTailDownWeight, weight = w.matrix,
			range1 = theta[npar.sofar + 2])
		V <- V + eval(taildnmod)*net.zero
		npar.sofar <- npar.sofar + 2
	}
	# create covariance matrix componenet for Euclidean models
	if(length(grep("Euclid",CorModels)) > 0){
		if(length(grep("Euclid",CorModels)) > 1)
			stop("Cannot have more than 1 Euclidean model")
		npar.parsil <- npar.sofar + 1
		if(use.anisotropy == FALSE) {
			dist.mat <- distGeo(x.row, y.row, x.col, y.col, 
				theta[npar.sofar + 2])
			npar.sofar <- npar.sofar + 2
		}
		else {
			dist.mat <- distGeo(x.row, y.row, x.col, y.col, 
				theta[npar.sofar + 2], theta[npar.sofar + 3], 
				theta[npar.sofar + 4])
			npar.sofar <- npar.sofar + 4
		}
		funname <- paste(tolower(substr(unlist(strsplit(CorModels,".", 
			fixed = T))[(1:length(unlist(strsplit(CorModels,".", 
			fixed = T))))[unlist(strsplit(CorModels,".", fixed = T)) == 
			"Euclid"] - 1], 1, 3)),".Euclid", sep = "")
		taileumod <- call(funname, distance.matrix = dist.mat,
			parsil = theta[npar.parsil])
		V <- V + eval(taileumod)
	}

	if(length(REs)) {
		for(ii in 1:length(REs)) {
			npar.sofar <- npar.sofar + 1
			V <- V + theta[npar.sofar]*REs[[ii]]
		}
	}

	# create diagonal covariance matrix component for nugget effect
	if(use.nugget == TRUE) {
		if(nRow != nCol) stop(		
			"covariancd matrix asymmetric -- cannot use nugget")
		npar.sofar <- npar.sofar + 1
		V <- V + diag(theta[npar.sofar], nrow = nRow, ncol = nCol)
	} else if(nRow == nCol){
		V + diag(1e-6, nrow = nRow, ncol = nCol)
	}
	  
	V

}

