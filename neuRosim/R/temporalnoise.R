temporalnoise <-
function(dim, nscan, sigma, rho=0.2, template, verbose=TRUE){

	if(length(dim)>3){
		stop("Image space with more than three dimensions is not supported.")
	}

	array.dim <- c(dim,rep(1,3-length(dim)),nscan)
        noise.array <- array(0, dim=array.dim)
 
	for(i in 1:array.dim[1]){
	  for(j in 1:array.dim[2]){
	    for(k in 1:array.dim[3]){
	      noise.array[i,j,k,] <- arima.sim(list(ar=rho, ma=0), sd=sigma, n=nscan)
	    }
	  }
	}

	noise <- array(c(noise.array), dim=c(dim, nscan))

        if(!missing(template)){
 		if(length(dim(template))>3){
			stop("Template should be a 2D or 3D array.")
		}
               template.time <- array(rep(template,nscan), dim=c(dim,nscan))
                ix <- which(template.time!=0)
                noise[-ix] <- 0
        }
 
       return(noise)
}

