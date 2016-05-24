systemnoise <-
function(dim, nscan, type=c("gaussian","rician"), sigma, vee=1, template, verbose=TRUE){

	if(length(dim)>3){
		stop("Image space with more than three dimensions is not supported.")
	}
	if(missing(type)){
		type="gaussian"
	}

	if(type=="gaussian"){
		noise <- array(rnorm(prod(dim)*nscan, 0, sigma), dim=c(dim,nscan))
	}else{ 
	   if(type=="rician"){
		noise <- array(rrice(prod(dim)*nscan, vee=vee, sigma=sigma), dim=c(dim,nscan))
	   } else {
		stop("Specified type of system noise is unknown. Type should be gaussian or rician.")
	   }
	}

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

