tasknoise <-
function(act.image, sigma, type=c("gaussian","rician"), vee=1){
	
	if(missing(act.image)){
		stop("An activation array is required.")
	}
	if(missing(type)){
	    type <- "gaussian"
	}
	
	dim <- dim(act.image)
	if(length(dim)>4){
		stop("The activation array has more than 4 dimensions")
	}
	

	if(length(dim)==0){
	    if(type=="gaussian"){
		noise <-rnorm(length(act.image), 0, sigma)
	    } else {
		noise <- rrice(length(act.image), vee, sigma)
	    }
	} else {
	    if(type=="gaussian"){
		noise <- array(rnorm(prod(dim), 0, sigma), dim=dim)
	    } else {
		noise <- array(rrice(prod(dim), vee, sigma), dim=dim)
	    }
	}
	ix <- which(zapsmall(act.image)!=0)
	noise[-ix] <- 0
	noise <- noise*sigma/sd(noise)

	return(noise)
}
