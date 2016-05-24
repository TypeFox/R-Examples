physnoise <-
function(dim, nscan, TR, sigma, freq.heart=1.17, freq.resp=0.2, template, verbose=TRUE){

        if(length(dim)>3){
                stop("Image space with more than three dimensions is not supported.")
        }
	HB <- 2*pi*freq.heart*TR
	RR <- 2*pi*freq.resp*TR
	t <- 1:nscan

	HRdrift <- sin(HB*t) + cos(RR*t)
	sigma.HR <- sd(HRdrift)
	HRweight <- sigma/sigma.HR

	noise <- array(rnorm(prod(dim)*nscan, 0, 1), dim=c(dim, nscan)) + HRweight*HRdrift

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

