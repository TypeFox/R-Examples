spatialnoise <-
function(dim, sigma, nscan, method=c("corr", "gammaRF", "gaussRF"), type=c("gaussian","rician"), rho=0.75, FWHM=4, gamma.shape=6, gamma.rate=1, vee=1, template, verbose=TRUE){

        if(length(dim)>3){
                stop("Image space with more than three dimensions is not supported.")
        }
	if(length(dim)==1){
		stop("Spatially correlated noise structures are not defined for vectors.")
	}
	if(missing(method)){
		method <- "corr"
	}
	if(missing(type)&&(method=="corr")){
		type <- "gaussian"
	}

	if(method=="corr"){
		if(length(dim)==2){
			noise <- array(0, dim=c(dim,nscan))
			for(z in 1:nscan){
				if(type=="gaussian"){
				  start <- array(rnorm(prod(dim), 0, sigma), dim=dim)
				}else{
				  start <- array(rrice(prod(dim), vee, sigma), dim=dim)
				}
				noise.scan <- array(0, dim=dim)
				noise.scan[1,1] <- start[1,1]
				for(i in 2:dim[1]){
					noise.scan[i,1] <-rho*noise.scan[(i-1),1] + sqrt(1-rho^2)*start[i,1]
				}
				for(j in 2:dim[2]){
					noise.scan[,j] <- rho*noise.scan[,(j-1)] + sqrt(1-rho^2)*start[,j]
				}
				for(i in 2:dim[1]){
					noise.scan[i,2:dim[2]] <- rho*noise.scan[(i-1),2:dim[2]] + sqrt(1-rho^2)*noise.scan[i,2:dim[2]]
				}
				noise[,,z] <- noise.scan
			}
		} else {
			noise <- array(0, dim=c(dim,nscan))
			for(z in 1:nscan){
				if(type=="gaussian"){
				  start <- array(rnorm(prod(dim), 0, sigma), dim=dim)
				}else{
				  start <- array(rrice(prod(dim), vee, sigma), dim=dim)
				}
				noise.scan <- array(0, dim=dim)
				noise.scan[1,1,1] <- start[1,1,1]
				for(i in 2:dim[1]){
					noise.scan[i,1,1] <- rho*noise.scan[(i-1),1,1] + sqrt(1-rho^2)*start[i,1,1]
				}
				for(j in 2:dim[2]){
					noise.scan[,j,1] <- rho*noise.scan[,(j-1),1] + sqrt(1-rho^2)*start[,j,1]
				}
				for(k in 2:dim[3]){
					noise.scan[,,k] <- rho*noise.scan[,,(k-1)] + sqrt(1-rho^2)*start[,,k]
				}
				for(i in 2:dim[1]){
					noise.scan[i,2:dim[2],2:dim[3]] <- rho*noise.scan[(i-1),2:dim[2],2:dim[3]] + sqrt(1-rho^2)*noise.scan[i,2:dim[2],2:dim[3]]
				}
				noise[,,,z] <- noise.scan
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
	}
	if(method=="gaussRF"){
		#require(AnalyzeFMRI, quietly=TRUE)
		s <- diag(FWHM^2, 3)/(8*log(2))
                if(length(dim)==2){
                        dim.RF <- c(dim,1)
                } else {
                        dim.RF <- dim
                }
		voxdim <- rep(1, length(dim.RF))
		if(!missing(template)){
			if(length(dim(template))>3){
				stop("Template should be a 2D or 3D array.")
			}
 			m <- array(ifelse(template!=0, 1, 0), dim=dim.RF)
		} else {
			m <- array(1, dim=dim.RF)
		}
		if(FWHM%%2==0){
			ksize <- FWHM + 1
		} else {
			ksize <- FWHM
		}
	
		noise <- array(0, dim=c(dim.RF,nscan))
		for(z in 1:nscan){
			noise[,,,z] <- array(c(Sim.3D.GRF(d=dim.RF,voxdim=voxdim,sigma=s,ksize=ksize,mask=m,type="field")$mat), dim=dim.RF)
		}
		noise <- array(c(noise), dim=c(dim,nscan))
	}
	if(method=="gammaRF"){
                #require(AnalyzeFMRI, quietly=TRUE)
                s <- diag(FWHM^2, 3)/(8*log(2))
                if(length(dim)==2){
                        dim.RF <- c(dim,1)
                } else {
                        dim.RF <- dim
                }
                voxdim <- rep(1, length(dim.RF))
                if(!missing(template)){
			if(length(dim(template))>3){
				stop("Template should be a 2D or 3D array.")
			}
                        m <- array(ifelse(template!=0, 1, 0), dim=dim.RF)
                } else {
                        m <- array(1, dim=dim.RF)
                }
                if(FWHM%%2==0){
                        ksize <- FWHM + 1
                } else {
                        ksize <- FWHM
                }

                noise <- array(0, dim=c(dim.RF,nscan))
                for(z in 1:nscan){
                        n <- Sim.3D.GRF(d=dim.RF,voxdim=voxdim,sigma=s,ksize=ksize,mask=m,type="field")$mat
			gamma.n <- qgamma(pnorm(c(n)), shape=gamma.shape, rate=gamma.rate)
			noise[,,,z] <- array(gamma.n, dim=dim.RF)
                }
                noise <- array(c(noise), dim=c(dim,nscan))
	}

	return(noise)
}

