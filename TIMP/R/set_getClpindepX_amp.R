    setMethod("getClpindepX", signature(model = "amp"), function(model, 
        multimodel, theta, returnX, rawtheta, dind) {
	if(returnX) 
		 theta <-  getThetaCl(rawtheta, multimodel)[[dind]]
        x <- matrix(0,nrow = length(theta@amps[[1]]),
                        ncol = length(theta@amps[[1]]))
        diag(x) <- theta@amps[[1]]
        for(i in 2:length(theta@amps)) {
          xn <- matrix(0,nrow = length(theta@amps[[i]]),
                       ncol = length(theta@amps[[i]]))
          diag(xn) <- theta@amps[[i]]
          x <- rbind(x,xn)
        }
	if(returnX) 
		    x <- as.vector(x) 
	x
		    
    })
