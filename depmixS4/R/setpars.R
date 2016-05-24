# depends on npar
setMethod("setpars","mix",
	function(object,values,which="pars",...) {
		if(!(length(values)==npar(object))) stop("Argument 'values' has incorrect length")
		bp <- npar(object@prior)
		if(bp > 0) {
		  switch(which,
			  "pars" = {
				  if(!all(getpars(object@prior,which=which) == values[1:bp])) {
					  object@prior=setpars(object@prior,values[1:bp],which=which)
					  # recompute init probabilities
					  object@init <- dens(object@prior)
				  }
			  },
			  "fixed" = {
				  object@prior <- setpars(object@prior,values[1:bp],which=which)
			  }
		  )
		  bp <- bp+1
		  values <- values[bp:npar(object)]
		}
		if(class(object)=="depmix"|class(object)=="depmix.fitted") {
			for(i in 1:object@nstates) {
				bp <- npar(object@transition[[i]])
				if(bp > 0) {
				  switch(which,
					  "pars"= {
						  if(!all(getpars(object@transition[[i]]) == values[1:bp])) {
							  object@transition[[i]] <- setpars(object@transition[[i]],values[1:bp])
							  # recompute transition densities if pars have changed
							  object@trDens[,,i] <- dens(object@transition[[i]])
						  }
					  },
					  "fixed" = {
						  object@transition[[i]] <- setpars(object@transition[[i]],values[1:bp],which="fixed")
					  }
				  )
				  bp <- bp+1
				  values <- values[bp:length(values)]
				}
			}
		}
		for(i in 1:object@nstates) {
			for(j in 1:object@nresp) {
				bp <- npar(object@response[[i]][[j]])
				if(bp > 0) {
				  switch(which,
					  "pars" = {
						  if(!all(getpars(object@response[[i]][[j]]) == values[1:bp])) {
							  object@response[[i]][[j]] <- setpars(object@response[[i]][[j]],values[1:bp])
							  # recompute observation densities if pars have changed
							  object@dens[,j,i] <- dens(object@response[[i]][[j]])
						  }
					  },
					  "fixed" = {
						  object@response[[i]][[j]] <- setpars(object@response[[i]][[j]],values[1:bp],which="fixed")
					  }
				  )	
				  bp <- bp+1
				  values <- values[bp:length(values)]
				}
			}
		}			
		return(object)
	}
)
