predict.mlcm <- function(object, newdata = NULL, type = "link", ...) {
	miss <- missing(newdata)
	if(object$method == "glm"){
		if(miss){
			ans <- predict(object$obj, type = type, ...)
			} else {
#			ans <- predict(object$obj, newdata = newdata,
#				type = type, ...)
			stop("No point in using newdata for glm method\n")
			}
		} else {# if formula
		if (miss){
			if (is.null(object$whichdim)){
				 dm <- dim(object$pscale)
				 alst <- vector("list", dm[2] + 1)
				 alst[[1]] <- object$par
				 ans <- sapply(seq_len(dm[2]), function(ix){
				 	blst <-lapply(seq_len(dm[2]),
				 	function(iy){
				 		if(diag(dm[2])[ix, iy] == 1) 
				 			object$stimulus else 1
					})
				for (iy in seq_len(dm[2])) alst[[iy + 1]] <- blst[[iy]]
				do.call(object$func, alst)
				})
				 
				} else{#if a spec dim
			ans <- with(object, func(par, stimulus))
			}
		  } else { #if new data # if all dims
		  	xx <- unlist(newdata)
		  	if (is.null(object$whichdim)){
				 dm <- dim(object$pscale)
				 alst <- vector("list", dm[2] + 1)
				 alst[[1]] <- object$par
				 ans <- sapply(seq_len(dm[2]), function(ix){
				 	blst <-lapply(seq_len(dm[2]),
				 	function(iy){
				 		if(diag(dm[2])[ix, iy] == 1) 
				 			xx else 1
					})
				for (iy in seq_len(dm[2])) alst[[iy + 1]] <- blst[[iy]]
				do.call(object$func, alst)		  			})
		  		} else { # if spec dim 
		  		ans <- with(object, func(par, xx))
		  		}
		  	}	
		}
	as.vector(ans)
}