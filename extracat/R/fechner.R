quickfechner = function(x, x.type = "diss", scale = "-", path.op = "+", sym.op = "+", rescale = FALSE, exclude.zero = FALSE, check = TRUE){
	 test = FALSE
	 
	diss <- x.type %in% c("d","diss","dissimilarity")
	sim <- x.type %in% c("s","sim","similarity")
	
	if(check){
		if(sim){
			stopifnot(regmax(x))
		}else{
			stopifnot(regmin(x))	
		}	
	}	
	stopifnot(xor(diss , sim))
	
	scale.minus <- scale %in% c("-","add","additive","+")
	scale.mult <- scale %in% c("*","mult","multiplicative","exp","expected","/","div")
	
	stopifnot( xor(scale.minus, scale.mult) )
	
	max.path <- path.op %in% c("-","add","additive","+","max","maximum")
	expected.path <- path.op %in% c("*","mult","multiplicative","exp","expected")
	
	stopifnot( xor(max.path, expected.path) )
	
	vs <- as.integer( ifelse(expected.path,1,0))


	sym.none <- sym.op %in% c(FALSE,NA,"none")
	sym.sum <-  sym.op %in% c("mean","+","sum")
	sym.min <- sym.op %in% c("min")

	stopifnot( xor( sym.none, (sym.min|sym.sum)) )
	
	
	
	# check matrix
	if(check){
		stopifnot( all(x >= 0) )
	}
	
		mx <- max(abs(x))
		
		#if(mx > 1){ # CHECK ME AGAIN
			 x <- x/mx
		#}
	# perform scaling
	if(!all(diag(x) %in% c(0,1))){
		if(scale.minus){
			if(diss){
				# diss
				dx <- diag(x)
				x <- x-dx	
				if(expected.path){
					# to sim
					x <- 1-x	
				}
			}else{
				# sim
				x <- 1 - x
				dx <- diag(x)
				if(expected.path){
					# back to sim
					x <- 1 - (x-dx)	
				}else{
					# keep diss
					x <-  x-dx
				}
			}	
		}else{
			if(diss){
				#to sim
				x <- 1-x
				dx <- diag(x)
				if(max.path){
					# back to diss
					x <- 1-x/dx
				}else{
					#keep sim
					x <- x/dx	
				}
				
			}else{
				#sim
				dx <- diag(x)
				if(max.path){
					# to diss
					x <- 1-x/dx
				}else{
					#keep sim
					x <- x/dx	
				}
			}	
		}	
	}
	
	
	# x is now a dissimilarity matrix with all(diag == 0)
	if(!test){
		fx <- .Call("quickfechner",x, as.integer(dim(x)),vs, as.integer(exclude.zero))
	}else{
		fx <- .Call("quickfechner2",x, as.integer(dim(x)),vs, as.integer(exclude.zero))
	}
	
	nm <- prod(dim(x))
	px <- fx[(nm+1):(2*nm)]
	fx <- fx[-c((nm+1):(2*nm))]
	dim(fx) <- dim(px) <- dim(x)
	if(expected.path){
		fx <- 1 - fx	
	}
	if(sym.min){
		fx <- mapply( function(y,z) min(y,z), x = fx, y = t(fx))
		dim(fx) <- dim(x)
	}
	if( sym.sum ){
		fx <- (fx+t(fx))
	}
	mx2 <- max(abs(fx))
	fx <- fx/mx2
	
	attr(fx,"path.matrix") <- px
	attr(fx,"scale.factor") <- mx*mx2

	return(fx)
	
}




regmax = function(x){
	v1 <- apply(x,1,which.max)
	v2 <- apply(x,2,which.max)
	return( all(order(v1)==v2) )
}
regmin = function(x){
	v1 <- apply(x,1,which.min)
	v2 <- apply(x,2,which.min)
	return( all(order(v1)==v2 ) )
}

getpath <- function(fm, pm = NULL, from = 1, to = nrow(fm)){
	
	v <- attr(fm,"path.matrix")[from,]
	if(is.null(v)&is.null(pm)){
		simpleWarning("No path matrix.")
		return(invisible(FALSE))
	}
	if(!is.null(pm)){
		v <- pm[from,]
	}
	path <- to
	k <- to
	while(v[k] != from){
		path <- c(v[k],path)
		k <- v[k]
	}
	path <- c(from,path)
	return(path)
}

