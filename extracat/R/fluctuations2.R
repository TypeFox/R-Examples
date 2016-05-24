fluctileXT <- function(x, hsplit = TRUE, split = "r", fun = NULL, col = NULL, gap.prop = 0.05, plotfun, plotdim = length(dim(x))-2, border = NULL, label = TRUE, add = FALSE, ...  ){
	
	dims <- as.list(dim(x))
	nd <- length(dims)
	
		#split directions init
	if(length(hsplit) == 1){
		hsplit <- rep( c(hsplit,!hsplit),nd )[1:nd] 
	}
	if(length(hsplit) < nd){
		hsplit <- rep( hsplit,nd )[1:nd] 
	}
	
	
	
	###### -------------------------------------------------------------------------- #####
	###### ----------------------------    LAYOUT ETC   ----------------------------- #####
	
	
	nx <- min(sum(hsplit*label),2)
	ny <- min(sum(!hsplit*label),2)
	
	if(is.null(border)){
		border <- 0.1
	}
	if(length(border) == 1){
		border <- 	border * c( nx/(nx+1),1/(nx+1),1/(ny+1), ny/(ny+1) )
	}
	if(length(border) == 2){
		border <- 	c( border[1] * c( nx/(nx+1),1/(nx+1)) , border[2] * c( ny/(ny+1),1/(ny+1)) )
	}

	
	
	###### -------------------------------------------------------------------------- #####
	###### -------------------------------------------------------------------------- #####
	
	
	
	
	
	
	
	hsplit <- as.list(hsplit)

	
	# prepare plotting for correct dimension
	plot <- as.list(rep(FALSE,nd))
	plot[[plotdim]] <- TRUE
	plotfun <- substitute(plotfun)
	
	print(plotfun)
	# transformations
	fun <- as.list(rep(fun,nd))[1:nd]
	
	# splitting rule
	split <- as.list(rep(split,nd))[1:nd]	
	
	# gap.prop
	gap.prop <- as.list(rep(gap.prop,nd))[1:nd]	

	# viewport names and environment
		e1 <- environment()
		ltrs <- expand.grid(letters,letters,letters)
		vpn <- paste(ltrs[,1],ltrs[,2],ltrs[,3],sep="")
		k <- 0
	
	
	# TODO: color palette
	
	#base viewport:
	if(!add){
		grid.newpage()
	}
	grandparent <- viewport(x = border[1] + (1-border[1]-border[2])/2 , y = border[3] + (1-border[3]-border[4])/2, width = 1-border[1]-border[2], height = 1-border[3]-border[4],name="base")
	
	
	dynasty <- recursile( z = x, parent = grandparent, dims = dims, hsplit = hsplit , split = split,
				 fun =  fun , col = col , gap.prop = gap.prop , plot = plot, env = e1)
	

	
	return(invisible(dynasty))
}




recursile <- function( z, parent ,dims, hsplit , split, fun , gap.prop ,  env, plot, col = NULL, ...  ){
	
		pushViewport(parent)
		nv <- dims[[1]]
		dims
	
		if(!is.null(fun[[1]])){
			z <- fun[[1]](z)
		}
		# col similar?
	
		
		
		if(plot[[1]]){
			# do the plotting
			#print(env$plotfun)
				#plotcall <- call(env$plotfun)
				#plotcall[[2]] <- z
				#plotcall[[3]] <- env$col
				#for(i in seq_along(env$args)){
				#	plotcall[[3+i]] <- env$args[[i]]
				#}
				if(!is.null(col)){
					env$plotfun[[ length(env$plotfun) + 1 ]] <- col
					names(env$plotfun)[length(env$plotfun)] <- "col"
				}
			try(
				subplot <- eval(env$plotfun, list(x = z, vp = parent))
			)
			popViewport()
			return(invisible(parent))
		}else{
			
				
		# relative or equal splits:
		
			if(split[[1]] %in% c("r","rel","relative")){
				if(is.null(dim(z))){
					w <- z/sum(z)
				}else{
					w <- apply(z,1,sum)
					w <- w/sum(w)
				}
				print(w)
			}else{
				# more options?
				w <- rep(1/nv,nv)
			}
			w <-   w * (1-gap.prop[[1]])   # rep((1-gap.prop[[1]])/nv,nv)
			x <-  (0:(nv-1)) * gap.prop[[1]]/(nv-1) + cumsum(w) - w/2
			h <- rep(1,nv)
			y <- rep(0.5,nv)
			
			children <- vpList()
			
			if(is.null(dim(z))){
				z <- as.list(z)
				
				for(i in 1:nv){
					env$k <- env$k+1
					if(!hsplit[[1]]){
						
						tmp <- viewport( y[i],1-x[i],h[i],w[i],just="centre" , name = env$vpn[env$k])
						children[[i]] <- recursile( z[[i]][[1]], parent = tmp, dims =  NULL, hsplit =  NULL , split =  NULL,
				 		fun =   NULL , gap.prop = NULL , plot = TRUE, env = env, col = col)
					}else{
						tmp <- viewport( x[i],y[i],w[i],h[i],just="centre" , name = env$vpn[env$k])
						children[[i]] <- recursile( z[[i]][[1]], parent = tmp, dims =  NULL, hsplit =  NULL , split =  NULL,
				 		fun =   NULL , gap.prop = NULL , plot = TRUE, env = env, col = col)
					}
				
					
				}
				
				
			}else{
				z <- apply(z,1,list)
				for(i in 1:nv){
				env$k <- env$k+1
				if(!hsplit[[1]]){
					tmp <- viewport( y[i],1-x[i],h[i],w[i],just="centre" , name = env$vpn[env$k])
					children[[i]] <- recursile( z[[i]][[1]], parent = tmp, dims = dims[-1], hsplit = hsplit[-1] , split = split[-1],
				 	fun =  fun[-1] , gap.prop = gap.prop[-1] , plot = plot[-1], env = env, col = col)
				}else{
					tmp <- viewport( x[i],y[i],w[i],h[i],just="centre" , name = env$vpn[env$k])
					children[[i]] <- recursile( z[[i]][[1]], parent = tmp, dims = dims[-1], hsplit = hsplit[-1] , split = split[-1],
				 	fun =  fun[-1] , gap.prop = gap.prop[-1] , plot = plot[-1], env = env, col = col)
				}
				
				
				}

			}
			
			
	}
		
		
		
		
		popViewport()
		
		return(invisible(vpTree(parent,children)))

	
	
}

simplerect <- function(x, dir = "b", just = "c", col = alpha(1,0.5), bg = "lightgrey", vp = NULL){
	fill <- col
	
	if(length(just) == 1){
		just <- rep(just,2)
	}
	if(just[1] %in% c("r","right")){
		just[1] <- "right"
		xc <- 1
	}	
	if(just[1] %in% c("l","left")){
		just[1] <- "left"
		xc <- 0
	}
	if(just[2] %in% c("t","top")){
		just[2] <- "top"
		yc <- 1
	}
	if(just[2] %in% c("b","bottom")){
		just[2] <- "bottom"
		yc <- 0
	}
	if(just[1] %in% c("c","centre", "center")){
		just[1] <- "centre"
		xc <- 0.5
	}
	if(just[2] %in% c("c","centre", "center")){
		just[2] <- "centre"
		yc <- 0.5
	}
	
	if(dir %in% c("b","both")){
        
        ht <- min(1,x)
        wt <- min(1,x)
	}
    if(dir %in% c("h","horizontal",2)){
        
        wt <- min(1,x)
        ht <- 1
    }
    if(dir %in% c("v","vertical",1)){
      
        wt <- 1
        ht <- min(1,x)
    }
	if(dir %in% c("n","none")){
   	
   	 	ht <- 1
   	 	wt <- 1
   	 	if(x==0){
   	 		fill <- "white"
   	 	}
  	}

		grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, just = "centre", default.units = "npc", 
name = NULL, gp = gpar(col = NA, fill = bg), 
draw = TRUE, vp = vp)

	grid.rect(x = xc, y = yc, width = wt, height = ht, just = just, default.units = "npc", 
name = NULL, gp = gpar(col = col, fill = fill), 
draw = TRUE, vp = vp)

	
	return(invisible(TRUE))	
}


scale1<-function(x) x/sum(x)

scale2<-function(x) x/max(x)

