
getIs <- function(biclust, dim,nstart=20, solver = "nn", adjust.dist = TRUE){
		
		if(inherits(biclust,"Biclust")){
		rxn <- biclust@RowxNumber
		nxc <- biclust@NumberxCol
		N <- biclust@Number
		
		Is <- list()
		for(i in 1:N){
			Is[[i]] <- list(  which(rxn[,i]), which(nxc[i,]))
		}
		}else{
			Is <- biclust	
		}
		Is <- OBC(Is,dim,nstart=nstart, solver = solver, adjust.dist = adjust.dist)
		return(Is)
}


getIs2 <- function(bic, dim, nstart=20, solver = "nn", cpr = FALSE, cpc = TRUE, adjust.dist = FALSE){
	
	if(inherits(bic,"Biclust")){
		rxn <- bic@RowxNumber
		nxc <- bic@NumberxCol
		N <- bic@Number
		
		Is <- list()
		for(i in 1:N){
			Is[[i]] <- list(  which(rxn[,i]), which(nxc[i,]))
		}
		}else{
			Is <- bic	
		}
	ns <- length(Is)
	
	n <- dim[1]
	m <- dim[2]
	
	M <- array(0, dim=c(n,m,ns))
	for(i in 1:ns){
		M[ Is[[i]][[1]], Is[[i]][[2]] ,i] <- 1
	}
	if(cpc){
		V <- t(apply(M,2,as.vector))
		V2 <- subtable(V,1:ncol(V))
		cFreq <- V2$Freq
		V2 <- as.matrix(V2[,-ncol(V2)])
		V2 <- optME(V2,dims=1, nstart = nstart, solver = solver, adjust.dist = adjust.dist)
				
		# match
		orig <- apply(V,1, function(z){
			paste(z,collapse = ":")
		})
		new <- apply(V2,1, function(z){
			paste(z,collapse = ":")
		})
		cord <- order(match(orig,new))
		
	}else{
		cord <- optME(M,dims = 2, solver = solver, nstart= nstart, return.table=FALSE, adjust.dist = adjust.dist)[[2]]
	}
	
	
	if(cpr){
		V <- t(apply(M,1,as.vector))
		V2 <- subtable(V,1:ncol(V))
		rFreq <- V2$Freq
		V2 <- as.matrix(V2[,-ncol(V2)])
		V2 <- optME(V2,dims=1, nstart = nstart, solver = solver, adjust.dist = adjust.dist)
				
		# match
		orig <- apply(V,1, function(z){
			paste(z,collapse = ":")
		})
		new <- apply(V2,1, function(z){
			paste(z,collapse = ":")
		})
		rord <- order(match(orig,new))
		
	}else{
		rord <- optME(M,dims = 1, solver = solver, nstart= nstart, return.table=FALSE, adjust.dist = adjust.dist)[[1]]
	}
	
	

	Is2 <- lapply(Is,function(z){
		 ret <- list(
				sort( which( rord %in% z[[1]]) ),
				sort( which( cord %in% z[[2]]) )
		   )
		   return(ret)
		})
	attr(Is2,"orders") <- list(rord,cord)
	
	
	return(Is2)
}


heattile <- function(x, biclust = NULL, Is = NULL, shape = "r", fluct = FALSE,
    gap.prop = 0.0, border = c(0.05,0.03,0.03,0.05), label = c(TRUE, FALSE), lab.opt = list(abbrev = 24, 
        lab.cex = 1, rot = 0), bg.col = "lightgrey", sym = FALSE, breaks = 20+10*sym, clust.col = NULL, clust.palette = "rgb", hm.palette = "div", clust.col.opt = list(), hm.col.opt = list(revert=TRUE)){
      
      # if biclust is NULL then check Is.
        
 if(!is.null(biclust)){
 		rxn <- biclust@RowxNumber
		nxc <- biclust@NumberxCol
		N <- biclust@Number
		
		Is <- list()
		for(i in 1:N){
			Is[[i]] <- list(  which(rxn[,i]), which(nxc[i,]))
		}
 	}else{
 		if(is.null(Is)){
 			# if Is is also NULL then perhaps x has Is as an attribute from sortMEbic
 			Is <- attr(x, "Is")
 		}
  		N <- length(Is)
 	} 
 	if(is.null(attr(Is,"orders")) & N > 0){
 		# needs to be checked ... 
 		Is <- rapply(Is,sort,how="list")
 	}
 	
 	
 if(N > 0){		
 	# cluster colors
 	if(is.null(clust.col)){
 		if(is.null(attr(Is,"colv"))){
 			clust.col <- getcolors(N,clust.palette,col.opt = clust.col.opt) 
 		}else{
 			clust.col <- attr(Is,"colv")
 		}
 	}
 }	
 
 	if(length(breaks) == 1){
 		if(sym){
 			maxv <- max(abs(x))
 			a <- -maxv
 			b <- maxv
 		}else{
 			a <- min(x)
 			b <- max(x)
 		}
 		breaks <- seq(a-1e-12,b+1e-12, (b-a+3e-12)/breaks)
 	}
 	nbreaks <- length(breaks)
 
 	# heatmap colors
 	hm.col <- getcolors(nbreaks,hm.palette,col.opt = hm.col.opt) 
 	
 	# color ids	
 
 	its <- as.integer(cut(x,breaks=breaks))
	
	#color matrix
	

	CM <- sapply( its, function(z) hm.col[z] )
	dim(CM) <- dim(x)

	
	
	
	# reorder x
	if(!is.null(attr(Is,"orders"))){
			o1 <- attr(Is,"orders")[[1]]#rev(attr(Is,"orders")[[1]])
			o2 <- attr(Is,"orders")[[2]]
	}else{
		o1 <- 1:nrow(x) #o1 <- nrow(x):1 
		o2 <- 1:ncol(x)
	} 
	if( fluct ){
		#o1 <- rev(o1) # why??? CHECK!
	}
	x <- x[o1, o2]
	CM <- CM[o1, o2]
	if(fluct){
		IM <- abs(x)
	}else{
		IM <- matrix(1,nrow=nrow(x),ncol=ncol(x) )
	}
		rownames(IM) <- rownames(x)
		colnames(IM) <- colnames(x)
		
	v1<-fluctile(IM,lab.opt=lab.opt, border = border, bg.col = bg.col, tile.col = CM, gap.prop = gap.prop, shape = shape, label = label)
		if(N > 0){
		Iss <- splitset(Is)
		colorit(Iss, dim = dim(IM),col = clust.col, vp = v1)
		}
	return(invisible(TRUE))

 
 }

 	      




        
 






OBC <- function(Is,dim, nstart=20, solver = "nn", adjust.dist = TRUE){
	
	
	ns <- length(Is)
	
	n <- dim[1]
	m <- dim[2]
	
	M <- array(0, dim=c(n,m,ns))
	for(i in 1:ns){
	
		M[ Is[[i]][[1]], Is[[i]][[2]] ,i] <- 1
	}

	M2 <- optME(M+1, nstart=nstart, solver = solver, adjust.dist = adjust.dist)
	ords <- attr(M2,"orders")
	
	Is2 <- lapply(Is,function(z){
		
		   ret <- list(
				sort( which( ords[[1]] %in% z[[1]]) ),
				sort( which( ords[[2]] %in% z[[2]]) )
		   )
		   return(ret)
		})
	attr(Is2,"orders") <- ords
	
	attr(Is2,"colv") <- attr(Is,"colv")#[ords[[3]]]
	return(Is2)
}


splitset <- function(ii){
	
	rr <- list()
	for(i in 1:length(ii)){
				
		rr[[i]]<-lapply(ii[[i]],function(z){
				
			   d <- diff(z)
			   cp <- which(d > 1)
					
			ret <- list()
			   for(j in seq_along(cp)){
					ret[[j]] <- z[1:cp[j]]
					z <- z[-c(1:cp[j])]
					cp <- cp-cp[j]
			   }
				
			   ret[[length(ret)+1]] <- z
			   ret[[length(ret)+1]] <- 0
			return(ret)
		})
		
		
	}
return(rr)

}


colorit = function(ii,dim,col,vp){
	
	ns <- length(ii)
	
	for(i in 1:ns){
		for(j in 1:( length(ii[[i]][[1]])-1 ) ){
			for(k in 1:( length(ii[[i]][[2]])-1 ) ){
				y0 <- (dim[1] - min(ii[[i]][[1]][[j]]) +1)/dim[1]
				y1 <- (dim[1] - max(ii[[i]][[1]][[j]]))/dim[1]
				x0 <- (min(ii[[i]][[2]][[k]])-1)/dim[2]
				x1 <- max(ii[[i]][[2]][[k]])/dim[2]
				grid.rect(x0, y0, x1-x0, y1-y0, gp = gpar(fill = alpha(col[i],0.3), col = col[i], 
														  lwd = 2, lty = 1), just = c("left", "bottom"),vp=vp)
			}
		}
	}
return(invisible(TRUE))
}





