


gsac <- function(x, nc = Inf, maxiter = 40, zero = TRUE, r0 = 0, force.cs = FALSE, force.rs = FALSE, resort = "complete", method = "Kendall", tau0 = 0.5, stack = "max" , clean = TRUE, clean.Is = TRUE, cutoff = -20, ... ){

	
	# x = matrix
	# zero = set sparse areas to 0 or r0*expectation
	# cs, rs = force row or column split (subtables do not overlap, rs&cs == TRUE is the TDP/SAC)
	# resort = method for reordering the subtables ("complete" matrix, individual "subtable"s, "none")
	# method = TDP method
	# cutoff = minimum acceptable residuals. rows/columns with less than ... avg residual are removed ??
	
	
	
	
	# list for the non-finished submatrices
	# indices 
	Is <- list(lapply(dim(x),function(z) 1:z))
	
	RM <- itab(x)
	R0 <- r0*RM
	
	x <- (1-r0)*x + R0
	
	clusters <- list()
	k <- 0
	
	# stop if max. number of clusters is reached
	
	while(length(Is) > 0){
		# remove redundant clusters during algorithm
		# better not?
		
		if(clean.Is & length(Is) > 2){
			ic <- getIsClusters(Is,dim(x))
			sc <- setcover(ic)
			Is <- Is[sc]
		}

		
	# reorder the current matrix
	# z/Is always refers to the original matrix x
	
	sz <- sapply(Is,function(y){
		length(y[[1]])*length(y[[2]])
	})
	##print("sizes")
	##print(sz)
	
	#id <- which.max(sz)
	Is <- Is[sz > 0]
	
	#id <- length(Is)
	if(stack %in% c(1,"first","f","fifo","FIFO","1st")){
		id <- 1
	}	
	if(stack %in% c("last","l","lifo","LIFO")){
		id <- length(Is)
	}
	if(stack %in% c("max","m","highest","h","high","size")){
		id <- which.max(sz)
	}
	if(stack %in% c("min","m","smallest","s","small")){
		id <- which.min(sz)
	}
	
	##print(id)
	
	z <- Is[[id]]
	Is <- Is[-id]
	
	##print(z)
	
	##print("check01")
	##print(any(is.na(x)))
	
if(any(sapply(z,length) == 1)){
	if(length(z[[1]]) > 1){
		z[[1]] <- z[[1]][ x[z[[1]],z[[2]]] > R0[z[[1]],z[[2]]] ]
	}
	if(length(z[[2]]) > 1){
		z[[2]] <- z[[2]][ x[z[[1]],z[[2]]] > R0[z[[1]],z[[2]]] ]
	}
	if(all(sapply(z,length) > 0)){
		clusters[[length(clusters)+1]] <- z 
	}else{
		# else the cluster is empty
		print("cluster empty!")
	}
}else{
	if(resort %in% c(FALSE,"none","n")){
		# do not reorder subtables
		x0 <- x[ z[[1]], z[[2]] , drop = FALSE]
	}
	if(resort %in% c("sub","s","subtable")){
		#reorder subtables independently
		##print(dim(xxx<-x[ z[[1]], z[[2]] ,drop=FALSE]))
		##print(table(is.na(xxx)))
		##print(any(  apply(xxx,1,sum) == 0  ))
		##print(any(  apply(xxx,2,sum) == 0  ))
				
		
		x0 <- optile(x[ z[[1]], z[[2]] ,drop=FALSE], ... )
		z[[1]] <- z[[1]][attr(x0,"orders")[[1]]]
		z[[2]] <- z[[2]][attr(x0,"orders")[[2]]]
	}
	if(resort %in% c("complete","c")){
		# reorder complete matrix
		ords <- optile(x, return.data = FALSE, ... )[[3]]
		z[[1]] <- z[[1]][order(match(z[[1]],ords[[1]]))]
		z[[2]] <- z[[2]][order(match(z[[2]],ords[[2]]))]
		x0 <- x[ z[[1]], z[[2]] ,drop=FALSE]
	}
##print("check02")	
	
	#remove zero rows / columns
	rx0 <- (x0-RM[ z[[1]], z[[2]] ,drop=FALSE]) > 1e-12
	rs <- apply(rx0,1,sum)
	cs <- apply(rx0,2,sum)
	
	if(any(rs == 0)){
		rs0 <- which(rs==0)
		z[[1]] <- z[[1]][-rs0]
	}
	if(any(cs == 0)){
		cs0 <- which(cs==0)
		z[[2]] <- z[[2]][-cs0]
	}
	


	
##print("check03")
	if(any(sapply(z,length)==1)){
		# reduce, add to clusters, next z
		# reduce subtable by all empty rows and columns
		if(length(z[[1]]) > 1){
			z[[1]] <- z[[1]][ x[z[[1]],z[[2]]] > R0[z[[1]],z[[2]]] ]
		}
		if(length(z[[2]]) > 1){
			z[[2]] <- z[[2]][ x[z[[1]],z[[2]]] > R0[z[[1]],z[[2]]] ]
		}
		# if the subtable has only 1 row and/or column accept it as a cluster
		if(all(sapply(z,length) > 0)){
			clusters[[length(clusters)+1]] <- z 
		}else{
			# else the cluster is empty
			print("cluster empty!")
		}
		##print("check04a")
	}else{
	
	x0 <- x[ z[[1]], z[[2]] ,drop=FALSE]
	##print(BCI(100*x0))
	##print(BCI(x0))
	cf <- cfluctile(x0, maxsplit = 1, plot = FALSE, tau0 = tau0, method = method)
		
	##print("check04b")	
	##	print(clusters)
		
	if(length(cf[[1]]) == 1){
		#print(z)
		#print(clusters)
			# no split >>>> reduce subtable by removing sparse rows and columns
					# rx0 <- ((x0-extracat:::itab(x0))>0)+0
					# rs <- apply(rx0,1,sum)
					# cs <- apply(rx0,2,sum)
					# if(any(rs == 0)){
						# rs0 <- which(rs==0)
						# x[z[[1]][rs0], z[[2]]] <- R0[z[[1]][rs0], z[[2]]]
						# z[[1]] <- z[[1]][-rs0]
					# }
					# if(any(cs == 0)){
						# cs0 <- which(cs==0)
						# x[z[[1]], z[[2]][cs0]] <- R0[z[[1]], z[[2]][cs0]]
						# z[[2]] <- z[[2]][-cs0]
					# }
					# remove all rows and columns with a negative avg residual
					
					# ix0 <- extracat:::itab(x0)
					# rx0 <- (x0-ix0)/sqrt(ix0)
					
					rx0 <- (x0-RM[ z[[1]], z[[2]] ,drop=FALSE])/sqrt(RM[ z[[1]], z[[2]] ,drop=FALSE])
					rm <- apply(rx0,1,mean)
					cm <- apply(rx0,2,mean)
					 if(any(rm < cutoff)){
						 rs0 <- which(rm < cutoff)
						 x[z[[1]][rs0], z[[2]]] <- R0[z[[1]][rs0], z[[2]]]
						 z[[1]] <- z[[1]][-rs0]
					 }
					 if(any(cm < cutoff )){
						 cs0 <- which(cm < cutoff)
						 x[z[[1]], z[[2]][cs0]] <- R0[z[[1]], z[[2]][cs0]]
						 z[[2]] <- z[[2]][-cs0]
					 }
		##	print("check05a")	
			
			#print(z)
			#print(clusters)
			#print("5a")	
			# if any useful rows and columns are left accept the subtable as a cluster
			
			
			if(all(sapply(z,length) > 0)){
				clusters[[length(clusters)+1]] <- z 
			}# else the cluster is empty
			##print(clusters)

	}else{
		# partition accepted: decide which subtables to keep
		##print(attr(cf,"tau.values"))
		r1 <- attr(cf, "orders")[[1]][[1]]
		r2 <- attr(cf, "orders")[[1]][[2]]
		c1 <- attr(cf, "orders")[[2]][[1]]
		c2 <- attr(cf, "orders")[[2]][[2]]
		#print(list(r1,r2,c1,c2))
		m1 <- mean( x0[r1,c2] )
		m2 <- mean( x0[r2,c1] )
		
		# the "corner" with the smallest avg is dropped
		
		
		sum1 <- sum( x0[r1,c1] )
		sum2 <- sum( x0[r2,c1] )
		sum3 <- sum( x0[r1,c2] )
		sum4 <- sum( x0[r2,c2] )
		##print(matrix(c(sum1,sum2,sum3,sum4),2,2))
		
		if(m1 < m2){
			
			z1 <- list( z[[1]] , z[[2]][c1] )
			z2 <- list( z[[1]][r2] , z[[2]] )
			
	
			
			if(force.rs){
				z1[[1]] <- z1[[1]][r1]
			}
			if(force.cs){
				z2[[2]] <- z2[[2]][c2]
			}
		
			
			if(zero){
				x[z[[1]][r1], z[[2]][c2]] <- R0[z[[1]][r1], z[[2]][c2]]
				##cat("killing ",length(z[[1]][r1]),"x", length(z[[2]][c2]))
				if(force.rs&force.cs){
					x[z[[1]][r2], z[[2]][c1]] <- R0[z[[1]][r2], z[[2]][c1]]
				}
			}
	
			
		}else{
			z1 <- list( z[[1]] , z[[2]][c2] )
			z2 <- list( z[[1]][r1] , z[[2]] )
			
			if(force.rs){
				z1[[1]] <- z1[[1]][r2]
			}
			if(force.cs){
				z2[[2]] <- z2[[2]][c1]
			}	
					
			if(zero){
				x[z[[1]][r2], z[[2]][c1]] <- R0[z[[1]][r2], z[[2]][c1]]
				##cat("killing ",length(z[[1]][r2]),"x", length(z[[2]][c1]))
				if(force.rs&force.cs){
					x[z[[1]][r1], z[[2]][c2]] <- R0[z[[1]][r1], z[[2]][c2]]
				}
			}
		}
		Is[[length(Is)+1]] <- z1
		Is[[length(Is)+1]] <- z2
		##print("check05b")
	}
	}# end ifelse dim = 1
}# end ifelse dim=1
	#print(Is)
	#print(clusters)
	#print("------------------")
	
	# cleaning section: use setcover to remove redundant clusters (if nc < Inf)
	
	tot <- length(Is) + length(clusters)
	
	if(clean & tot >= nc & length(clusters) > 1){
	##	print("c1")
		# remove redundant clusters
		ic <- getIsClusters(clusters,dim(x))
		sc <- setcover(ic)
		##print(sc)
		clusters <- clusters[sc]
		tot <- length(Is) + length(clusters)
	}
	if(clean & tot >= nc & length(Is) > 1){
	##	print("c2")
		# remove redundant Is-entries before accepting them as clusters
		ic <- getIsClusters(Is,dim(x))
		sc <- setcover(ic)
		##print(sc)
		Is <- Is[sc]
		tot <- length(Is) + length(clusters)
	}
	if(clean & tot >= nc & length(clusters)*length(Is) > 0){
	##	print("c3")
		# remove redundant clusters in clusters+Is
		ic <- getIsClusters(c(clusters,Is),dim(x))
		sc <- setcover(ic)
		##print(sc)
		sc1 <- which(sc > length(clusters))
		Is <- Is[sc[sc1] - length(clusters)]
		sc2 <- which(sc <= length(clusters))
		##print(sc2)
		clusters <- clusters[sc[sc2]]
		tot <- length(Is) + length(clusters)
	}

	if(tot >= nc){
		# no redundancies left (clean == TRUE)
		clusters <- c(clusters,Is)
		Is <- list()
		print("Maximum number of clusters (nc) reached.")
	}
	if(k >= maxiter){
		clusters <- c(clusters,Is)
		Is <- list()
		print("Maximum number of iteriations (maxiter) reached.")
	}
		
		
	
	
	}# end while
	##print("finishing")
	clusters <- rapply(clusters,f=base::sort,how="list")
	##print(clusters)
	
	# final "cleaning": remove redundant clusters
	if(clean){
		ic <- getIsClusters(clusters,dim(x))
		sc <- setcover(ic)
		clusters <- clusters[sc]
	}
	

	
	return(clusters)	
}


getIsClusters <- function(Is,dim){
	n <- length(Is)
	idv <- sapply(Is,function(z){
		M <- matrix(0,dim[1],dim[2])
		M[z[[1]],z[[2]]] <- 1
		return(as.vector(M))
	})
	return(idv)
}

setcover <- function(x, k = NULL, rat = 1, s = NULL, w = NULL, check = TRUE){
	if(is.null(dim(x))){ dim(x) <- c(length(x),1)}
	if(check){
		x0 <- x
	}
	if(is.null(w)){
		w <- rep(1,nrow(x))
	}
	w <- w/sum(w)
	
	if(is.null(s)){
		s <- sum( apply(x,1,max)*w )
	}
		
	if(is.null(k)) k <- ncol(x)
	i <- 0
	ids <- NULL
	cov <- 0
	while(i < k & cov/s < rat){
		i <- i+1
		#rs <- apply(x,2,sum)
		rs <- w %*% x
		ids <- c(ids,j<-which.max(rs))
		cov <- cov + rs[j]
		drp <- which(x[,j]==1)
		w <- w[-drp]
		x <- x[-drp,,drop=FALSE]
		if(any(dim(x) < 2)){
			i <- k
		}
	}
	#cat("coverage = ", cov/s)
	attr(ids,"coverage") <- cov/s
	
	if(check){
		# test for redundancy
		# are all entries in a cluster in at least one other cluster?
		all.ok <- FALSE
		while(!all.ok){
			cs <- apply(x0[,ids,drop=FALSE],1,sum)==1
			un <- apply(x0[,ids,drop=FALSE],2,function(z){
				table(cs,z)[2,2] == 0
			})
			if(any(un)){
				ids <- ids[-which(un)[1]]
			}else{
				all.ok <- TRUE
			}
		}
		
	
		
		
	}
	return(ids)
	
}

# cseqmat <- function(x){
	# # matrix to compute cluster sequences
	# blc<-apply(x,2,function(z) which.max(z==1))
	
	# return(outer(blc,blc,">")+0)
	
# }


# aggregate.table = function(x,compress = c(TRUE,FALSE), aggregate = sum , order = NULL){
	# nd <- length(dim(x))
	# #aggregate = as.list(rep(aggregate,nd)[1:nd])
	# compress = as.list(rep(compress,nd)[1:nd])
	
	# eq <- sapply(compress,"==",TRUE)
	
	# if(is.null(order)){
		# order <- order(eq)
	# }
	
	# x <- as.data.table(as.data.frame(x))
	
	# # setnames(x,"N","Freq")
	
	# x <- x[Freq > 0]
	
	# fi <- ncol(x)
	
	# xc <- x[,list(Freq)]
	
	# #for(i in 1:nd){
	# #	x[,i,with=FALSE] <- as.factor(unlist(x[,i,with=FALSE]))
	# #}
	
	# for(i in order){
		# if( is.logical(compress[[i]])){
			# if(compress[[i]][1] == TRUE){
				# # compress identical profiles
				# f <- paste("Freq ~ ",names(x)[i],sep="")
				# #setkeyv(x, names(x)[-c(i,fi)])
				# #nmz <- names(x)[-c(i,fi)]
				# #CM <- x[,as.list(xtabs(as.formula(f))),by=key(x)]
				# #CM <- x[,as.list(xtabs(Freq~.,data=.SD)),by=nmz]
				# a <- unlist(x[,i,with=FALSE])
				# b <- apply(x[,-1,with=FALSE],1,paste,collapse=":")
				# CM <- xtabs( x[,Freq] ~ a + b)
				
				# #tt <- t(CM[,-c(1:(nd-1)),with=FALSE])
			
				# require(amap)
				# st <- subtree( suppressWarnings(hcluster(CM)),h=0)$data
			
				# xc[,eval(quote(names(x)[i])) := as.factor(paste("C",i,st,sep=":")[unlist(x[,i,with=FALSE])])]
			# }else{
				# xc[,eval(quote(names(x)[i])) := x[,i,with=FALSE]]
			# }
			
		# }else{
			# lvl <- paste("C",i,levels(x[,i,with=FALSE]),sep=":")
			# stopifnot(length(lvl) == length(compress[[i]]))
			# st <- compress[[i]]
			# #yy <- sapply(unique(st),function(z) lvl[which(st==z)])
			# xc[,eval(quote(names(x)[i])) := paste("C",i,st,sep=":")[as.integer(x[,i,with=FALSE])]] 
		# }
		
		
		
	# }
	# #mapply(function(y,z){
	# #	if(is.logical(z)){
	# #		if(z){}
	# #	}
	# #	
	# #		
	# #},y = aggregate, z = compress)
	
# }






