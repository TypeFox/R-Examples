LPmerge <- function(Maps,max.interval=1:3,weights=NULL) {
	n.maps <- length(Maps)	
	if (n.maps < 2) {
		print("Error.  Must have at least two maps.")
        stop()
	}
	if (is.null(weights)) {weights <- rep(1,n.maps)}
	stopifnot(length(weights)==n.maps)

	map.names <- attributes(Maps)$names
	if (is.null(map.names)) {map.names <- 1:n.maps}
	
	#sort and round maps and convert factors to characters
	num.mark <- rep(NA,n.maps)
	num.unique.mark <- rep(NA,n.maps)
	for (i in 1:n.maps) {
		map <- Maps[[i]]
		map <- map[order(map[,2]),]
		map[,2] <- round(map[,2],2)  #round to 0.01 cM
		map[,1] <- as.character(map[,1])
		num.unique.mark[i] <- length(unique(map[,1]))
		num.mark[i] <- nrow(map)
		Maps[[i]] <- map
	}
	errs <- which(num.unique.mark!=num.mark)
	if (length(errs)>0) {
		print("Error. Redundant markers present in following maps:")
		print(paste(map.names[errs],collapse=" "))
		stop()
	}
	
	map <- Maps[[1]]
	markers <- map[,1]
	bins <- unique(map[,2])
	mark.bins <- match(map[,2],bins)
	for (i in 2:n.maps) {
		map <- Maps[[i]]
		mark.i <- map[,1]
		bins.i <- unique(map[,2])
		mark.bins.i <- match(map[,2],bins.i)
		
		mo <- length(markers)
		m <- length(union(mark.i,markers))
		new.mark.bins <- rep(0,m)
		shared <- which(is.element(markers,mark.i))
		not.shared <- setdiff(1:mo,shared)
		shared.i <- match(markers[shared],mark.i)
		not.shared.i <- setdiff(1:length(mark.bins.i),shared.i)
		if (length(shared)>0) {
			new.mark.bins[shared] <- paste(mark.bins[shared],mark.bins.i[shared.i],sep=".")
		}
		if (length(not.shared)>0) {
			new.mark.bins[not.shared] <- paste(mark.bins[not.shared],0,sep=".")
		}
		if (length(not.shared.i)>0) {
			if (m > mo) {
				new.mark.bins[(mo+1):m] <- paste(0,mark.bins.i[not.shared.i],sep=".")
			}
			markers <- c(markers,mark.i[not.shared.i])
		}
		mark.bins <- new.mark.bins
	}
	
	n.mark <- length(markers)
	bins <- unique(mark.bins)
	n.bin <- length(bins)
	bin.list <- list()
	
	for (j in 1:n.bin) {
		bin.list[[j]] <- markers[which(mark.bins==bins[j])]	
	}
	
	constraints <- numeric(0)
	for (i in 1:n.maps) {
		map <- Maps[[i]]
		map <- data.frame(bin=match(mark.bins[match(map[,1],markers)],bins),pos=map[,2])
		uniq <- unique(map$bin)
		map <- map[match(uniq,map$bin),]
		m <- nrow(map)

		#construct constraints
		j <- which(map[1,2]==map[,2])
		while (max(j) < m) {
			k <- which(map[max(j)+1,2]==map[,2])
			d <- map[k[1],2] - map[j[1],2]
			z <- expand.grid(j,k)
			constraints <- rbind(constraints,t(apply(z,1,function(x){return(c(i,map[x,1]))})))
			j <- k
		}		
	}
	
	n.constraint <- nrow(constraints)
	print(paste("# markers:",n.mark))
	print(paste("# bins:",n.bin))
	print(paste("# constraints:",n.constraint))

	#generate maxFS subsystem by removing constraints
	A <- Matrix(0,nrow = 0,ncol=n.bin,sparse=TRUE)
	for (i in 1:n.constraint) {
		v <- rep(0,n.bin)
		v[constraints[i,2]] <- -1
		v[constraints[i,3]] <- 1
		A <- rBind(A,v)
	}
	
	#######
	maxFS <- function(A) {
		#Algorithm 1 in Chinneck (2001)
		print("Finding maximum feasible subsystem.")
		n.constraint <- nrow(A)
		n.mark <- ncol(A)
		B <- cBind(A,Diagonal(n.constraint))
				
		#Step1
		f <- c(rep(0,n.mark),rep(1,ncol(B)-n.mark))
		ans <- Rglpk_solve_LP(f,B,rep(">=",nrow(B)),rep(1,nrow(B)))
		if (ans$status != 0) {
			stop("Error in LP solver.")	
		} else {
			if (ans$optimum < 1e-6) {
				return(integer(0))  #system is feasible
			}
			elastic.var <- ans$solution[n.mark+1:nrow(B)]
			sorted <- sort(elastic.var,decreasing=TRUE,index.return=TRUE)
			HoldSet <- sorted$ix[which(sorted$x > 1e-4)]
			eliminate <- integer(0)
			if (length(HoldSet)==1) {
				return(HoldSet)
			} else {

				#Step2
				repeat {
					candidates <- HoldSet
					min.SINF <- Inf
					for (i in 1:length(candidates)) {
						#delete constraint
						B <- cBind(A[-c(eliminate,candidates[i]),],Diagonal(n.constraint - length(eliminate)-1))
						constraint.id <- (1:n.constraint)[-c(eliminate,candidates[i])]
						f <- c(rep(0,n.mark),rep(1,ncol(B)-n.mark))
						ans <- Rglpk_solve_LP(f,B,rep(">=",nrow(B)),rep(1,nrow(B)))
						if (ans$status != 0) {
							stop("Error in LP solver.")	
						} else {
							if (ans$optimum < 1e-6) {
								eliminate <- c(eliminate,candidates[i])
								return(eliminate)
							} else {
								if (ans$optimum < min.SINF) {
									winner <- candidates[i]
									min.SINF <- ans$optimum
									elastic.var <- ans$solution[n.mark+1:nrow(B)]
									sorted <- sort(elastic.var,decreasing=TRUE,index.return=TRUE)
									HoldSet <- constraint.id[sorted$ix[which(sorted$x > 1e-4)]]
									if (length(HoldSet)==1) {
										next.winner <- HoldSet
									} else {
										next.winner <- NULL
									}
								}		
							}
						}
					} #for i 	
					
					#Step 3
					eliminate <- c(eliminate,winner)
					if (!is.null(next.winner)) {return(c(eliminate,next.winner))}
				} #repeat
			}
		}
	}
	#########
	
	eliminate <- maxFS(A)
	n.bad <- length(eliminate)
	if (n.bad>0) {
		print("Eliminated following constraints to resolve marker order conflicts: ")
		for(i in 1:n.bad) {
			mark1 <- paste(bin.list[[constraints[eliminate[i],2]]],collapse=" ")
			mark2 <- paste(bin.list[[constraints[eliminate[i],3]]],collapse=" ")
			print(paste("Map ",map.names[constraints[eliminate[i],1]],": ",mark1," < ",mark2,sep=""))
		}
		A <- A[-eliminate,]
	} else {
		print("Linkage maps had no ordering conflicts.")
	}

	n.composite.maps <- length(max.interval)
	
	result <- list()
	
	for (p in 1:n.composite.maps) {
	
		error.terms <- numeric(0)
		n.terms <- rep(0,n.maps)

		for (i in 1:n.maps) {
			map <- Maps[[i]]
			map <- data.frame(bin=match(mark.bins[match(map[,1],markers)],bins),pos=map[,2])
			uniq <- unique(map$bin)
			map <- map[match(uniq,map$bin),]
			m <- nrow(map)
	
			#make error terms
			for (q in 1:max.interval[p]) {
				for (j in 1:(m-q)) {
					n.terms[i] <- n.terms[i]+1
					d <- map[j+q,2] - map[j,2]
					error.terms <- rbind(error.terms,c(map[j,1],map[j+q,1],d,i))
				}	

				#wrap-around 
				for (j in (m-q+1):m) {
					n.terms[i] <- n.terms[i]+1
					d <- map[j,2] - map[(j+q)%%m,2]
					error.terms <- rbind(error.terms,c(map[(j+q)%%m,1],map[j,1],d,i))
				}
			}
		}

		n.error.terms <- nrow(error.terms)
			
		print("Generating consensus map.")
		N <- n.bin + n.error.terms
		G <- cBind(A,Matrix(0,nrow=nrow(A),ncol=n.error.terms,sparse=TRUE))
	
		b <- rep(0,nrow(G))
		for (i in 1:n.error.terms) {
		v <- rep(0,N)
		v[error.terms[i,1]] <- -1
		v[error.terms[i,2]] <- 1
		v[n.bin+i] <- 1
		b <- c(b,error.terms[i,3])
		G <- rBind(G,v)

		v <- rep(0,N)
		v[error.terms[i,1]] <- 1
		v[error.terms[i,2]] <- -1
		v[n.bin+i] <- 1
		b <- c(b,-error.terms[i,3])
		G <- rBind(G,v)
	}
	
		f <- c(rep(0,n.bin),weights[error.terms[,4]]/n.terms[error.terms[,4]])  
		ans <- Rglpk_solve_LP(f,G,rep(">=",nrow(G)),b)	
		if (ans$status != 0) {
		stop("Error in LP solver.")	
		} else {
		map <- data.frame(marker=markers,position=ans$solution[match(mark.bins,bins)],stringsAsFactors=F)
		composite.map <- map[order(map$position),]
		row.names(composite.map) <- NULL

        print(paste("Max.Interval = ",max.interval[p],sep=""))
		print(paste("Consensus map length:",max(composite.map$position)))

		link.maps <- numeric(0)
		for (k in 1:n.maps) {
			map <- Maps[[k]]
			ix <- match(composite.map$marker,map[,1])
			link.maps <- cbind(link.maps,ifelse(is.na(ix),NA,map[ix,2]))
		}
		RMSE <- apply(link.maps,2,function(x){sqrt(mean((composite.map$position-x)^2,na.rm=TRUE))})
		print(data.frame(map=c(map.names,"mean","sd"),RMSE=round(c(RMSE,mean(RMSE),sd(RMSE)),2)))
		colnames(link.maps) <- map.names
		result[[p]] <- data.frame(marker=composite.map$marker,consensus=composite.map$position,link.maps,stringsAsFactors=F)
		}
	}
	return(result)
}
