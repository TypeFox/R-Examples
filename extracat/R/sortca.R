
tspsort = function(x, dims, perm.cat, vs = 1, comp = 2, iter = Inf, tsp.method = "cheapest_insertion", ...){
	# version sortca2
	if(! "ca" %in% .packages(all.available = TRUE) ){
		install.packages("ca")	
	}
	#require(ca)
	if (requireNamespace("ca", quietly = TRUE)) {
	dx <- as.data.frame(as.table(x))
	dx2 <- dx
	B <- Burt(dx)
	c1 <- ca::ca(B)
	nlvl <- dim(x)
	
	stopifnot(any(perm.cat))
	
	nd <- length(nlvl)
	#coords <- c1$rowcoord[,1:comp]*rep(c1$sv[1:comp],each=nrow(c1$rowcoord))
	coords <- c1$rowcoord %*% diag(sqrt(c1$sv))
	comp <- max(2,comp)
	if(comp == 2){
		iter <- 1	
	}
	crit <- 0
	
if(comp == 2){
	angs <- apply(coords[,1:comp],1,function(z){
		if(z[1] < 0){
			360-angle(c(0,1),z)
		}else{
			angle(c(0,1),z)
		}
	})
	
	grp <- rep(c(1:nd),nlvl)
	ord <- order(angs)
	angs <- angs[ord]
	grp <- grp[ord]
	ids <- unlist(lapply(nlvl,function(z) 1:z))[ord]
	
	diffs <- diff(c(angs[length(angs)]-360,angs))
	
	firstid <- which.max(diffs)
	
	angs <- angs - angs[firstid]
	angs[angs < 0] <- 360 + angs[angs < 0]
		
	orders <- tapply(angs,grp,order)
	ids <- tapply(ids,grp,I)
	
	for(i in 1:nd){
		dx[,i]	<- factor(dx[,i], levels=levels(dx[,i])[   ids[[i]] [ orders[[i]] ]  ])
		ids[[i]] <- ids[[i]] [ orders[[i]] ]
	}
	x <- xtabs(Freq~., data = dx)
	
		attr(x,"orders") <- ids
}else{
	#require(TSP)
		angs <- apply(coords[,1:comp],1,function(y){
			apply(coords[,1:comp],1,function(z){
				angle(y,z)	
			})	
		})
		(angs <- angs+t(angs))/2
		ts1 <- TSP(angs)
		m <- dim(B)[1]
		if(iter < m){
			starts <- sample(1:m,iter)
		}else{
			iter <- m
			starts <- 1:m
		}
	for(i in 1:iter){
		
		dx2 <- dx
		st1 <- solve_TSP(ts1, method = tsp.method,control = list(start = as.integer(i)))
		ord <- as.integer(st1)
		grp <- rep(c(1:nd),nlvl)
		ids <- unlist(lapply(nlvl,function(z) 1:z))
		orig.ids <- seq_along(ids)
	
		tmp<-cbind(ord,c(ord[length(ord)],ord[-length(ord)]))
		ord.angs <- apply(tmp,1,function(z){
			angs[z[1],z[2]]
		})
		firstid <- which.max(ord.angs)
		grp <- grp[ord]
		ids <- ids[ord]
		orig.ids <- orig.ids[ord]
		
		ord <- seq_along(ord)
		if(firstid > 1){
			ord[firstid:length(ord)] <- ord[firstid:length(ord)]-firstid+1
			ord[1:(firstid-1)] <- ord[1:(firstid-1)]+length(ord)-firstid+1
		}
		grp <- grp[ord]
		ids <- ids[ord]
		orig.ids <- orig.ids[ord]
		
		orders <- tapply(ids,grp,I)
		
    for(i in 1:nd){
			dx2[,i]	<- factor(dx2[,i], levels=levels(dx2[,i])[    orders[[i]]   ])
		}
		x2 <- xtabs(Freq~.,data=dx2)
		 crit2<-classcrit(x2)

		if( crit2 > crit){
			crit <- crit2
			x <- x2
			attr(x,"orders") <- orders
			attr(x,"tl") <- attr(st1,"tour_length")
			attr(x,"ids") <- orig.ids
			attr(x,"grp") <- grp
			attr(x,"ds") <- ord.angs[ord]
		}
		
	}# iter
	for(i in 1:nd){
			dx[,i]	<- factor(dx[,i], levels=levels(dx[,i])[    attr(x,"orders")[[i]]   ])
		}
}
		
	
attr(x,"criterion") <- BCI(x)

	#return(x)
	return(c( unlist(attr(x,"orders"))-1,attr(x,"criterion")))
	}else{
		stop("Please install package 'ca' to use this function.")
	}
}
