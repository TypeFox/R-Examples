##Script generated in:
# 2011
# 11:10:49 AM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################


bin2N <- function(V.binary)
{
	n <- length(V.binary)
	for(i in 1 : n) 
		V.binary[i] <- V.binary[i] * 2 ^ (n-i)
	return(sum(V.binary))
}

buildLsw <- function(points3D, ends, rquad, Xsw, e)
{
	n <- nrow(points3D)
	nends <- length(ends)	
	Lsw <- rbind(points3D[1 : e, ],
			rquad[1, ],
			as.numeric(Xsw),
			rquad[2, ],
			points3D[(e + 1) : n, ])
	ends[ends > e] <- ends[ends > e] + 3
	return(list(Lsw = Lsw, ends = ends))
}


buildLzero <- function(points3D, ends, rquad, e1, e2)
{
	findComplement <- function(set, subset) {
		#if(length(set) < length(subset))
		#	warning("set longer than subset\n")
		complement <- set[!is.element(set, subset)]
		return(complement)
	}
	
	n <- nrow(points3D)
	nends <- length(ends)
	temp <- sort(c(ends, e1, e2))
	icomp <- which(temp == e1)
	jcomp <- which(temp == e2) - 1
	samecomp <- abs(icomp - jcomp) == 0
	
	if(samecomp) {	
		Lzero <- rbind(points3D[1 : e1, ],
				rquad[c(1,4), ],
				points3D[(e2 + 1) : ends[icomp], ],
				rquad[2, ], 
				points3D[(e1 + 1) : e2, ],
				rquad[c(3,2), ],
				if(icomp < nends) points3D[(ends[icomp] + 1) : n, ]
		)
		ends <- c(if((icomp - 1) >= 1)	ends[1 : (icomp - 1)],
				ends[icomp] - e2 + e1 + 2,
				if(nends >= 1) ends[icomp : nends] + 5)
	}
	else {
		if(ends[jcomp - 1] + 1 == ends[jcomp]) {
			points3D <- points3D[-ends[jcomp]]
			ends[jcomp : nends] <- ends[jcomp : nends] - 1
		}
		Lzero <- rbind(points3D[1 : e1, ],
				rquad[c(1,4), ],
				points3D[(e2 + 1) : ends[jcomp], ],
				points3D[(ends[jcomp - 1] + 1) : e2, ],
				rquad[c(3,2), ],
				points3D[(e1 + 1) : ends[icomp], ],
				points3D[findComplement((ends[icomp] + 1) : ends[nends],
								(ends[jcomp - 1] + 1) : ends[jcomp]), ])
		shift <- ends[jcomp] - ends[jcomp - 1] + 4
		ends <- c(if((icomp - 1) >= 1)	ends[1 : (icomp - 1)],
				ends[icomp : (jcomp - 1)] + shift,
				if((jcomp + 1) <= nends) ends[(jcomp + 1) : nends] + 4)
	}
	return(list(Lzero = Lzero, ends = ends))
}



checkVertex <- function(points3D, ends, ij, quadrilateral, k)
{
	vertices <- c(k, (k %% 4) + 1)
	intersections <- c()
	for(i in 1 : (nrow(points3D) - 1))
	{
		pointsij <- rbind(quadrilateral[vertices, ], points3D[i : (i + 1), ])
		intersections[i] <- singleIntersection(pointsij, "binary")
	}
	todel<-c(ij[, 1], ends)
	intersections <- intersections[-todel]
	evaluation <- sum(as.numeric(unique(intersections) != 0))
	return(evaluation)
}


xClean <- function(points3D, ends, quadrilateral, ij, S, toevaluate)
{
	
	torusShift <- function(n, shift) {
		temp <- (n + shift) %% 4
		if(temp == 0) temp <- 4
		return(temp)
	}
	
	for(i in 1 : length(toevaluate))
	{
		in.out <- bin2N(S[toevaluate[i],])
		switch(in.out,
				"1" = {k <- toevaluate[i];
					pos <- 2;
					neigh <- matrix(c(torusShift(k[1], 1), 1), nrow = 1)},
				"2" = {k <- torusShift(toevaluate[i], -1)
					pos<-1;
					neigh <- matrix(c(k[1],2), nrow = 1)},
				"3" = {k <- c(torusShift(toevaluate[i], -1), toevaluate[i]);
					pos <- c(1,2);
					neigh <- matrix(c(k[1], 2, torusShift(k[2], 1), 1), nrow = 2, byrow = TRUE)}
		)
		for(j in 1 : length(k))
		{
			S[toevaluate[i], pos[j]] <- checkVertex(points3D, ends, ij, quadrilateral, k[j])
			S[neigh[j, 1], neigh[j, 2]] <- S[toevaluate[i], pos[j]]
		}
	}
	return(S)
}


buildQuadrilateral <- function(points3D, ends, pointsij, ij, Xs, k)
{
	normV <- function (v)	return(sqrt(sum(v ^ 2)))
	
	ks <- rep(k, 4)
	S <- matrix(1, nrow = 4, ncol = 2)
	perm <- c(1, 3, 2, 4)
	distances <- rep(Inf, 4)
	vertices <- pointsij[perm, ]
	
	for(i in 1 : nrow(vertices))
	{
		vertices[i, ] <- Xs[i, ] + ks[i] * (vertices[i, ] - Xs[i, ])
		distances[i] <- normV(c(vertices[i, 1:2], 0) - c(Xs[i, 1:2], 0))
	}
	S <- xClean(points3D, ends, vertices, ij, S, 1:4)
	if(sum(S) == 0) 
		return(vertices[perm, ])
	
	repeat
	{
		ones <- unique(which(S == 1, arr.ind = TRUE)[, 1])
		tohalf <- which(distances == max(distances[ones]))      
		ks[tohalf] <- ks[tohalf] * 0.5
		vertices[tohalf, ] <- Xs[tohalf, ] + 
				ks[tohalf] * (pointsij[perm, ][tohalf, ] - Xs[tohalf, ])
		distances[tohalf] <- normV(c(vertices[tohalf, 1:2], 0) - c(Xs[tohalf, 1:2], 0))
		S[tohalf, ] <- c(1, 1)
		S <- xClean(points3D, ends, vertices, ij, S, tohalf)
		if(sum(S) == 0) 
			break
	}
	return(vertices[perm, ])
}


rotateQuadrilateral <- function(points3D, ends, quadrilateral, ij, Xs, ialpha)
{
	perm <- c(1, 3, 2, 4)
	vertices <- quadrilateral[perm, ]
	rquad <- vertices
	alpha <- ialpha;
	repeat
	{
		R <- t(rMatrix(alpha))
		for(i in 1 : nrow(vertices)) 
			rquad[i, ] <- Xs[i, ] + (R %*% (vertices[i, ] - Xs[i, ]))
		S <- matrix(1, nrow = 4, ncol = 2)
		S <- xClean(points3D, ends, rquad, ij, S, 1 : 4)
		if(sum(S) == 0) 
			break
		else 
			alpha <- alpha * 0.5
	}
	return(rquad[perm, ])
}


buildTriangles <- function(points3D, ends, vertices, quad, rquad, ij)
{
	intersections <- rep(list(rep(list(list(list(sign = c(), ks = c()))), 3)), 4)	
	extendedij <- sort(rep(ij, 2))
	nedge <- nrow(points3D) - 1
	triangles <- array(c(vertices, quad, rquad), dim = c(4, 3, 3))
	for(i in 1 : 4)
	{
		for(j in 1 : 3)
		{
			for(k in 1 : nedge)
			{
				pointsij <- rbind(triangles[i, , j], 
						triangles[i, , (j %% 3) + 1],
						points3D[k : (k + 1), ])
				temp <- singleIntersection(pointsij, "k")
				
				if(length(temp) > 1) {
					intersections[[i]][[j]][[1]]$sign <- 
							c(intersections[[i]][[j]][[1]]$sign, temp$sign)	
					intersections[[i]][[j]][[1]]$ks <- 
							c(intersections[[i]][[j]][[1]]$ks, temp$ks)
				}
				else { 
					intersections[[i]][[j]][[1]]$sign <- 
							c(intersections[[i]][[j]][[1]]$sign, 0)	
					intersections[[i]][[j]][[1]]$ks <- 
							c(intersections[[i]][[j]][[1]]$ks, 0)
				}
			}
			intersections[[i]][[j]][[1]]$sign <- 
					intersections[[i]][[j]][[1]]$sign[-c(extendedij[i], ends)]
			intersections[[i]][[j]][[1]]$ks <- 
					intersections[[i]][[j]][[1]]$ks[-c(extendedij[i], ends)]	
		}	
	}
	return(intersections)
}


checkTriangles <- function(intersections)
{
	check <- rep(1, 4)
	for(i in 1 : 4)
	{
		triangle <- intersections[[i]]	
		toedit <- which(triangle[[3]][[1]]$sign != 0)
		triangle[[3]][[1]]$ks[toedit] <- 1 - triangle[[3]][[1]]$ks[toedit]
		intposition <- lapply(triangle, function(x) which(x[[1]]$sign != 0))
		
		if(length(intposition[[2]]) == 0 && identical(intposition[[1]], intposition[[3]])) {
			signcheck <- identical(triangle[[3]][[1]]$sign[toedit],
					triangle[[1]][[1]]$sign[toedit])
			ordercheck <- identical(order(triangle[[3]][[1]]$ks[toedit]),
					order(triangle[[1]][[1]]$ks[toedit]))
			
			if(signcheck && ordercheck) 
				check[i] <- 0	
			else
				break
		}
	}
	return(check)
}


auxiliaryPoints <- function(points3D, ends, ancestor.signs, M) {
	uM <- M
	uM[lower.tri(uM)] <- 0
	#uM[ends, ] <- uM[ends, ] * 0
	#uM[, ends] <- uM[, ends] * 0	
	moves <- which(uM == -1, arr.ind = TRUE)
	temp <- which(uM == 1, arr.ind = TRUE)
	n.under <- nrow(moves)
	n.over <- nrow(temp)
	
	if(n.under == 0 || n.over == 0) {
		warning("No undercross or overcross left")
		return(c())
	}
	
	npoints <- nrow(points3D)
	track.length.Lsw <- rep(Inf, n.under)
	cached.Lsw <- vector("list", n.under)
	
	for(c in 1 : n.under) {
		es <- as.numeric(moves[c, ])
		e1 <- es[1]
		e2 <- es[2]
		ij <- matrix(rep(es, each = 2) + c(0,1), ncol = 2, byrow = TRUE) 
		pointsij <- points3D[t(ij), ] 
		skeinsign <- skeinSign(points3D, es)
		skeinsigns <- c(ancestor.signs, skeinsign)
		temp<- matrix(c(pointsij[2, ] - pointsij[1, ], pointsij[4, 
						] - pointsij[3, ]), nrow = 2, byrow = TRUE)
		temp[, 3] <- temp[, 3] * 0
		v <- temp[1, ]
		w <- temp[2, ]
		Xij <- singleIntersection(pointsij, "lk")$Qs
		Xs <- rbind(Xij, Xij)
		quad <- buildQuadrilateral(points3D, ends, pointsij, ij, Xs, k = 0.8)
		alpha.temp <- acos(sum(v * w) / (sqrt(sum(v^2)) * sqrt(sum(w^2))))
				
		alpha <- min(0.99 * alpha.temp, 0.99 * (pi - alpha.temp), pi / 8)
		alpha <- alpha * skeinsign * -1
		
		repeat {	
			rquad <- rotateQuadrilateral(points3D, ends, quad, ij, Xs, alpha)
			triangles<- buildTriangles(points3D, ends, pointsij, quad, rquad, ij[, 1])
			checked <- checkTriangles(triangles)
			if(sum(checked) != 0) 
				alpha <- alpha * 0.5
			else
				break
		}	
		k.sw <- 0.9
		k.z <- 1
		repeat {
			k.z <- k.z * 2
			Xsw <- Xij[1, ] + k.z * (Xij[2, ] - Xij[1, ])
			Xsw <- rquad[1, ] + k.sw * (Xsw - rquad[1, ])
			pointsij.test <- rbind(rquad[2, ], Xsw, rquad[c(4,3), ])

			if(singleIntersection(pointsij.test, "k")$sign == -1 && k.z < 2^50) 
				k.sw <- k.sw + ((1 - k.sw) / 2)
			else
				break
		}
		Lsw <- buildLsw(points3D, ends, rquad, Xsw, e1)
		track.length.Lsw[c] <- nrow(msrFast(Lsw$Lsw, Lsw$ends)$points3D)
		cached.Lsw[[c]] <- list(rquad = rquad, Xsw = as.numeric(Xsw), e1 = e1, e2 = e2, skeinsigns = skeinsigns)
	}	
	selected <- which.min(track.length.Lsw)
	return(cached.Lsw[[selected]])
}

skeinIterator <- function(points3D, ends, M = c()) {
	
	nextSet <- function(bin.list) {
		bin.list <- rep(bin.list, each = 2)
		n <- length(bin.list)
		tojoin <- rep(c(0, 1), (n / 2))
		for(i in 1 : n) bin.list[[i]] <- c(bin.list[[i]], tojoin[i])
		return(bin.list)
	}
	
	if(is.null(M)) 
		M <- intersectionMatrix(points3D)
	
	tree <- list(0)
	leaves <- list() 
	innerv <- list(1)
	tree[[innerv[[1]]]] <- list(points = points3D,
			ends = ends,
			signs = c(),
			M = M)
	
	ncomp <- length(ends) + 1
	if(nrow(points3D) == (2 * ncomp)) {
		return(list(leaves = leaves, tree = tree))
	}
	
	while(length(innerv) != 0) {
		nalive <- length(innerv)
		descendance <- nextSet(innerv)
		for(k in 1 : nalive) {
			index <- bin2N(innerv[[k]])
			dpoints <- tree[[index]]$points 
			n <- nrow(dpoints)
			dextends <- c(tree[[index]]$ends, n)
			dsigns <- tree[[index]]$signs
			dM <- tree[[index]]$M 
			param <- auxiliaryPoints(dpoints, head(dextends, -1) , dsigns, dM)
			if(is.null(param)) {
				warning("no -1 left")
				return(list(leaves = leaves, tree = tree))
			}
			newd <- bin2N(c(innerv[[k]], 0)) ##il primo e' un Lsw
			Lsw <- buildLsw(dpoints, dextends, param$rquad, param$Xsw, param$e1)
			temp <- msrFast(Lsw$Lsw, head(Lsw$ends, -1))
			tree[[newd]] <- list(points = temp$points3D,
					ends = temp$ends,
					signs = param$skeinsigns,
					M = temp$M)
			newd <- bin2N(c(innerv[[k]], 1)) ##il secondo e' un L0
			Lzero <- buildLzero(dpoints, dextends, param$rquad, param$e1, param$e2)			
			temp <- msrFast(Lzero$Lzero, head(Lzero$ends, -1))
			tree[[newd]] <- list(points = temp$points3D,
					ends = temp$ends,
					signs = param$skeinsigns,
					M = temp$M)
			
		}
		
		dMlist <- vector("list", length(descendance))
		toleaves <- rep(Inf, length(descendance))
		
		for(i in 1 : length(descendance)) {
			index <- bin2N(descendance[[i]])	
			tempM <- tree[[index]]$M
			tempends <- tree[[index]]$ends
			tempM[lower.tri(tempM)] <- 0
			dMlist[[i]] <- tempM
			toleaves[i] <- ((length(which(dMlist[[i]] == -1)) == 0) || (length(which(dMlist[[i]] == 1)) == 0))
		}	
		up <- which(toleaves == FALSE)
		down <- which(toleaves == TRUE)
		innerv <- lapply(up, function(x) descendance[[x]])
		leaves <- c(leaves, lapply(down, function(x) descendance[[x]]))			
	}
	return(list(leaves = leaves, tree = tree))
}
