##Script generated in:
# 2011
# 11:41:38 AM
#by: 
# Author: Maurizio Rinaldi @ UPO, Novara, Italy
###############################################################################

listVertices <- function(points3D, ord.int)
{
	vertexlist <- ord.int
	for (i in 1 : length(ord.int)) {
		indicesij <- c(ord.int[[i]][[1]][1], ord.int[[i]][[1]][2])
		sign <- ord.int[[i]][[2]]
		oriented.sign <- orientedSign(points3D, indicesij) * sign
		vertexlist[[i]] <- list(ord.int[[i]][[4]], sign, oriented.sign)}
	return(vertexlist)
}


vertexPresentation <- function(points3D, ends = c()) {
	npoints <- nrow(points3D)
	nends <- length(ends)
	nedge <- npoints - 1
	M <- intersectionMatrix(points3D, ends)
	undercross <- which(M == -1, arr.ind = TRUE)
	if( nrow(undercross) < 2) #no need to go ahead. ##
		return( c() ) ##
	nunder <- nrow(undercross)
	radii <- c()
	Qp.list <- e.under <- e.over <- list()
	
	for (i in 1 : nunder) {
		Qp.list[[i]] <- singleIntersection(rbind(points3D[undercross[i, 1] + c(0, 1), ],
						points3D[undercross[i, 2] + c(0, 1),]), "ks")
		e.under[[i]] <- list(undercross[i, ], Qp.list[[i]][[1]],
				Qp.list[[i]][[2]][1], pt = Qp.list[[i]][[3]][1, ])
		Qs2D <- e.under[[i]][["pt"]]
		e.over[[i]] <- list(undercross[i, 2 : 1], -Qp.list[[i]][[1]],
				Qp.list[[i]][[2]][2], pt = Qp.list[[i]][[3]][2, ])
	}
	a <- unordered.edges <- c(e.under, e.over)
	n <- length(unordered.edges)
	B <- matrix(nrow = n, ncol = 2)
	
	for  (i in 1 : n) {
		B[i, ] <- a[[i]][[1]]
	}
	orderb <- lexOrder(B)
	ordered.edges <- a[orderb]
	B <- matrix(nrow = n, ncol = 2)
		
	for  (i in 1 : n) {
		B[i, ] <- ordered.edges[[i]][[1]]
	}
	par.list <- vector("list", nedge)
	for (i in 1 : nrow(B)) 
		par.list[[B[i, 1]]] <- c(par.list[[B[i, 1]]], ordered.edges[[i]][[3]])
	for (i in 1 : nedge) {
		if(length(par.list[[i]]) > 0) 
			par.list[[i]] <- order(par.list[[i]])
	}

	accumulate <- c(0, cumsum(lapply(par.list, length)))
	ordered.int <- list()
	for (i in 1 :nedge) 
		ordered.int[[i]] <-  par.list[[i]] + accumulate[i]
	ordered.int <- unlist(ordered.int)
	ordered.int <- ordered.edges[ordered.int]
	vertex.presentation <- listVertices(points3D, ordered.int)
	vertex.range <- sapply(ordered.int, "[", 1)
	points.range <- (1 : npoints) - 0.5
	temp <- as.vector(sapply(vertex.range, "[", 1))
	abs.position <- sort(c(points.range, temp)) 
	position.to.vertex <- isWholeNumber(abs.position)
	position.not.vertex <- (1 : length(abs.position))[isWholeNumber(abs.position) == FALSE]
	points3Dout <- abs.position
	points3Dout[which(position.to.vertex == TRUE)] <- vertex.presentation
	manipulated <- list()
	for (i in 1 : npoints)
		manipulated[[i]] <- list(points3D[i, ], 0, 0)
	if (nends > 0) {
		for (i in 1 : nends)  
			manipulated[[ends[i]]] <- c(manipulated[[ends[i]]], "end")
	}
	points3Dout[position.not.vertex] <- manipulated
	lp3 <- sapply(points3Dout, length)
	endsout <- which(lp3 == 4)
	zeroes <- which(lp3 == 0) 
	if (length(zeroes) > 0) 
		points3Dout <- points3Dout[-zeroes]
	doubles <- c()
	for (i in 2 : length(points3Dout)) {
		if (identical(points3Dout[[i]], points3Dout[[i - 1]])) 
			doubles <- c(doubles,i)
	}
	shift <- length(which(doubles < endsout[1]))
#if(length(unique(points3Dout)) != length(points3Dout))
#		warning("delete(points) not equal to points.")
	endsout <- endsout - shift
	return(list(points3Dout = points3Dout,
					endsout = endsout,
					M = M))
}


auxiliaryAlexander <- function(points2D, ends) {
	undercross <- which(sapply(points2D, "[", 2) == -1)
	comp.ends <- 0
	for (i in 1 : length(ends)) {
		if(length(ends > 0)) 
			comp.ends <- c(comp.ends, length(which(undercross < ends[i]) == TRUE))
	}
 #unique inserted
	comp.ends <- unique(c(comp.ends, length(undercross)))
	nunder <- length(undercross)
	pos.Qs.over <- which(sapply(points2D, "[" , 2) == 1)
	nover <- length(pos.Qs.over)
	p22 <- sapply(points2D[pos.Qs.over], "[", 1)
	p2 <- sapply(points2D[undercross], "[", 1)
	pos.Q.over <- c()
	for (i in 1 : nover)
		pos.Q.over[i] <- pos.Qs.over[which(lapply(p22, identical, p2[[i]]) == TRUE)]
	indices <- list()
	for (i in 1 : length(pos.Q.over))
		indices[[i]] <- c(i, endComponentCorrection(pos.Q.over[i], ends, undercross, length(points2D)))
	outgoing.under <- sapply(sapply(indices, "[", 1), findNextE, comp.ends)
	triples <- list()
	for (i in 1 : length(outgoing.under))
		triples[[i]] <- unlist(list(outgoing.under[i], indices[[i]]))
	temp <- order(outgoing.under)
	signs <- unlist(sapply(points2D[which(sapply(points2D, "[", 2) == -1)], "[" , 3))[temp]
	triples <- triples[temp]
	return(list(nunder = nunder, 
				triples = triples, 
				comp.ends = comp.ends, 
				signs = signs,
				pos.Q.over = pos.Q.over,
				ends = ends,
				undercross = undercross))
}


endComponentCorrection <- function(n, ends, undercross, npoints) {
	extends <- c(0, ends, npoints) 
	temp <- extends[which(extends >= n)[1]]
	join <- c(temp,
			if(n > tail(undercross, 1))
				extends[length(extends)] + 1
			else
				undercross[localize(n, undercross)][1]) 
	if(order(join)[1] == 2)
		return(localize(n, undercross))
	else
		return(localize(extends[(which(extends == temp[1]) - 1)][1], undercross))
}

singleIntersectionS <- function (pointsij, kind = "binary") 
{
    sint <- 0
    pointsij2D <- pointsij[, c(1, 2)]
    vectors3D <- matrix(c(pointsij[2, ] - pointsij[1, ], pointsij[4, 
        ] - pointsij[3, ]), nrow = 2, byrow = TRUE)
    vectors2D <- vectors3D[, c(1, 2)]
    if (nrow(unique(pointsij)) <= 3) 
        return(0)
    detv <- vectors2D[1, 1] * vectors2D[2, 2] - vectors2D[1, 
        2] * vectors2D[2, 1]
    if (identical(detv, 0)) {
        if (pointsij[1, ] == pointsij[2, ] || pointsij[3, ] == 
            pointsij[4, ]) {
            warning("collapsed edge")
            return(0)
        }
        else {
            E <- t(rbind(vectors2D, pointsij2D[3, ] - pointsij2D[1, 
                ]))
            E[2, ] <- E[2, ] - (E[2, 1]/E[1, 1]) * E[1, ]
            if (identical(E[2, ], c(0, 0, 0))) {
                ksi <- c(0, 0)
                segment <- pointsij2D[2, ] - pointsij2D[1, ]
                for (i in 1:2) ksi[i] <- segment %*% (pointsij2D[i + 
                  2, ] - pointsij2D[1, ])/(segment %*% segment)
                if ((0 < ksi[1] && ksi[1] < 1) || (0 < ksi[2] && 
                  ksi[2] < 1)) {
                  warning("superimposed edges")
                  return(0)
                }
            }
        }
    }
    else {
        inverse <- (1 / detv) * matrix(c(vectors2D[2, 2], vectors2D[1, 
            2], -vectors2D[2, 1], -vectors2D[1, 1]), nrow = 2)
        temp <- pointsij2D[3, ] - pointsij2D[1, ]
        ks <- inverse %*% temp
        conditionI <- ((0 < ks[1] && ks[1] < 1) && (0 < ks[2] && 
            ks[2] < 1))
        if (!conditionI) 
            return(0)
		else
            return(c(1, detv))
    }
    return(sint)
}


#this function inputs a (closed) polygonal link and returns the (counterclockwise) rotation factor 
#of the components
  
rotFactor <- function(points3D, ends) {
	points2D <- points3D[, 1 : 2]
	lp <- c(0, ends, nrow(points2D)) 
	rot <- list()
	for (j in seq(length(lp) - 1)) {
			rot[[j]] <- list()
    		for (i in (lp[j] + 1) : (lp[j + 1] - 1)) {
				s1 <- i + 1 - lp[j + 1]
       			if (s1 > 0) 
					r1 <- lp[j] + s1 + 1 
				else 
					r1 <- i + 1
 	 			a1 <- points2D[r1, ] - points2D[i, ]
  				s2 <- i + 2 - lp[j + 1] 
  				if(s2 > 0) 
					r2 <- lp[j] + s2 + 1 
				else 
					r2 <- i + 2
				a2 <- points2D[r2, ] - points2D[r1, ]
 				rot[[j]][[i - lp[j]]] <- rbind(a1, a2)
    		}
	}
rotf <- c()
for (i in seq(length(rot))) 
	rotf[i] <- round(sum(sapply(rot[[i]], cpr)) / (2 * pi))
return(rotf)
}


# this function determines the direction of the escape path
escapeDir <- function(x1, x2)   
	cp(rbind(x1, x2)) * (x1 - x2)


buildRis <- function(points3D, ends, compute.points2D = FALSE) {
  	rv <- vertexPresentation(points3D, ends)
	if(is.null( rv )) #meaning there is no cross.
		return( c() ) #no need to go ahead.
  	endsout <- rv[[2]]
  	points3Dout <- rv$points3Dout
	if(compute.points2D)
  		points2D <- t(sapply(sapply(points3Dout, "[" , 1) , "[" , c(1, 2)))
	point2D <- rid2D(rv$points3Dout)
	for (i  in 1 : length(points3Dout)) 
		points3Dout[[i]][[1]] <- point2D[[i]]
	ris <- auxiliaryAlexander(points3Dout, endsout)
	if(compute.points2D)
		return(list(ris, points2D))
	else
		return(ris)
}


escapePath <- function(points2D, undercross, positionQsovercross, endout) {
	nunder <- length(undercross)
	x1 <- undercross[nunder - 1]
	x2 <- positionQsovercross[nunder - 1]
	a1 <- points2D[x2 + 1, ] - points2D[x2, ] 
	a2 <- points2D[x1 + 1, ] - points2D[x1, ]
	np <- nrow(points2D)
	diam <- sqrt(max(apply((points2D - matrix(apply(points2D, 2, mean), 
		nrow = np, ncol = 2, byrow = TRUE)) ^ 2, 1, sum))) * 2
	epsilon <- .5 *distancepstart( points2D[x1, ], points2D, 
			residualIndices(c(x1, x2), endout, nrow(points2D)) )
	p <- points2D[x2, ] + epsilon * escapeDir(a1, a2)
	p.out <- p + diam *(p - points2D[x2, ]) / norM(p - points2D[x2, ])
	set.int <- list()
	for(i in 1 : (nrow(points2D) - 1))  
		set.int[[i]] <- singleIntersectionS(rbind(p, p.out, points2D[i : (i + 1), ]), "binary")  
	temp <- setdiff(which(sapply(set.int, length) > 1), endout)
	word <- c()
	sign <- c()
	exponent <- rep(0, 1 + length(endout))
	if (length(temp) > 0) {
		for (i in 1 : length(temp)) {
			word[i] <- findComponent(temp[[i]], endout) + 1
			sign[i] <- sign(set.int[[temp[i]]][2])
		}
		exponent <- c()
		pf <- cbind(word, sign)
		for(i in 1 : (1 + length(endout)))
			exponent[i] <- sum(pf[which(word == i), 2])
	}
	return(list(exponent = exponent,
			p = p,
			p.out = p.out,
			set.int = set.int))
}

computeFactors <- function (points3D, ends) {
	ris <- buildRis(points3D, ends, TRUE)
	temp <- ris[[1]]
	points2D <- ris[[2]]
	endout <- temp[[6]]
	overf <- c()
	for (i in seq(length = length(temp[[5]]))) 
		overf[i] <- findComponent(temp[[5]][i], endout)
	overf <- as.vector(table(overf + 1))
	undercross <- temp[[7]]
	positionQsovercross <- temp[[5]]
	exponent <- escapePath(points2D, undercross, positionQsovercross, endout)	
	rotf <- rotFactor(points3D, ends)
	toeval7 <- paste("endout=[", paste(endout, collapse = ","), "]", sep = "")
	toeval8 <- paste("esponente=[", paste(exponent$exponent, collapse = ","), "]", sep = "") #edited
	toeval9 <- paste("rotf=[", paste(rotf, collapse = ","), "]", sep = "")
	toeval10 <- paste("overf=[", paste(as.vector(overf), collapse = ","), "]", sep = "")
	sympy(toeval7)
	sympy(toeval8)
	sympy(toeval9)
	sympy(toeval10)
	sympy("t=zeros((1,1+len(endout)));
i=-1
while (i<len(endout)):
	i=i+1;
	t[i]=Symbol('t'+str(i+1))
z=zeros((1,len(endout)+1))
word=1
i=-1
while(i<len(endout)):
	i=i+1 
	z[i]=t[i]**(esponente[i])
	word=word/z[i]*t[i]**((rotf[i]/2.-overf[i]/2.))"
)
	#polnorm <- sympy("factor(archdet)*word/(t[len(endout)]-1)")
	word <- sympy("word")
	return(list(word, rotf, overf, exponent))
}


mVA <- function (points3D, ends = c(), normalized = FALSE, return.A = FALSE)
{
	ris <- buildRis(points3D, ends, FALSE)
	if(is.null( ris )) {								##unknot/unlink for sure
		return(ifelse( ( is.null(ends) || identical(ends, numeric(0))), '1', '0') )		##
	}												##
	endout <- ris[[6]]
	FC <- sapply(1:ris[[1]], findComponent, ris[[3]])
	triples <- matrix(0, nrow = ris[[1]], ncol = 3)
	for (k in 1:nrow(triples)) triples[k, ] <- ris[[2]][[k]] - 1
	row.str <- apply(triples, 1, paste, collapse = ",")
	matrix.str <- paste("[[", paste(row.str, collapse = "] , ["), "]]", sep = "")
	toeval1 <- paste("nunder=", as.character(ris[[1]], sep = ""))
	toeval2 <- paste("triple=Matrix(", matrix.str, ")", sep = "")
	toeval3 <- paste("eends=[", paste(ris[[3]], collapse = ","), "]", sep = "")
	toeval4 <- paste("signs=[", paste(ris[[4]], collapse = ","), "]", sep = "")
	toeval5 <- paste("FC=[", paste(FC, collapse = ","), "]", sep = "")
	toeval6 <- paste("endout=[", paste(endout, collapse = ","), "]", sep = "")
	sympy(toeval1)
	sympy(toeval2)
	sympy(toeval3)
	sympy(toeval4)
	sympy(toeval5)
	sympy(toeval6)
	sympy("t=zeros((1,len(eends)))
i=-1
while (i<len(eends)-1):
	i=i+1
	t[i]=Symbol('t'+str(i+1))
Arch=zeros(nunder)
k=0
while(k<nunder):
	if (signs[k]==+1):
		Arch[triple[k,0],triple[k,1]]=Arch[triple[k,0],triple[k,1]]-1
		Arch[triple[k,0],triple[k,0]]=Arch[triple[k,0],triple[k,0]]+t[-1+FC[int(triple[k,2])]]
		Arch[triple[k,0],triple[k,2]]=Arch[triple[k,0],triple[k,2]]+1-t[-1+FC[int(triple[k,1])]]
	else:
		Arch[triple[k,0],triple[k,1]]=Arch[triple[k,0],triple[k,1]]-t[-1+FC[int(triple[k,2])]]
		Arch[triple[k,0],triple[k,0]]=Arch[triple[k,0],triple[k,0]]+1
		Arch[triple[k,0],triple[k,2]]=Arch[triple[k,0],triple[k,2]]+t[-1+FC[int(triple[k,1])]]-1
	k=k+1
archdet=Arch[0:(nunder-1),0:(nunder-1)].det()
p = Poly(archdet, t[0])
q=zeros((1,p.degree+1))
for i in range(0,p.degree+1):
	q[i]=p.coeff(i)
a=0
i=-1
while (a==0 and p.degree>0):
	i=i+1
	a=q[i]
r=p.degree
polyB=factor(archdet)/(t[len(endout)]-1)"
)

if(( is.null(ends) || identical(ends, numeric(0)) )) { #if it is a knot
	sympy( "poly=(factor(archdet)/t[0]**(int((r+i)/2))).expand()" ) ##factor added
	polynomial <- sympy("poly")
	}
else { #it is a link
	if(!normalized) { #and normalized = FALSE
		polynomial <- sympy('polyB')
	}
		else { #and normalized = TRUE
				word <- computeFactors(points3D, ends)[[1]]
				polynomial <- sympy("polyB * word")
			}
}
polynomial <- parseToR(polynomial)

	if (normalized & !return.A) return(polynomial)
	if (normalized & return.A) {
		A <- sympy('Arch')
		return(list( polynomial = polynomial, A = A) )
	}
	if (!normalized & !return.A) return(polynomial)
	if (!normalized & return.A) {
		A <- sympy('Arch')
		return(list( polynomial = polynomial, A = A) )
	}
}