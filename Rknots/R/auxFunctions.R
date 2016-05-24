##Script generated in:
# 2011
# 8:00:08 PM
#by: 
# Authors: Federico Comoglio @ D-BSSE, ETH Zurich
# 		   and Maurizio Rinaldi @ University of Piemonte Orientale
###############################################################################

norM <- function(x) sqrt(sum(x ^ 2))

isWholeNumber <- function(x, tol = .Machine$double.eps ^ 0.5)  
	abs(x - round(x)) < tol

localize <- function(i, set) 
	which(sort(c(i, set)) == i)[1]

findComponent <- function(edge, i)  
	length(which(i < edge))

findNextE <- function(k, ends) {
	if  (k %in% ends == FALSE)
		return(k + 1)
	else 
		return(ends[which(ends == k) - 1][1] + 1)
}

rid2D <- function(points3D) 
	lapply(lapply(points3D, "[[", 1), "[", 1 : 2)

get2D <- function(points3D.l) {
	points3Dout <- points3D.l
	point2D <- rid2D( points3D.l )
	for (i in 1 : length(points3D.l)) 
		points3Dout[[i]][[1]] <- point2D[[i]]
	return(points3Dout)   
}

cyclicP <- function(multiins, ends)  {
	cyc <- function(ins) c(ins[2 : length(ins)], ins[1])
	k <- c()
	for (i in 1 : (length(ends) - 1)) 
		k <- c( k, cyc(multiins[(ends[i] + 1) : ends[i + 1]]) )
	return(k)
}

cp <- function(A) 
	sign(A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])

cross_product <- function(P1, P2) {
	A <- matrix(c(P1,P2), ncol = 2, byrow = TRUE)
	cp <- A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
	return(cp)
}

cpr <- function(A)  
	cp(A) * acos( A[1, ] %*% A[2, ] / (norM(A[1, ] * norM(A[2, ]))) )

#P = point, E1 = first point of the edge (segment), E2 = second one.
distancePS <- function(P, E1, E2) {
	h <- function(C, A, B) return( abs(cross_product(B - A, C - A) / norM(A - B)) )
	r <- function(C, A, B) return( sum( (C - A) * (B - A) ) / ( sum((A - B) ^ 2) ) )
	r.temp <- r(P, E1 , E2)
	if (r.temp > 1) 
		dist <- norM(P - E2)
	else
	if (r.temp < 0)  
		dist <- norM(P - E1)
	else
		dist <- h(P, E1, E2)
	return(dist)
}

distancepstart <- function(P, points3D, start) {
	points2D <- points3D[, 1:2]
	dmin <- Inf
	for(i in 1 : length(start))
		dmin <- min(dmin, distancePS(P, points2D[start[i], ], points2D[start[i] + 1, ]))
	return(dmin)
}

residualIndices <- function(indices, ends, npoints) {
	extends <- c(0, ends, npoints);
	dropstart <- sort( c(extends, indices, indices - 1) )
	(1 : npoints)[ -unique(dropstart) ]
}

orientedSign <- function (points3D, edge.indices) {
	x <- diff(points3D[edge.indices[1] + c(0, 1), ])[1:2]
	y <- diff(points3D[edge.indices[2] + c(0, 1), ])[1:2]
	A <- matrix(c(x, y), ncol = 2, byrow = TRUE)
	cp <- A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
	return(sign(cp))
}

skeinSign <- function (points3D, edge.indices) 
{	
	cp <- orientedSign( points3D, edge.indices )
	pointsij <- points3D[rep(edge.indices, each = 2) + c(0, 1), ]
	return(sign(cp) * singleIntersection(pointsij, "sign"))
}

lexOrder <- function(B) {
	nc <- max(nchar(B))
	temp <- apply(B, c(1,2), fillToN, nc)
	temp <- apply(temp, 1, paste, sep = "", collapse = "")
	return(order(temp))
}

fillToN <- function(string, n){
	nc <- nchar(string)
	if (nchar(string) < n)
		temp <- paste(
				paste(rep("0", n - nc), sep = "", collapse = ""),
				string, sep = "")
	else
		temp <- string
	return(temp)
}

parseToSympy <- function(polynomial) {
	parsed <- paste(unlist(strsplit(polynomial, split='\\^')), collapse = '**')
	return(parsed)
}

parseToR <- function(polynomial) {
	parsed <- paste( unlist( strsplit(polynomial, split = '[*][*]') ), collapse = '^' )
	return(parsed)
}






