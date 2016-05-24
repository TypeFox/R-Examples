# Script comments and history
# 2011
# Feb 28, 2011
# 8:04:50 AM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

cross <- function(v1, v2) {
	cx <- v1[2] * v2[3] - v2[2] * v1[3]
	cy <- v1[1] * v2[3] - v1[3] * v2[1]
	cz <- v1[1] * v2[2] - v2[1] * v1[2]
	return(c(cx, -cy, cz))
}


interiorTriangle <- function(triangle2D, point) { #triangle2D is a 3x2 Matrix
	v1 <- triangle2D[2, ] - triangle2D[1, ] #the two edges
	v2 <- triangle2D[3, ] - triangle2D[1, ]
	#vectors2D <- rbind(triangle2D[2, ] - triangle2D[1, ], triangle2D[3, ] - triangle2D[1, ])
	#det.v <- vectors2D[1, 1] * vectors2D[2, 2] - vectors2D[1, 2] * vectors2D[2, 1]
	det.v <- v1[1] * v2[2] - v1[2] * v2[1]
	Ct <- matrix(c(v2[2], -v1[2], -v2[1], v1[1]), ncol = 2)
	xs <- (1 / det.v) * Ct %*% (point - triangle2D[1, ])
	inside <- (xs[1] > 0 & xs[2] > 0 & sum(xs) < 1)
	return(inside)
}


triangleIntersection <- function(triangle3D, segment) {
	v1 <- triangle3D[3, ] - triangle3D[1, ]
	v2 <- triangle3D[2, ] - triangle3D[1, ]
	n <- cross(v1, v2)	
	t <- (n %*% (segment[1, ] - triangle3D[1, ])) / (n %*% (segment[1, ] - segment[2, ]))
	P <- segment[1, ] + t * (segment[2, ] - segment[1, ])
	inside <- interiorTriangle(triangle3D[, 1 : 2], P[1 : 2])
	return(inside)
} 


segmentSet <- function(comp, i, extends)
{
	start <- 1 : extends[length(extends)]
	indextr <- extends[comp] + i + 0:3
	set <- setdiff(setdiff(start, indextr), extends + 1) - 1
	return(set)
}
	
move3D <- function(points3D, ends) {
	ncomp <- length(ends) + 1
	repeat { 
		extends <- c(0, ends, nrow(points3D))
		npoints <- nrow(points3D)
		npoints.in <- npoints
		component <- 1 
		while (component <= ncomp){
			perm <- c(3, 1, 2, 3)
			points3D.comp <- points3D[(extends[component] + 1) : extends[component + 1], ]
			npoints.comp <- nrow(points3D.comp)
			while (max(perm) <= npoints.comp) 
			{
				triangle <- points3D.comp[perm, ]
				intersections <- c()
				segments <- segmentSet(component, perm[2], extends)
				for (i in segments) { 
					segment <- rbind(points3D[i, ], points3D[i + 1, ])
					intersections <- c(intersections, triangleIntersection(triangle, segment)) 
				}
				if (length(intersect(intersections, TRUE)) != 0) {
					perm <- perm + 1
					next
				}
				#total <- unique( intersections )
				#if (identical(total, FALSE) | identical(total, c())) 
				if( !all(intersections) || is.null(intersections) )
					points3D <- points3D[-(extends[component]+perm[3]), ]
				perm <- perm + 1
				npoints <- nrow(points3D)
				npoints.comp <- npoints.comp - 1
				extends[(component + 1) : length(extends)] <- extends[(component + 1) : length(extends)] - 1
				ends <- extends[-c(1, length(extends))]
				points3D.comp <- points3D[(extends[component] + 1) : extends[component + 1],]
				npoints.comp <- nrow(points3D.comp)
			}
			component <- component + 1
		}
		if (npoints == npoints.in) 
			break
	}
	return(list(points3D = points3D, ends = extends[-c(1, length(extends))]))
}



AlexanderBriggs <- function(points3D, ends = c()) {
	points3D <- move3D(points3D, ends)
	return(points3D)
}


