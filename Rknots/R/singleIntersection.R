# Script comments and history
# 2011
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

singleIntersection <-
function (pointsij, kind = "binary") 
{
	sint <- 0
	pointsij2D <- pointsij[, c(1, 2)]
	vectors3D <- matrix(c(pointsij[2, ] - pointsij[1, ], pointsij[4, 
					] - pointsij[3, ]), nrow = 2, byrow = TRUE)
	vectors2D <- vectors3D[, c(1, 2)]
	if (nrow(unique(pointsij)) <= 3) 
		return(0)
	
	detv <- vectors2D[1, 1] * vectors2D[2, 2] -
			vectors2D[1,2] * vectors2D[2,1]
	if (identical(detv, 0)) {
		if (pointsij[1, ] == pointsij[2, ] || pointsij[3, ] == pointsij[4, ]) {
			warning("collapsed edge")
			return(0)
		}
		else {
			E <- t(rbind(vectors2D, pointsij2D[3, ] - pointsij2D[1, ]))
			E[2, ] <- E[2, ] - (E[2, 1] / E[1, 1]) * E[1, ]
			if(identical(E[2, ], c(0, 0, 0)))
			{
				ksi <- c(0, 0)
				segment <- pointsij2D[2, ] - pointsij2D[1, ]
				for (i in 1 : 2) 
					ksi[i] <- segment %*% (pointsij2D[i + 2, ] - pointsij2D[1, ]) / (segment %*% segment)
				if ((0 < ksi[1] && ksi[1] < 1) || (0 < ksi[2] && ksi[2] < 1)) {
					warning("superimposed edges")
					return(0)
				}
			}
		}
	}
	else {
		inverse <- (1 / detv) * matrix(c(vectors2D[2, 2], vectors2D[1, 2], -vectors2D[2, 1],
						-vectors2D[1, 1]), nrow = 2) 
		temp <- pointsij2D[3, ] - pointsij2D[1, ]
		ks <- inverse %*% temp
		conditionI <- ((0 < ks[1] && ks[1] < 1) && (0 < ks[2] && ks[2] < 1))
		if (!conditionI) 
			return(0)
		binary <- (kind == "binary")
		if (conditionI && binary)
			return(1)
		else {
			zs <- c(pointsij[1, 3] + ks[1] * vectors3D[1, 3], 
					pointsij[3, 3] + ks[2] * vectors3D[2, 3])
			sign <- sign(zs[1] - zs[2])
			Qs <- vectors3D
			ifelse (sign == 1, {i <- 1; j <- 2}, {i <- 2; j <- 1})
			Qs[i, ] <- pointsij[2 * i - 1, ] + ks[i] * vectors3D[i, ]
			Qs[j, ] <- Qs[i, ]
			Qs[j, 3] <- zs[j]
			sint <- switch(kind, 
					sign = sign, 
					q = list(sign = sign, Qs = Qs), 
					k = list(sign = sign, ks = ks[1]), 
					ks = list(sign = sign, ks = ks, Qs = Qs), 
					lk = list(sign = sign, ks = ks[1], Q2D = Qs[1, 1:2], Qs = Qs), 
			warning("unmatched call option"))
		}
	}
	return(sint)
}
