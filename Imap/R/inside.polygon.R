inside.polygon <- function(pts, h)
{

# Revised from the the function "Inside()" by Joseph S. Verducci (Snews: 09 Feb 1999)
# See <http://www.biostat.wustl.edu/archives/html/s-news/2002-07/msg00020.html>

# Input:
#    pts -- (n x 2) matrix of test points
#      h -- A hull or polygon defined by a matrix ([k+1] x 2) of (ordered) vertices, 
#                 [with last row = first row]
# Output:
#     (n x 1) logical vector indicating inclusion of tests points
#           in the convex polygon formed by h 

	n <- nrow(pts)
	k <- nrow(h) - 1
	if(!all(h[1,  ] == h[k + 1,  ])) {
		h <- rbind(h, h[1,  ])
		k <- k + 1
	}

# compute slopes and intercepts of boundary lines

        if(any(diff(h[,1]) == 0)) {
           E  <- (2:nrow(h))[diff(h[,1]) == 0]
           h[E, 1] <- h[E, 1] + E*1e-10 
        }
	b <- diff(h[, 2])/diff(h[, 1])
	a <- diag(as.matrix(h[1:k,  ]) %*% as.matrix(t(cbind( - b, 1))))
	
#   compute centroid and its residuals from boundary lines

	m <- t(as.matrix(rep(1/k, k))) %*% as.matrix(h[1:k,  ])
	r0 <- m[2] - (a + b * m[1])

#  compute residuals of test points from boundary lines

	r <- outer(pts[, 2], rep(1, k)) - outer(pts[, 1], b) - outer(rep(1, n), a)
	signs <- (r %*% diag(r0)) > 0
	ans <- as.vector(signs %*% rep(1, k) == k)
	ans
}

