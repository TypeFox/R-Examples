"tacvfARMA"   <-  
function(phi = numeric(0), theta = numeric(0), maxlag = 20, useCt = TRUE, sigma2 = 1)
{

	if(useCt) {
		dd <- function(x, y, ml) .C("tacvfARMA_C", as.double(x), as.integer(length(x)), as.double(y), as.integer(length(y)), as.integer(ml), res = double(ml+1))$res
		return(sigma2*dd(phi, theta, maxlag))
	}
	p   <-   length(phi)
	q   <-   length(theta)
	maxlagp1   <-   maxlag + 1

	if(max(p, q) == 0) {
			return(c(sigma2, numeric(maxlag)))
	}
	r   <-   max(p, q) + 1
	b   <-   numeric(r)
	C   <-   numeric(q + 1)

	C[1]   <-   1
	theta2   <-   c(-1, theta)
	phi2   <-   numeric(3 * r)
	phi2[r]   <-   -1
	if(p > 0) {
		phi2[r + 1:p]   <-   phi
	}
	if(q > 0) {
		for(k in 1:q) {
			C[k + 1]   <-    - theta[k]
			if(p > 0) {
				for(i in 1:min(p, k)) {
				  C[k + 1]   <-   C[k + 1] + phi[i] * C[k + 1 - i]
				}
			}
		}
	}

	for(k in 0:q) {
		for(i in k:q) {
			b[k + 1]   <-   b[k + 1] - theta2[i + 1] * C[i - k + 1]
		}
	}

	if(p == 0) {
		g   <-   c(b, numeric(maxlagp1))[1:maxlagp1]
		
		return(g)
	}
	else if(p > 0) {
		a   <-   matrix(numeric(r^2), ncol = r)
		for(i in 1:r) {
			for(j in 1:r) {
				if(j == 1) {
				  a[i, j]   <-   phi2[r + i - 1]
				}
				else if(j != 1) {
				  a[i, j]   <-   phi2[r + i - j] + phi2[r + i + j - 
				    2]
				}
			}
		}


		g   <-   solve(a,  - b)
		
		if(length(g) <= maxlag) {
			g   <-   c(g, numeric(maxlagp1 - r)) 

			for(i in (r + 1):maxlagp1) {
				g[i]   <-   phi %*% g[i - 1:p]
			}
		
			return(sigma2*g[1:maxlagp1])
		}
		else if(length(g) >= maxlagp1) {
			
			return(sigma2*g[1:maxlagp1])
		}
	}
}
