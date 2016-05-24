BUUHWT_2D = function(im){
BUUHWE_2D(im)
}

SHAH = function(im){
BUUHWE_2D(im)
}

BUUHWE_2D <- function(im) {

# Initialisation
      
	d <- dim(im)
	m <- d[1]
	n <- d[2]

	if ((m == 1) | (n == 1)) stop("both dimensions should be at least 2") 
	
	noe <- 2 * m * n - m - n     # number of edges

# Weights

	weights <- im - im + 1 

			#1 partout, de la meme taille que im


# NEW EDGES
edges <- matrix(0, 2*(n*m-1), 2)

	# Columns: From, to.
grid=matrix(1:(n*m),nrow=m)
tgrid=t(grid)


for (k in 1:(m*n-1)){
    edges[k,] =c(grid[k],grid[k+1])
    edges[(m*n-1)+k,] =c(tgrid[k],tgrid[k+1]) 
	}

RM_V = seq(m,m*n-1, by=m)
RM_H = seq(n,m*n-1, by=n)	+(m*n-1)

RM=c(RM_V,RM_H)
edges =edges[-RM,]

## modify order of the edges so as to coincide with Buuhwe_2D_Step

lengthH = (n*m-1)-length(RM_H)
selectHinput= edges[(nrow(edges)-lengthH+1):nrow(edges),1]
reorderingH = order(selectHinput)
lengthV = (n*m-1)-length(RM_V)
edges = edges[c((1:lengthV), lengthV+reorderingH ) ,]






  
# Decomposition history

	decomp.hist <- array(0, dim=c(3, 2, n*m-1))

	# Rows of each 3 x 2 matrix comprising decomp.hist: 
	# nodes of edge
	# filter coeffs used to decompose the nodes
	# (detail coeff, detail weight) of the decomposition


# Vectorised image and weights

	vec.im <- as.numeric((im))
	vec.weights <- as.numeric((weights))




# Decomposition

	for (s in 1:(n*m-1)) {

		a <- vec.weights[edges[,1]]
		b <- vec.weights[edges[,2]]

		h1 <- (1 + a^2 / b^2)^(-1/2)
		h2 <- -(1 + b^2 / a^2)^(-1/2)
		l1 <- -h2
		l2 <- h1

		details <- h1 * vec.im[edges[,1]] + h2 * vec.im[edges[,2]]

details =round(details,7)
## on evite les erreurs numeriques en arrondissant a la 7e decimale...


    		details.min.ind <- (which(abs(details) == min(abs(details))))[1]
		smooth.at.min <- l1[details.min.ind] * vec.im[edges[details.min.ind,1]] + 
				l2[details.min.ind] * vec.im[edges[details.min.ind,2]]

		det.weight.at.min <- h1[details.min.ind] * vec.weights[edges[details.min.ind,1]] + 
				h2[details.min.ind] * vec.weights[edges[details.min.ind,2]]
		sm.weight.at.min <- l1[details.min.ind] * vec.weights[edges[details.min.ind,1]] + 
				l2[details.min.ind] * vec.weights[edges[details.min.ind,2]]

		decomp.hist[1,,s] <- edges[details.min.ind,]
		decomp.hist[2,,s] <- c(h1[details.min.ind], h2[details.min.ind])
		decomp.hist[3,,s] <- c(details[details.min.ind], det.weight.at.min)



    eating.up <- edges[details.min.ind,][1]
    eaten.up <- edges[details.min.ind,][2]


		vec.im[eating.up] <- smooth.at.min
		vec.im[eaten.up] <- details[details.min.ind]

		vec.weights[eating.up] <- sm.weight.at.min
		vec.weights[eaten.up] <- det.weight.at.min

		edges[edges == eaten.up] <- eating.up
		edges <- edges[edges[,1] != edges[,2],]

		if (length(edges) == 2) edges <- matrix(edges, 1, 2)

	}

	return(list(d = dim(im), decomp.hist=decomp.hist))
#return(list(d = dim(im), decomp.hist=decomp.hist, vec.im=vec.im, vec.weights=vec.weights))

}
