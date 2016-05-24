fwtnpperm <-function (input, f, LocalPred = LinearPred, neighbours = 1,
intercept = TRUE,closest = FALSE, nkeep = 2, initboundhandl = "reflect", mod = 
    sample(1:length(input), (length(input) - nkeep), FALSE), 
    do.W = FALSE, varonly = FALSE) 
{

X <- input
I <- intervals(X, initboundhandl)
lengths <- lengthintervals(X, I, type = "midpoints", neighbours, closest)
X <- as.row(X)
f <- as.row(f)
nkeep <- max(nkeep, 1)
n <- length(X)
removelist <- NULL
lengthsremove <- NULL
neighbrs <- list()
gamlist <- list()
alphalist <- list()
schemehist <- NULL
interhist <- NULL
clolist <- NULL
pointsin <- matrix(1:n, 1, n)
pointsin <- pointsin[order(X)]
coeff <- f
matno <- n - nkeep
W <- v<-NULL
   if ((do.W == 1) & (varonly == 1)) {
        varonly <- FALSE
    }
    ex <- do.W + varonly
    if (ex == 1) {
        W <- diag(n)
    }
    if (varonly) {
        v <- rep(1, times = n)
    }

for (j in 1:matno) {
        remove <- mod[j]
        removelist[j] <- remove
        out <- getnbrs(X, remove, pointsin, neighbours, closest)
        nbrs <- out$n
        index <- out$index
        res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours)
        if (length(res) == 2) {
        	l <- res[[1]]
        	clolist[j] <- res[[2]][[1]]
        	nbrs <- res[[2]][[2]]
        	index <- res[[2]][[3]]
        }
        else {
        	l <- res
        }
        neighbrs[[j]] <- nbrs
        weights <- l[[1]]
        pred <- l[[2]]
        if (length(l) == 3) {
        	scheme <- NULL
        	int <- NULL
        	details <- NULL
        }
        else {
        	scheme <- l[[5]]
        	int <- l[[4]]
        	details <- l[[6]]
        }
        coeff[remove] <- coeff[remove] - pred
        l1 <- PointsUpdate(X, coeff, nbrs, index, remove, pointsin, weights, lengths)
        coeff <- l1$coeff
        lengths <- l1$lengths
        r <- l1$r
        weights <- l1$weights
        N <- l1$N
        alpha <- l1$alpha
        if (ex) {
	        if(varonly){
	        	W[r, ] <- W[r, ] - colSums(as.vector(weights) * matrix(W[index, ], nrow = length(nbrs)))
	        	W[index, ] <- W[index, ] + matrix(alpha) %*% W[r, ]
	            	v[remove]<-sum(W[r,]^2)
	            	np<-setdiff(1:length(pointsin),r)
		    	W<-W[np,]	
	        }
	        else{
	            	W[remove, ] <- W[remove, ] - colSums(as.vector(weights) * matrix(W[nbrs, ], nrow = length(nbrs)))
	            	W[nbrs, ] <- W[nbrs, ] + matrix(alpha) %*% W[remove, ]
	        }
        }
        lengthsremove[j] <- lengths[r]
        gamlist[[j]] <- weights
        alphalist[[j]] <- alpha
        schemehist[j] <- scheme
        interhist[j] <- int
        lengths <- lengths[setdiff(1:length(pointsin), r)]
        pointsin <- setdiff(pointsin, remove)
}
if(varonly){
	v[pointsin]<-rowSums(W^2)
    	W<-NULL
}
    
    N <- length(pointsin)

    return(list(x=input,coeff = coeff,lengths = lengths, lengthsremove = lengthsremove, pointsin = pointsin, 
        removelist = removelist, neighbrs = neighbrs, neighbours = neighbours, 
        schemehist = schemehist, interhist = interhist, clolist = clolist, 
        gamlist = gamlist, alphalist = alphalist, W = W,v=v))
}

