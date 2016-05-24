transmatdual2 <-
function (x, f, Pred = 5, neigh = 1, int = TRUE, clo = TRUE, 
    keep = 2, varonly = FALSE) 
{

Wnew <- NULL
if (length(x) == length(unique(x))) {
	out <- fwtnp2(x, f, LocalPred = Pred, neighbours = neigh, intercept = TRUE, closest = clo, nkeep = keep)
}
else {
        out <- fwtnpmp(x, f, LocalPred = Pred, neighbours = neigh, intercept = TRUE, closest = clo, nkeep = keep, mpdet = "min")
}

n <- length(x)
Adual <- NULL
Tdual <- NULL
matno <- n - keep
lca <- out$lca
pointsin <- out$pointsin
remlist <- lca[, 1]
newpoints <- (c(pointsin, rev(remlist)))
if (matno > 0) {
	nn<-lca[matno,2]
	Tdual <- Amatdual2(matno, pointsin, remlist,lca[matno,3:(2+nn)], lca[matno,(3+2*nn):(2+3*nn)],lca[matno,(3+nn):(2+2*nn)])
        if (matno > 1) {
		for (i in 2:matno) {
            		nn <- lca[matno-i+1, 2]
            		nbrs <- lca[matno-i+1, 3:(2 + nn)]
            		alpha <- lca[matno-i+1, (3 + nn):(2 + 2 * nn)]
            		weights <- lca[matno-i+1, (3 + 2 * nn):(2 + 3 * nn)]
            		Adual <- Amatdual2(matno-i+1, pointsin, remlist, nbrs, weights, alpha)
      		        augment <- rbind(cbind(Tdual, 0),  0)
      			augment[nrow(augment), nrow(augment)] <- 1
      			Tdual <- augment %*% Adual
        	}
      	}

        W <- Tdual
        Wnew <- matrix(0, length(x), length(x))
        reo <- NULL
        reo <- match(x, x[newpoints])
        Wnew<-W[reo,reo]
}
if (varonly) {
	return(diag(Wnew%*%t(Wnew)))
}
else {
	return(list(out = out, Wnew = Wnew))
}

}

