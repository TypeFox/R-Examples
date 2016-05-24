`transmatdualperm` <-
function (x, f, Pred = AdaptNeigh, neigh = 1, int = TRUE, clo = TRUE, 
    keep = 2, perm=sample(1:length(x),(length(x)-keep),FALSE), varonly = 
    FALSE) 
{

Wnew <- NULL
vars <- matrix(1, 1, length(x))
if (length(x) == length(unique(x))) {
        out <- fwtnpperm(x, f, LocalPred = Pred, neighbours = neigh, 
            intercept = int, closest = clo, nkeep = keep,mod=perm)
}
else {
        out <- fwtnpmp(x, f, LocalPred = Pred, neighbours = neigh, 
            intercept = TRUE, closest = clo, nkeep = keep, mpdet = "min")
}

n <- length(x)
Adual <- NULL
Tdual <- NULL
matno <- n - keep
pointsin <- out$pointsin
remlist <- out$removelist
newpoints <- (c(pointsin, rev(remlist)))
if (matno > 0) {
        Tdual <- Amatdual(matno, pointsin, remlist, out$neighbrs[[matno]], out$gamlist[[matno]], out$alphalist[[matno]])
        if (matno > 1) {
        	for (i in 2:matno) {
                	Adual <- Amatdual(matno - i + 1, pointsin, remlist, out$neighbrs[[matno - i + 1]], out$gamlist[[matno - i + 1]], out$alphalist[[matno - i + 1]])
                	augment <- rbind(cbind(Tdual, 0), 0)
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
        return(diag(Wnew))
}
else {
    return(list(out = out, Wnew = Wnew))           
}

}

