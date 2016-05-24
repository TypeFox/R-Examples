niche.overlap.pair <- 
function(vectA, vectB, method = c("pianka", "schoener","petraitis","czech","morisita", "levins") ){ 
      method <- match.arg(method)
	  nij <- vectA
	  pij <- nij/sum(nij)
      
	  nkj <- vectB
	  pkj <- nkj/sum(nkj)
      switch(method,
        levins = {
            upper <- sum(pij * pkj)
	        lower <- sum(pij^2)
	        result <- upper/lower
        },
        schoener = {
            result <- 1 - sum(abs(pij - pkj))/2
        },
        petraitis = {
            nij <- nij[(vectA > 0) & (vectB > 0)]
            pij <- nij/sum(nij)
            nkj <- nkj[(vectA > 0) & (vectB > 0)]
            pkj <- nkj/sum(nkj)
            Eik <- sum(pij*log(pkj))-sum(pij*log(pij))
            result <- exp(Eik)
        },
        pianka = {
            result <- sum(pij*pkj)/sqrt((sum(pij^2)) * (sum(pkj^2)))
        },
        czech = {
            result <- 1 - (sum(abs(pij - pkj)))/2
        }, 
        morisita = {
            result <- 2*sum(pij*pkj)/((sum(pij^2)) + (sum(pkj^2)))
        }
      )
    return(result)
}

