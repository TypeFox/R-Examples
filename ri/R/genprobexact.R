genprobexact <-
function(    Z,
                        blockvar = NULL,
                        clustvar = NULL) {

    desmat.out <- desmat.sanitize(Z,blockvar,clustvar)
    
    Z <- desmat.out$Z
    clustvar <- desmat.out$clustvar
    blockvar <- desmat.out$blockvar

    probs <- ave(Z, blockvar)
    probs.index <- c(1:length(probs))
    probsclus <- probs[match(clustvar, probs.index)]
    return(probsclus)
	}
