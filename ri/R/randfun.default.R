randfun.default <-
function(desmat.out) {

    desmat <- desmat.out$desmat
    Z <- desmat.out$Z
    clustvar <- desmat.out$clustvar
    blockvar <- desmat.out$blockvar

    perm <- rep(0, length(Z))
    for(b in 1:nrow(desmat)){
                block.tr <- sample(desmat[b,]$n,desmat[b,]$m)
                perm[blockvar==b][block.tr] <- 1
        }
    perm.index <- c(1:length(perm))
    permclus <- perm[match(clustvar, perm.index)]
    return(permclus)
}
