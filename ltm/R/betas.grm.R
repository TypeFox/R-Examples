betas.grm <-
function (thetas, constrained, ind1, ind2, p, trasform = TRUE) {
    betas <- vector("list", p)
    for (i in 1:p) {
        betas[[i]] <- if (constrained) {
            thets <- thetas[seq(ind1[i], ind2[i])]
            if (trasform)
                c(cumsum(c(thets[1], exp(thets[-1]))), abs(thetas[length(thetas)]))
            else
                c(thets, abs(thetas[length(thetas)]))
        } else {
            thets <- thetas[seq(ind1[i], ind2[i] - 1)]
            if (trasform) 
                c(cumsum(c(thets[1], exp(thets[-1]))), thetas[ind2[i]])
            else
                c(thets, thetas[ind2[i]])
        }
    }
    betas
}
