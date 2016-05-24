Mstep.pois <- function(x, cond, pm, pn){
    lambda <- as.numeric(matrix(x, nrow=1) %*% cond$u)/
                  apply(cond$u, MARGIN=2, FUN=sum)
    return(list(lambda=lambda))
}

