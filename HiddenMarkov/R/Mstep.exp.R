Mstep.exp <- function(x, cond, pm, pn){
    rate <- apply(cond$u, MARGIN=2, FUN=sum)/
          as.numeric(matrix(x, nrow=1) %*% cond$u)
    return(list(rate=rate))
}

