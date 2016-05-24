`cltoSim` <-
function(cl){
        n <- length(cl)
        clvec <- c(rep(cl, n))
        clcomp <- rep(cl,each=n)
        matrix((clvec == clcomp)*1, ncol=n)
}

