`DemoInfo` <-
function (mat) 
{
    eigs.mat <- eigen(mat)
    dom.pos <- which.max(abs(eigs.mat$values))
    L1 <- Re(eigs.mat$values[dom.pos])
    w.tmp <- eigs.mat$vectors[, dom.pos]
    w <- Re(w.tmp/sum(w.tmp))
    v.tmp <- eigen(t(mat))$vectors[, dom.pos]
    v <- Re(v.tmp/v.tmp[1])
    sens <- (v %*% t(w))/as.numeric(v %*% w)
    elas <- (mat/L1) * sens
    return(list(lambda = L1, SSD = w, RV = v, Sensitivities = sens, 
        Elasticities = elas, PPM = mat))
}
