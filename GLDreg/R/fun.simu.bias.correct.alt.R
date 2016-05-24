fun.simu.bias.correct.alt <-
function (r, fit.obj) 
{
    n.simu <- nrow(r)
    index.var <- 2:(match("L1", dimnames(r)[[2]]) - 1)
    r.adj <- data.matrix(r[, index.var])
    r.adj.val <- do.call("rbind", lapply(1:n.simu, function(i, 
        r.adj, fit.obj, index.var) t(matrix(as.numeric(r.adj[i, 
        ])) - matrix(fit.obj[[3]][index.var - 1])), r.adj, fit.obj, 
        index.var))
    r.adj.val <- colMeans(r.adj.val)
    mode(r.adj) <- "numeric"
    r.adj.f <- r.adj - matrix(rep(r.adj.val, n.simu), ncol = ncol(r.adj), 
        byrow = T)
    return(r.adj.f)
}
