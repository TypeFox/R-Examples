gplot.sglasso <- function(object, rhoid, layout = layout.circle, ...){
    nrho <- object$nrho
    if(missing(rhoid)) rhoid <- 1:nrho
	if(!all(is.element(rhoid, 1:nrho))) stop("rhoid is not a valid index")
    rho <- object$rho[rhoid]
    Theta.hat <- Kh(object, rho)
    adjmat <- lapply(Theta.hat, function(x) as.matrix(abs(x) > 1.0e-13))
    rhoid.max <- which.max(unlist(lapply(adjmat, sum)))
    g <- graph.adjacency(adjmat[[rhoid.max]], mode = "undirected", diag = FALSE)
    if(is.function(layout)) g.layout <- layout(g)
    p <- dim(Theta.hat[[1]])[1]
    if(is.matrix(layout)){
        if((dim(layout)[1] != p) | (dim(layout)[2] != 2))
        stop("wrong dimension in layout matrix")
        g.layout <- layout
    }
    g <- graph.adjacency(adjmat[[1]], mode = "undirected", diag = FALSE)
    g$layout <- g.layout
    if(!is.element("vertex.label", names(list(...)))) V(g)$label <- 1:p
    sub <- substitute(paste(rho[i], " = ", v), list(i = rhoid[1], v = rho[1]))
    plot.igraph(g, sub = sub, ...)
    op <- par(ask=dev.interactive())
    if(length(rhoid) > 1){
        for(i in 2:length(rhoid)){
            if(!all(as.vector(adjmat[[i - 1]] == adjmat[[i]]))){
                g <- graph.adjacency(adjmat[[i]], mode = "undirected", diag = FALSE)
                g$layout <- g.layout
                if(!is.element("vertex.label", names(list(...)))) V(g)$label <- 1:p
                sub <- substitute(paste(rho[i], " = ", v), list(i = rhoid[i], v = rho[i]))
                plot.igraph(g, sub = sub, ...)
                op <- par(ask=dev.interactive())
            }
        }
    }
    par(op)
}
