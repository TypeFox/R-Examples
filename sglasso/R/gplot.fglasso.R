gplot.fglasso <- function(object, rhoid, tp = c(1, 2), sub.tp1, sub.tp2, cex.sub = 1, k = 1.5, layout = layout.circle, ...){
	p <- object$p
	nrho <- object$nrho
    if(missing(rhoid)) stop("rhoid is not specified")
    if(length(rhoid) != 1) stop("rhois is not a scalar")
	if(!is.element(rhoid, 1:nrho)) stop("rhoid is not a valid index")
	rho <- object$rho[rhoid]
	if(length(tp) != 2) stop("length(tp) is not equal to 2")
	tp <- sort(tp)
	tp1 <- tp[1]
	tp2 <- tp[2]
	if(tp1 == tp2) stop("tp1 is equal to tp2")
	if(!is.element(tp1, 1:object$tp)) stop("tp1 is not a correct time point")
	if(!is.element(tp2, 1:object$tp)) stop("tp2 is not a correct time point")
	if(missing(sub.tp1)) sub.tp1 <- paste("Network at time point ", tp1)
	if(missing(sub.tp2)) sub.tp2 <- paste("Network at time point ", tp2)
	Theta.hat <- Kh(object, rho)[[1]]
	id.tp1 <- 1:p + p * (tp1 - 1)
	id.tp2 <- 1:p + p * (tp2 - 1)
	adjmat.tp1 <- as.matrix(abs(Theta.hat[id.tp1, id.tp1]) > 0)
	g.tp1 <- graph.adjacency(adjmat.tp1, mode = "undirected", diag = FALSE)
	adjmat.tp2 <- as.matrix(abs(Theta.hat[id.tp2, id.tp2]) > 0)
	g.tp2 <- graph.adjacency(adjmat.tp2, mode = "undirected", diag = FALSE)
	adjmat.net <- diag(2 * p)
	adjmat.net[1:p, p + 1:p] <- as.matrix(abs(Theta.hat[id.tp1, id.tp2]) > 0)
	g.net <- graph.adjacency(adjmat.net, mode = "directed", diag = FALSE)
	if(is.function(layout)) g.layout <- layout(g.tp1)
	if(is.matrix(layout)){
		if((dim(layout)[1] != p) | (dim(layout)[2] != 2))
		stop("wrong dimension in layout matrix")
		g.layout <- layout
	}
	g.tp1$layout <- g.layout
	g.tp2$layout <- g.layout
	g.tp1$layout[, 1] <- g.tp1$layout[, 1] - k
	g.tp2$layout[, 1] <- g.tp2$layout[, 1] + k
	if(!is.element("vertex.label", names(list(...)))) V(g.tp1)$label <- V(g.tp2)$label <- 1:p
	E(g.tp1)$color <- E(g.tp2)$color <- "black"
	g.net$layout <- rbind(g.tp1$layout, g.tp2$layout)
	V(g.net)$label <- NA
	E(g.net)$color <- "darkgray"
	E(g.net)$lty <- 2
	xlim <- range(c(g.tp1$layout[, 1], g.tp2$layout[, 1]))
	ylim <- range(c(g.tp1$layout[, 2], g.tp2$layout[, 2]))
	plot.igraph(g.net, rescale = FALSE, xlim = xlim, ylim = ylim, edge.curved = 0.1, ...)
	plot.igraph(g.tp1, rescale = FALSE, xlim = xlim, add = TRUE, ...)
	plot.igraph(g.tp2, rescale = FALSE, xlim = xlim, add = TRUE, ...)
	x <- mean(range(g.tp1$layout[, 1]))
	y <- min(g.tp1$layout[, 2]) - 0.5
	text(x = x, y = y, labels = sub.tp1, cex = cex.sub)
	x <- mean(range(g.tp2$layout[, 1]))
	text(x = x, y = y, labels = sub.tp2, cex = cex.sub)
	out <- list(graph.tp1 = g.tp1, graph.tp2 = g.tp2, graph.net = g.net, layout = g.layout, ...)
	invisible(out)
}
