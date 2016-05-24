#' Function to visualise an igraph object via arc diagram
#'
#' \code{visNetArc} is supposed to visualise a graph object of class "igraph" via arc diagram in one-dimensional layout. More precisely, it displays vertices (nodes) along an axis, with edges linked by arcs. With proper ordering of vertices (e.g. according to communities and degrees), arc diagram is able to identify clusters and bridges (as effective as two-dimensional layout). One advantage of using arc diagram is to allow for easy annotations along vertices.
#'
#' @param g an object of class "igraph"
#' @param orientation the orientation of the plots. It can be either "vertical" (default) or "horizontal"
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param ordering a numeric vector about the ordering of vertices. It is optional. It is highly recommend to order vertices according to communities and degrees
#' @param labels the label of the vertices. The default vertex labels are the name attribute of the nodes
#' @param vertex.label.color the color of vertex labels
#' @param vertex.label.cex the font size of vertex labels
#' @param vertex.color the fill color of the vertices. The default vertex colors are transparent
#' @param vertex.frame.color the color of the frame of the vertices. The default vertex frame colors are black
#' @param vertex.size the size of each vertex. By default, it is decided according to node degrees
#' @param vertex.pch the shape of each vertex. Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \url{http://www.statmethods.net/advgraphs/parameters.html}
#' @param vertex.lwd line width for the vertices (default 1)
#' @param edge.color the color of the edges (default "grey")
#' @param edge.width line width for the edges (default 1)
#' @param edge.lty line type for the edges (default 1)
#' @param ... additional graphic parameters associated with 'mtext'
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{visNet}}
#' @include visNetArc.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/80)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' g <- dNetInduce(g, V(g), knn=0)
#'
#' # 3) color nodes according to communities identified via a spin-glass model and simulated annealing
#' com <- spinglass.community(g, spins=4)
#' vgroups <- com$membership
#' palette.name <- visColormap(colormap="rainbow")
#' vcolors <- palette.name(length(com))[vgroups]
#' 
#' # 4) size nodes according to degrees
#' vdegrees <- igraph::degree(g)
#' 
#' # 5) sort nodes: first by communities and then degrees
#' tmp <- data.frame(ind=1:vcount(g), vgroups, vdegrees)
#' ordering <- tmp[order(vgroups,vdegrees),]$ind
#'
#' # 6) visualise graph using 1-dimensional arc diagram
#' visNetArc(g, ordering=ordering, labels=V(g)$name, vertex.label.color=vcolors, 
#' vertex.color=vcolors, vertex.frame.color=vcolors, vertex.size=log(vdegrees)+0.1)
#'
#' # 7) as comparison, also visualise graph on 2-dimensional layout 
#' visNet(g, colormap="bwr", layout=layout.kamada.kawai(g), vertex.label=V(g)$name, 
#' vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere")

visNetArc <- function(g, orientation=c('vertical','horizontal'), newpage=T, ordering=NULL, labels=V(g)$name, vertex.label.color="black", vertex.label.cex=1, vertex.color="transparent", vertex.frame.color="black", vertex.size=log(degree(g))+0.1, vertex.pch=21, vertex.lwd=1, edge.color="grey", edge.width=1, edge.lty=1, ...)
{

    orientation <- match.arg(orientation)
    
    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
    }
    if(is.null(V(g)$name)){
        V(g)$name <- as.character(V(g))
    }    
    nsize <- vcount(g)
    nedge <- ecount(g)

    edgelist <- get.edgelist(g)
    
    nodes <- V(g)$name
    nums <- seq_along(nodes)
    if (!is.null(labels)) {
        if (length(labels) != nsize){
            stop("\nLength of 'labels' differs from number of nodes")
        }
    }else{
        labels <- nodes
    }
    aux_ord <- 1:nsize
       
    if (!is.null(ordering)) {
        if (length(ordering) != nsize) 
            stop("\nLength of 'ordering' differs from number of nodes")
        nodes <- nodes[ordering]
        nums <- nums[ordering]
        labels <- labels[ordering]
        aux_ord <- ordering
    }
       
    ## for vertices (nodes)
    
    if (length(vertex.label.color) != nsize){
        vertex.label.color <- rep(vertex.label.color, length=nsize)
    }
    vertex.label.color <- vertex.label.color[aux_ord]
    
    if (length(vertex.label.cex) != nsize){
        vertex.label.cex <- rep(vertex.label.cex, length=nsize)
    }
    vertex.label.cex <- vertex.label.cex[aux_ord]
    
    if (length(vertex.color) != nsize){
        vertex.color <- rep(vertex.color, length=nsize)
    }
    vertex.color <- vertex.color[aux_ord]
    
    if (length(vertex.frame.color) != nsize){
        vertex.frame.color <- rep(vertex.frame.color, length=nsize)
    }
    vertex.frame.color <- vertex.frame.color[aux_ord]

    if (length(vertex.size) != nsize){
        vertex.size <- rep(vertex.size, length=nsize)
    }
    vertex.size <- vertex.size[aux_ord]

    if (length(vertex.pch) != nsize){
        vertex.pch <- rep(vertex.pch, length=nsize)
    }
    vertex.pch <- vertex.pch[aux_ord]
    
    if (length(vertex.lwd) != nsize){
        vertex.lwd <- rep(vertex.lwd, length=nsize)
    }
    vertex.lwd <- vertex.lwd[aux_ord]
    
    ## for edges (arcs)
    if (length(edge.color) != nedge){
        edge.color <- rep(edge.color, length=nedge)
    }
    if (length(edge.width) != nedge){
        edge.width <- rep(edge.width, length=nedge)
    }
    if (length(edge.lty) != nedge){
        edge.lty <- rep(edge.lty, length=nedge)
    }
    
    ####################################################
    
    nf <- rep(1/nsize, nsize)
    fin <- cumsum(nf)
    ini <- c(0, cumsum(nf)[-nsize])
    centers <- (ini + fin)/2
    tmp <- nodes
    e_num <- matrix(0, nrow(edgelist), 2)
    for (i in 1:nedge) {
        e_num[i, 1] <- centers[which(tmp == edgelist[i, 1])]
        e_num[i, 2] <- centers[which(tmp == edgelist[i, 2])]
    }
    radios <- abs(e_num[, 1] - e_num[, 2])/2
    max_radios <- which(radios == max(radios))
    max_rad <- unique(radios[max_radios]/2)
    locs <- rowSums(e_num)/2
    
    ######################################################################################
    ## Visualisation
    if (newpage){
        grDevices::dev.new()
    }
    if(orientation == "horizontal"){
        plot(0.5, 0.5, xlim = c(-0.015, 1.015), ylim = c(-0.01, 1*max_rad*2), type="n", xlab = "", ylab = "", axes=F)
        z <- seq(0, pi, l=100)
        for (i in 1:nedge){
            radio <- radios[i]
            x <- locs[i] + radio*cos(z)
            y <- radio*sin(z)
            graphics::lines(x, y, col=edge.color[i], lwd=edge.width[i], lty=edge.lty, lend=1, ljoin=2, lmitre=1)
        }
        if (!is.null(vertex.pch)){
            graphics::points(x=centers, y=rep(0, nsize), pch=vertex.pch, col=vertex.color, bg=vertex.frame.color, cex=vertex.size, lwd=vertex.lwd)
        }
        if (!is.null(labels)){
            graphics::mtext(labels, side=1, line=0, at=centers, cex=vertex.label.cex, outer=F, col=vertex.label.color, las=2, font=1, ...)
        }
    
    }else if(orientation == "vertical"){
        plot(0.5, 0.5, ylim = c(-0.015, 1.015), xlim = c(-0.01, 1*max_rad*2), type="n", xlab = "", ylab = "", axes=F)
        z <- seq(0, pi, l=100)
        
        for (i in 1:nedge){
            radio <- radios[i]
            y <- locs[i] + radio*cos(z)
            x <- radio*sin(z)
            graphics::lines(x, y, col=edge.color[i], lwd=edge.width[i], lty=edge.lty, lend=1, ljoin=2, lmitre=1)
        }
        
        if (!is.null(vertex.pch)){
            graphics::points(y=centers, x=rep(0, nsize), pch=vertex.pch, col=vertex.color, bg=vertex.frame.color, cex=vertex.size, lwd=vertex.lwd)
        }
        if (!is.null(labels)){
            graphics::mtext(labels, side=2, line=0, at=centers, cex=vertex.label.cex, outer=F, col=vertex.label.color, las=2, font=1, ...)
        }    
    }
    
    invisible()
}
