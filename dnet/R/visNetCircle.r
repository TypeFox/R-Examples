#' Function to visualise an igraph object via circle diagram
#'
#' \code{visNetCircle} is supposed to visualise a graph object of class "igraph" via circle diagram. For better visualisation, ordering of vertices is determined according to communities and degrees.
#'
#' @param g an object of class "igraph"
#' @param com an object of class "communities" (see \url{http://igraph.org/r/doc/communities.html})
#' @param circles how circles are drawn in the plot. It can be either "single" for all communities being drawn in a single circle (by default) or "multiple" for communities being drawn in the different circles (i.e. one circle per community)
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param ordering a numeric vector about the ordering of vertices. It is optional. It is highly recommend to order vertices according to communities and degrees
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param vertex.label the label of the vertices. The default vertex labels are the name attribute of the nodes
#' @param vertex.size the size of each vertex. By default, it is decided according to node degrees
#' @param vertex.label.color the color of vertex labels
#' @param vertex.label.cex the font size of vertex labels
#' @param vertex.label.dist the distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex.
#' @param vertex.shape the shape of each vertex. It can be one of "circle", "square", "csquare", "rectangle", "crectangle", "vrectangle", "pie" (\url{http://igraph.org/r/doc/vertex.shape.pie.html}), "sphere", and "none". If it sets to NULL, these vertices with negative will be "csquare" and the rest "circle". 
#' @param edge.width line width for the edges (default 1)
#' @param edge.lty line type for the edges (default 1)
#' @param edge.color.within the color for edges within a community (default "grey")
#' @param edge.color.crossing the color for edges between communities (default "black")
#' @param mark.shape a numeric scalar or vector controlling the smoothness of the vertex group marking polygons. Its possible values are between -1 (fully polygons) and 1 (fully smoothness)
#' @param mark.expand a numeric scalar or vector, the size of the border around the marked vertex groups
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{visNet}}
#' @include visNetCircle.r
#' @examples
#' \dontrun{
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
#' mcolors <- palette.name(length(com))
#' vcolors <- mcolors[vgroups]
#' 
#' # 4) size nodes according to degrees
#' vdegrees <- igraph::degree(g)
#' 
#' # 5) sort nodes: first by communities and then degrees
#' tmp<-data.frame(ind=1:vcount(g), vgroups, vdegrees)
#' ordering <- tmp[order(vgroups,vdegrees),]$ind
#'
#' # 6) visualise graph using circle diagram
#' # 6a) drawn into a single circle 
#' visNetCircle(g=g, colormap="bwr", com=com, ordering=ordering)
#'
#' # 6b) drawn into multlpe circles (one circle per community) 
#' visNetCircle(g=g, colormap="bwr", com=com, circles="multiple", ordering=ordering)
#'
#' # 7) as comparison, also visualise graph on 2-dimensional layout 
#' mark.groups <- communities(com)
#' mark.col <- visColoralpha(mcolors, alpha=0.2)
#' mark.border <- visColoralpha(mcolors, alpha=0.2)
#' edge.color <- c("grey", "black")[crossing(com,g)+1]
#' visNet(g, colormap="bwr", glayout=layout.fruchterman.reingold, vertex.color=vcolors, 
#' vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, 
#' mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)
#' }

visNetCircle <- function(g, com, circles=c("single","multiple"), newpage=T, ordering=NULL, colormap=c("rainbow", "bwr","jet","gbr","wyr","br","yr","wb"), vertex.label=V(g)$name, vertex.size=log(igraph::degree(g))+2, vertex.label.color="black", vertex.label.cex=0.6, vertex.label.dist=0.75, vertex.shape="sphere", edge.width=1, edge.lty=1, edge.color.within="grey", edge.color.crossing="black", mark.shape=1, mark.expand=10, ...)
{

    circles <- match.arg(circles)

    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
    }
    if(is.null(V(g)$name)){
        V(g)$name <- as.character(V(g))
    }
    nsize <- vcount(g)
     
    if (class(com) != "communities"){
        stop("The function must apply to 'communities' object.\n")
    }
    
    if(!is.null(ordering)){
        ig <- graph.data.frame(d=get.data.frame(g,what="edges"), directed=F, vertices=as.data.frame(get.data.frame(g,what="vertices")[ordering,]))
        com$membership <- com$membership[ordering]
        com$names <- com$names[ordering]
        
        if (length(vertex.label) != nsize){
            vertex.label <- rep(vertex.label, length=nsize)
        }
        vertex.label <- vertex.label[ordering]
        
        if (length(vertex.size) != nsize){
            vertex.size <- rep(vertex.size, length=nsize)
        }
        vertex.size <- vertex.size[ordering]

        if (length(vertex.label.color) != nsize){
            vertex.label.color <- rep(vertex.label.color, length=nsize)
        }
        vertex.label.color <- vertex.label.color[ordering]
        
        if (length(vertex.label.cex) != nsize){
            vertex.label.cex <- rep(vertex.label.cex, length=nsize)
        }
        vertex.label.cex <- vertex.label.cex[ordering]
        
        if (length(vertex.label.dist) != nsize){
            vertex.label.dist <- rep(vertex.label.dist, length=nsize)
        }
        vertex.label.dist <- vertex.label.dist[ordering]
        
        if (length(vertex.shape) != nsize){
            vertex.shape <- rep(vertex.shape, length=nsize)
        }
        vertex.shape <- vertex.shape[ordering]
        
    }else{
        ig <- g
    }

    vgroups <- com$membership
    palette.name <- visColormap(colormap=colormap)
    mcolors <- palette.name(length(com))
    vcolors <- mcolors[vgroups]
    
    mark.groups <- communities(com)
    mark.col <- supraHex::visColoralpha(mcolors, alpha=0.2)
    mark.border <- supraHex::visColoralpha(mcolors, alpha=0.8)
    
    if(length(edge.color.within) & length(edge.color.crossing)){
        within_crossing <- c(edge.color.within,edge.color.crossing)
    }else{
        within_crossing <- c("grey", "black")
    }
    edge.color <- within_crossing[crossing(com,ig)+1]

    ## direction = 1 for clockwise; -1 for anti-clockwise
    radian.rescale <- function(x, start=0, direction=1) {
        c.rotate <- function(x){
            (x + start) %% (2 * pi) * direction
        }
        rescale <- function(x, to=c(0, 1), from=range(x, na.rm = TRUE)){
            (x-from[1])/diff(from)*diff(to)+to[1]
        }
        y <- rescale(x, c(0, 2*pi), range(x))
        c.rotate(y)
    }
    vertex.label.degree <- radian.rescale(x=1:vcount(ig), start=0, direction=-1)
    
    ######################################################################################
    if(circles == "single"){
        glayout <- layout.circle(ig)
    }else if(circles == "multiple"){
        graphs <- lapply(communities(com), function(x){
            dNetInduce(ig, nodes_query=V(ig)[unlist(x)]$name, knn=0, remove.loops=F, largest.comp=T)
        })
        layouts <- lapply(graphs, function(x) {
            layout.circle(x)
        })
        layouts_name <- lapply(graphs, function(x) {
            V(x)$name
        })
        lay <- layout.merge(graphs, layouts)
        rownames(lay) <- unlist(layouts_name) 
        glayout <- lay[V(ig)$name,]
    }
    ######################################################################################
    ## Visualisation
    
    visNet(ig, glayout=glayout, newpage=newpage, 
        vertex.label=vertex.label, 
        vertex.size=vertex.size,
        vertex.label.color=vertex.label.color,
        vertex.label.cex=vertex.label.cex,
        vertex.label.dist=vertex.label.dist,
        vertex.shape=vertex.shape,
        vertex.color=vcolors,
        vertex.frame.color=vcolors,
        vertex.label.degree=vertex.label.degree,
        edge.width=edge.width,
        edge.lty=edge.lty,
        edge.color=edge.color,
        mark.groups=mark.groups,
        mark.col=mark.col,
        mark.border=mark.border,
        mark.shape=mark.shape,
        mark.expand=mark.expand,
        ...)

    invisible()
}
