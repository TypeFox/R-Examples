#' Function to visualise a graph object of class "igraph"
#'
#' \code{xVisNet} is supposed to visualise a graph object of class "igraph". It also allows vertices/nodes color-coded according to the input pattern. 
#'
#' @param g an object of class "igraph"
#' @param pattern a numeric vector used to color-code vertices/nodes. Notably, if the input vector contains names, then these names should include all node names of input graph, i.e. V(g)$name, since there is a mapping operation. After mapping, the length of the patern vector should be the same as the number of nodes of input graph; otherwise, this input pattern will be ignored. The way of how to color-code is to map values in the pattern onto the whole colormap (see the next arguments: colormap, ncolors, zlim and colorbar)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If pattern is null, it always sets to false
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.kamada.kawai"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param vertex.frame.color the color of the frame of the vertices. If it is NA, then there is no frame
#' @param vertex.size the size of each vertex. If it is a vector, each vertex may differ in size
#' @param vertex.color the fill color of the vertices. If it is NA, then there is no fill color. If the pattern is given, this setup will be ignored
#' @param vertex.shape the shape of each vertex. It can be one of "circle", "square", "csquare", "rectangle", "crectangle", "vrectangle", "pie" (\url{http://igraph.org/r/doc/vertex.shape.pie.html}), "sphere", and "none". If it sets to NULL, these vertices with negative will be "csquare" and the rest "circle". 
#' @param vertex.label the label of the vertices. If it is NA, then there is no label. The default vertex labels are the name attribute of the nodes
#' @param vertex.label.cex the font size of vertex labels.
#' @param vertex.label.dist the distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex.
#' @param vertex.label.color the color of vertex labels.
#' @param edge.arrow.size the size of the arrows for the directed edge. The default value is 1.
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{xSubneterGenes}}, \code{\link{xSubneterSNPs}}
#' @include xVisNet.r
#' @examples
#' # 1) generate a ring graph
#' g <- make_ring(10, directed=TRUE)
#'
#' # 2) visualise the graph
#' # 2a) visualise in one go
#' xVisNet(g=g, vertex.shape="sphere", glayout=layout_with_kk)
#' # 2b) visualise the graph with layout first calculated
#' glayout <- layout_(g, with_kk(), normalize(), component_wise())
#' xVisNet(g=g, vertex.shape="sphere", glayout=glayout)
#' # 2c) visualise the graph with layout appended to the graph itself
#' g <- add_layout_(g, with_kk(), normalize(), component_wise())
#' xVisNet(g=g, vertex.shape="sphere")
#'
#' # 4) visualise the graph with vertices being color-coded by the pattern
#' pattern <- runif(vcount(g))
#' names(pattern) <- V(g)$name
#' xVisNet(g=g, pattern=pattern, colormap="bwr", vertex.shape="sphere")

xVisNet <- function(g, pattern=NULL, colormap=c("yr","jet","gbr","wyr","br","bwr","rainbow","wb"), ncolors=40, zlim=NULL, colorbar=T, newpage=T, glayout=layout_nicely, vertex.frame.color=NA, vertex.size=NULL, vertex.color=NULL, vertex.shape=NULL, vertex.label=NULL, vertex.label.cex=NULL, vertex.label.dist=NULL, vertex.label.color="black", edge.arrow.size=0.8, ...)
{
    
    if (class(g) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }else{
    	ig <- g
    }
    
    if(is_directed(ig)){
    	#edge.arrow.size <- 0.8
    }
    
    if(!newpage){
    	if(dev.cur()>1){
    		#dev.off()
    	}
    }
    
    par_old <- graphics::par()
    
	dnet::visNet(g=ig, pattern=pattern, colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=colorbar, newpage=newpage, glayout=glayout, vertex.frame.color=vertex.frame.color, vertex.size=vertex.size, vertex.color=vertex.color, vertex.shape=vertex.shape, vertex.label=vertex.label, vertex.label.cex=vertex.label.cex, vertex.label.dist=vertex.label.dist, vertex.label.color=vertex.label.color, edge.arrow.size=edge.arrow.size, ...)
    
    suppressWarnings(graphics::par(par_old))
    
    invisible()
}
