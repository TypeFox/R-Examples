#' Function to visualise a graph object of class "igraph" or "graphNEL"
#'
#' \code{visNet} is supposed to visualise a graph object of class "igraph" or "graphNEL". It also allows the color-coding of vertices by providing the input pattern. 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param pattern a numeric vector used to color-code vertices/nodes. Notably, if the input vector contains names, then these names should include all node names of input graph, i.e. V(g)$name, since there is a mapping operation. After mapping, the length of the patern vector should be the same as the number of nodes of input graph; otherwise, this input pattern will be ignored. The way of how to color-code is to map values in the pattern onto the whole colormap (see the next arguments: colormap, ncolors, zlim and colorbar)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If pattern is null, it always sets to false
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout.auto", "layout.random", "layout.circle", "layout.sphere", "layout.fruchterman.reingold", "layout.kamada.kawai", "layout.spring", "layout.reingold.tilford", "layout.fruchterman.reingold.grid", "layout.lgl", "layout.graphopt", "layout.svd" and "layout.norm". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param vertex.frame.color the color of the frame of the vertices. If it is NA, then there is no frame
#' @param vertex.size the size of each vertex. If it is a vector, each vertex may differ in size
#' @param vertex.color the fill color of the vertices. If it is NA, then there is no fill color. If the pattern is given, this setup will be ignored
#' @param vertex.shape the shape of each vertex. It can be one of "circle", "square", "csquare", "rectangle", "crectangle", "vrectangle", "pie" (\url{http://igraph.org/r/doc/vertex.shape.pie.html}), "sphere", and "none". If it sets to NULL, these vertices with negative will be "csquare" and the rest "circle". 
#' @param vertex.label the label of the vertices. If it is NA, then there is no label. The default vertex labels are the name attribute of the nodes
#' @param vertex.label.cex the font size of vertex labels.
#' @param vertex.label.dist the distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex.
#' @param vertex.label.color the color of vertex labels.
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{dNetFind}}
#' @include visNet.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 3) visualise the subg with vertices being color-coded by the pattern
#' pattern <- runif(vcount(subg))
#' names(pattern) <- V(subg)$name
#' visNet(g=subg, pattern=pattern, colormap="bwr", vertex.shape="sphere")

visNet <- function(g, pattern=NULL, colormap=c("bwr","jet","gbr","wyr","br","yr","rainbow","wb"), ncolors=40, zlim=NULL, colorbar=T, newpage=T, glayout=layout.fruchterman.reingold, vertex.frame.color=NA, vertex.size=NULL, vertex.color=NULL, vertex.shape=NULL, vertex.label=NULL, vertex.label.cex=NULL, vertex.label.dist=NULL, vertex.label.color="black", ...)
{
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }    
    nsize <- vcount(ig)
    
    if(is.null(vertex.label)){
        if ("geneSymbol" %in% list.vertex.attributes(ig)){
            vertex.label <- V(ig)$geneSymbol
        }else{
            vertex.label <- V(ig)$name
        }
    }
    
    shapes <- rep("circle", vcount(ig))
    if(is.null(vertex.shape)){
        if ("score" %in% list.vertex.attributes(ig)){
            shapes[V(ig)$score<0] <- "csquare"
        }
        vertex.shape <- shapes
    }
    
    ######################################################################################
    if (!is.null(pattern)){
    
        flag <- 0
        if(!is.null(names(pattern))){
            pattern <- pattern[V(ig)$name]
        }
        if(length(pattern)==nsize){
            flag <- 1
        }
                
        if(flag==1){
        	
        	pattern <- as.numeric(pattern)
        	
            if(is.null(zlim)){
                vmin <- floor(stats::quantile(pattern, 0.05))
                vmax <- ceiling(stats::quantile(pattern, 0.95))
                if(vmin < 0 & vmax > 0){
                    vsym <- abs(min(vmin, vmax))
                    vmin <- -1*vsym
                    vmax <- vsym
                }
                zlim <- c(vmin,vmax)
            }
            
            ## A function to map a vector to colors
            vec2color <- function(vec, colormap=colormap, ncolors=ncolors, zlim=zlim){
                palette.name <- visColormap(colormap=colormap)
                colors <- palette.name(ncolors)
                scale <- length(colors)/(max(zlim)-min(zlim))
                sapply(1:length(vec), function(x){
                    ind <- floor(1+(vec[x]-min(zlim))*scale)
                    colors[max(1,min(ncolors,ind))]
                })
            }
            vertex.color <- vec2color(pattern, colormap=colormap, ncolors=ncolors, zlim=zlim)
        }else{
            warning("The input 'pattern' is ignored. Please check the help for enabling your input")
            pattern <- NULL
            if(is.null(vertex.color)){
                vertex.color <- "SkyBlue2"
            }
        }
    }else{
        if(is.null(vertex.color)){
            vertex.color <- "SkyBlue2"
        }
    }
    ######################################################################################
    
    max.labels <- max(nchar(vertex.label))
    vertex.size2 <- 15
    vertex.label.dist2 <- 0 
    vertex.label.cex2 <- 0.6
    if (nsize < 50){
        if (max.labels > 2){
            vertex.size2 <- 8
            vertex.label.dist2 <- 0.5   
        }
    }
    if (nsize < 100 && nsize >= 50){
        if (max.labels > 2){
            vertex.size2 <- 8
            vertex.label.dist2 <- 0.5   
        }
        vertex.label.cex2 <- 0.5
    }
    if (nsize >= 100){
        if (max.labels > 3){
            vertex.size2 <- 8
            vertex.label.dist2 <- 0.5
        }
        vertex.label.cex2 <- 0.4
    }
    
    if(is.null(vertex.size)){
        vertex.size <- vertex.size2
    }
    if(is.null(vertex.label.dist)){
        vertex.label.dist <- vertex.label.dist2
    }
    if(is.null(vertex.label.cex)){
        vertex.label.cex <- vertex.label.cex2
    }
    
    ######################################################################################
    ## Visualisation
    if (newpage){
        grDevices::dev.new()
    }
    
    plot.igraph(ig, layout=glayout, 
        vertex.frame.color=vertex.frame.color,
        vertex.size=vertex.size, 
        vertex.color=vertex.color,
        vertex.shape=vertex.shape,
        vertex.label=vertex.label,
        vertex.label.cex=vertex.label.cex, 
        vertex.label.dist=vertex.label.dist, 
        vertex.label.color=vertex.label.color, 
        vertex.label.family="sans",
        ...)
        
    ## colorbar
    if (!is.null(pattern) && length(pattern)==nsize){
        if(colorbar){
            graphics::par(fig=c(0,0.1,0.5,1), new=TRUE)
            palette.name <- visColormap(colormap=colormap)
            colors <- palette.name(ncolors)
            lab.scale <- length(colors)/(zlim[2]-zlim[1])
            for (i in 1:length(colors)) {
                yValue <- (i-1)/ncolors
                hValue <- 1/ncolors
                xValue <- 0
                wValue <- 0.5
        
                ## for rect
                xleft <- xValue
                ybottom <- yValue
                xright <- xValue+wValue
                ytop <- yValue+hValue
                graphics::rect(xleft,ybottom,xright,ytop, col=colors[i], border="transparent")
        
                if(i == 1 | i == 1+length(colors)/2){
                    tx <- (i-1)/lab.scale + zlim[1]
                    graphics::text(x=xright*2, y=ybottom, labels=tx, cex=0.6)
                }else if(i==length(colors)){
                    tx <- i/lab.scale + zlim[1]
                    graphics::text(x=xright*2, y=ytop, labels=tx, cex=0.6)
                }
            }
        }
    }
    
    invisible()
}