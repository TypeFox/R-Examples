#' Function to visualise a direct acyclic graph (DAG) with node colorings according to a named input data vector (if provided)
#'
#' \code{visDAG} is supposed to visualise a direct acyclic graph (DAG) with node colorings according to a named input data vector (if provided)
#'
#' @param g an object of class "igraph"
#' @param data a named input data verctor used to color-code vertices/nodes. The input data vector must have names, and these names should include all node names of input graph, i.e. V(g)$name, since there is a mapping operation. The way of how to color-code is to map values in the data onto the whole colormap (see the next arguments: colormap, ncolors, zlim and colorbar)
#' @param height a numeric value specifying the height of device
#' @param width a numeric value specifying the width of device
#' @param margin margins as units of length 4 or 1
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/data values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If data is null, it always sets to false
#' @param colorbar.fraction the relative fraction of colorbar block against the device size
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param layout.orientation the orientation of the DAG layout. It can be one of "left_right" for the left-right layout (viewed from the DAG root point), "top_bottom" for the top-bottom layout, "bottom_top" for the bottom-top layout, and "right_left" for the right-left layout
#' @param node.info tells the ontology term information used to label nodes. It can be one of "none" for no node labeling, "term_id" for using Term ID, "term_name" for using Term Name (the first 15 characters), "both" for using both of Term ID and Name (the first 15 characters), and "full_term_name" for using the full Term Name
#' @param graph.node.attrs a list of global node attributes. These node attributes will be changed globally. See 'Note' below for details on the attributes
#' @param graph.edge.attrs a list of global edge attributes. These edge attributes will be changed globally. See 'Note' below for details on the attributes
#' @param node.attrs a list of local edge attributes. These node attributes will be changed locally; as such, for each attribute, the input value must be a named vector (i.e. using Term ID as names). See 'Note' below for details on the attributes
#' @return 
#' An object of class 'Ragraph'
#' @note
#' A list of global node attributes used in "graph.node.attrs":
#' \itemize{
#' \item{"shape": the shape of the node: "circle", "rectangle", "rect", "box" and "ellipse"}
#' \item{"fixedsize": the logical to use only width and height attributes. By default, it sets to true for not expanding for the width of the label}
#' \item{"fillcolor": the background color of the node}
#' \item{"color": the color for the node, corresponding to the outside edge of the node}
#' \item{"fontcolor": the color for the node text/labelings}
#' \item{"fontsize": the font size for the node text/labelings}
#' \item{"height": the height (in inches) of the node: 0.5 by default}
#' \item{"width": the width (in inches) of the node: 0.75 by default}
#' \item{"style": the line style for the node: "solid", "dashed", "dotted", "invis" and "bold"}
#' }
#' A list of global edge attributes used in "graph.edge.attrs":
#' \itemize{
#' \item{"color": the color of the edge: gray by default}
#' \item{"weight": the weight of the edge: 1 by default}
#' \item{"style": the line style for the edge: "solid", "dashed", "dotted", "invis" and "bold"}
#' }
#' A list of local node attributes used in "node.attrs" (only those named Term IDs will be changed locally!):
#' \itemize{
#' \item{"label": a named vector specifying the node text/labelings}
#' \item{"shape": a named vector specifying the shape of the node: "circle", "rectangle", "rect", "box" and "ellipse"}
#' \item{"fixedsize": a named vector specifying whether it sets to true for not expanding for the width of the label}
#' \item{"fillcolor": a named vector specifying the background color of the node}
#' \item{"color": a named vector specifying the color for the node, corresponding to the outside edge of the node}
#' \item{"fontcolor": a named vector specifying the color for the node text/labelings}
#' \item{"fontsize": a named vector specifying the font size for the node text/labelings}
#' \item{"height": a named vector specifying the height (in inches) of the node: 0.5 by default}
#' \item{"width": a named vector specifying the width (in inches) of the node: 0.75 by default}
#' \item{"style": a named vector specifying the line style for the node: "solid", "dashed", "dotted", "invis" and "bold"}
#' }
#' @export
#' @importFrom graph edgeData
#' @importFrom Rgraphviz getDefaultAttrs agopen plot
#' @seealso \code{\link{dDAGreverse}}, \code{\link{dDAGroot}}, \code{\link{dDAGinduce}}, \code{\link{dDAGlevel}}
#' @include visDAG.r
#' @examples
#' \dontrun{
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) randomly select vertices as the query nodes
#' # the more common, the query nodes can be term id
#' nodes_query <- V(g)[sample(V(g),5)]$name
#'
#' # 3) obtain the induced subgraph based on all possible paths
#' subg <- dDAGinduce(g, nodes_query, path.mode="all_paths")
#'
#' # 4) just visualise the induced subgraph
#' visDAG(g=subg, node.info="both")
#'
#' # 5) color-code nodes/terms according to its level
#' data <- dDAGlevel(subg)
#' visDAG(g=subg, data=data, node.info="both")
#' # 5a) globally change the node and edge attributes
#' visDAG(g=subg, data=data, layout.orientation="top_bottom", node.info="both", graph.node.attrs=list(fixedsize=FALSE,shape="box",color="transparent"), graph.edge.attrs=list(color="black"))
#' # 5b) locally highlight the root by changing its shape into "box"
#' root <- dDAGroot(subg)
#' root.shape <- "box"
#' names(root.shape) <- V(subg)[root]$name
#' visDAG(g=subg, data=data, node.info="both", node.attrs=list(shape=root.shape))
#' # 5c) further locally remove the root labelling
#' root.label <- ""
#' names(root.label) <- V(subg)[root]$name
#' visDAG(g=subg, data=data, node.info="both", node.attrs=list(shape=root.shape,label=root.label))
#' }

visDAG <- function (g, data=NULL, height=7, width=7, margin=rep(0.1,4), colormap=c("yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.1, newpage=T, layout.orientation=c("left_right","top_bottom","bottom_top","right_left"), node.info=c("none", "term_id", "term_name", "both", "full_term_name"), graph.node.attrs=NULL, graph.edge.attrs=NULL, node.attrs=NULL)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    layout.orientation <- match.arg(layout.orientation)
    node.info<- match.arg(node.info)
    
    ## check input graph
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }

    ######################################################################################
    
    if (!is.null(data) & !is.null(names(data))){
    
        ## check input data
        data <- data[!is.na(names(data))]
    
        ## check mapping between input data and graph
        flag <- 1
        ind <- match(names(data), V(ig)$name)
        nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
        if(length(nodes_mapped)!=vcount(ig)){
            #stop("The function must require that the row names of input data contain all those in the input graph.\n")
            #flag <- 0
        }
        data <- data[nodes_mapped]
        
        if(flag==1){
            ## determine the color range
            if(is.null(zlim)){
                vmin <- floor(stats::quantile(data, 0.05))
                vmax <- ceiling(stats::quantile(data, 0.95))
                if(vmin < 0 & vmax > 0){
                    vsym <- abs(min(vmin, vmax))
                    vmin <- -1*vsym
                    vmax <- vsym
                }
                if(!is.null(zlim)){
                    if(zlim[1] < floor(min(data)) | zlim[2] > ceiling(max(data))){
                        #zlim <- c(vmin,vmax)
                    }
                }else{
                    zlim <- c(vmin,vmax)
                }
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
            node.fillcolor <- vec2color(data, colormap=colormap, ncolors=ncolors, zlim=zlim)
            names(node.fillcolor) <- names(data)
        }else{
            warning("The input 'data' is ignored. Please check the help for enabling your input")
            data <- NULL

        }
    }else{
        data <- NULL
    }
    
    ######################################################################################
    if(layout.orientation=="bottom_top" || layout.orientation=="right_left"){
        ig <- dDAGreverse(ig)
    }
    
    ## convert from igraph object into graph object
    dag <- igraph.to.graphNEL(ig)
    
    ## define node labels
    getTermInfo <- function(g, vids, numChar=15, mulLines=F){
        fullNames <- V(g)[vids]$term_name
        names(fullNames) <- V(g)[vids]$name
    
        if(mulLines==F){
            shortNames <- paste(substr(fullNames,1,numChar), ifelse(nchar(fullNames)>numChar, '...', ''), sep='')
        }else{
            shortNames <- sapply(fullNames,function(x){
                return(paste(strwrap(x, numChar), sep="", collapse = "\\\n"))
            })
        }
    
        names(shortNames) <- names(fullNames)
        return(shortNames)
    }
    
    termNames <- getTermInfo(ig, vids=dag@nodes, numChar=15, mulLines=F)
    nodeInfo <- switch(node.info,
                       none = NULL,
                       term_id = dag@nodes,
                       term_name = termNames,
                       both = paste(dag@nodes, termNames, sep="\\\n"),
                       full_term_name = getTermInfo(ig, vids=dag@nodes, numChar=30, mulLines=T)
                       )

    ########################################################################

    ## a function to update the parameters
    update_attributes <- function(default, update){
        for(i in 1:length(names(default))){
            item <- names(default)[i]
            if(item %in% names(update)){
                default[i] <- update[which(names(update)==item)]
            }
        }
        return(default)
    }


    ## a function to update the parameters (vectorised)
    update_attributes_vectorised <- function(default, update){
        for(i in 1:length(names(default))){
            item <- names(default)[i]
            if(item %in% names(update)){
                new <- update[[which(names(update)==item)]]
                old <- default[[i]]
                old[names(new)] <- new
                default[[i]] <- old
            }
        }
        return(default)
    }

    ################
    #graph.node.attrs <- list(shape="circle",color="transparent")
    ## the default parameters for "graph.node.attrs"
    graph.node.attrs.default <- list(
                    shape = c("ellipse","box","circle","plaintext")[1],
                    fixedsize = "TRUE",
                    fillcolor = "transparent",
                    color = "gray",
                    fontcolor = "black",
                    fontsize = 14,
                    height = 0.5,
                    width = 0.75,
                    style = "solid"
                    )
    ## update parameters for "graph.node.attrs"
    graph.node.attrs.default <- update_attributes(graph.node.attrs.default, graph.node.attrs)
    
    #graph.edge.attrs <- list(color="black")
    ## the default parameters for "graph.edge.attrs"
    graph.edge.attrs.default <- list(
                    color = "gray",
                    weight = 1,
                    style = c("solid", "dashed", "dotted", "invis", "bold")[1]
                    )
    ## update parameters for "graph.edge.attrs"
    graph.edge.attrs.default <- update_attributes(graph.edge.attrs.default, graph.edge.attrs)
    
    ########################################################################

    flag <- 0
    if(length(grDevices::dev.list())==0){
        flag <- 1
    }
    ## global Graphviz attributes
    #opar <- graphics::par() # make a copy of current settings
    graphics::par("fin"=c(0.69,0.69))
    graphAttrs <- Rgraphviz::getDefaultAttrs(curAttrs=list(), layoutType='dot')
    
    if(flag==1){
        grDevices::dev.off(which=grDevices::dev.cur())
    }
    
    graphAttrs$cluster <- NULL
    # http://svitsrv25.epfl.ch/R-doc/library/Rgraphviz/html/GraphvizAttributes.html
    graphAttrs$graph$rankdir <- ifelse(layout.orientation=="left_right" || layout.orientation=="right_left", "LR", "TB")
    ## set the node attributes
    graphAttrs$node$shape <- graph.node.attrs.default$shape
    graphAttrs$node$fixedsize <- graph.node.attrs.default$fixedsize
    graphAttrs$node$fillcolor <- graph.node.attrs.default$fillcolor
    graphAttrs$node$color <- graph.node.attrs.default$color
    graphAttrs$node$fontcolor <- graph.node.attrs.default$fontcolor
    graphAttrs$node$fontsize <- graph.node.attrs.default$fontsize
    graphAttrs$node$height <- graph.node.attrs.default$height
    graphAttrs$node$width <- graph.node.attrs.default$width
    graphAttrs$node$style <- graph.node.attrs.default$style
    ## set the edge attributes
    graphAttrs$edge$color <- graph.edge.attrs.default$color
    graphAttrs$edge$weight <- graph.edge.attrs.default$weight
    graphAttrs$edge$style <- graph.edge.attrs.default$style
    
    #################################
    ## local node attributes lists
    node.attrs.default <- list()
    
    ########################
    ## set the local node attributes    
    if(is.null(data)){
        node.fillcolor <- rep(graph.node.attrs.default$fillcolor, length(dag@nodes))
    }else{
        ## those nodes without data will be assigned with default colors
        tmp <- rep(graph.node.attrs.default$fillcolor, length(dag@nodes))
        names(tmp) <- dag@nodes
        tmp[names(node.fillcolor)] <- node.fillcolor
        node.fillcolor <- tmp
    }
    if(is.null(nodeInfo)) {
        node.label <- character(length(dag@nodes))
    }else{
        node.label <- nodeInfo
    }
    node.shape <- rep(graph.node.attrs.default$shape, length(dag@nodes))
    node.fixedsize <- rep(graph.node.attrs.default$fixedsize, length(dag@nodes))
    node.color <- rep(graph.node.attrs.default$color, length(dag@nodes))
    node.fontcolor <- rep(graph.node.attrs.default$fontcolor, length(dag@nodes))
    node.fontsize <- rep(graph.node.attrs.default$fontsize, length(dag@nodes))
    node.height <- rep(graph.node.attrs.default$height, length(dag@nodes))
    node.width <- rep(graph.node.attrs.default$width, length(dag@nodes))
    node.style <- rep(graph.node.attrs.default$style, length(dag@nodes))
    
    names(node.fillcolor) <- names(node.label) <- names(node.shape) <- names(node.fixedsize) <- names(node.color) <- names(node.fontcolor) <- names(node.fontsize) <- names(node.height) <- names(node.width) <- names(node.style) <- dag@nodes
    
    #user.node.shape <- "box"
    #nodeRoot <- V(ig)[which(V(g)$term_distance==0)]$name
    #names(user.node.shape) <- nodeRoot
    #node.attrs <- list(shape=user.node.shape)
    ## the default parameters for "node.attrs"
    node.attrs.default <- list(
                    fillcolor = node.fillcolor,
                    label = node.label,
                    shape = node.shape,
                    fixedsize =node.fixedsize,
                    color = node.color,
                    fontcolor = node.fontcolor,
                    fontsize = node.fontsize,
                    height = node.height,
                    width = node.width,
                    style = node.style
                    )
    ## update parameters for "node.attrs.default"
    node.attrs.default <- update_attributes_vectorised(node.attrs.default, node.attrs)
    
    ########################
    ## local edge attributes lists
    edgeAttrs <- list()
    if(class(suppressWarnings(try(graph::edgeData(dag,attr="relation"), T)))!="try-error"){
        edgeAttrs$color <- ifelse(graph::edgeData(dag,attr="relation")=="is_a", "black", "red")
        edgeAttrs$style <- ifelse(graph::edgeData(dag,attr="relation")=="is_a", "dashed", "dashed")
    }
    
    # graph: An object of class 'graphNEL'
    # name: The name of the graph
    # attrs: A list of graphviz attributes
    # nodeAttrs: A list of specific node attributes
    # edgeAttrs: A list of specific edge attributes
    agDAG <- Rgraphviz::agopen(graph=dag, name="Ragraph", attrs=graphAttrs, nodeAttrs=node.attrs.default, edgeAttrs=edgeAttrs)

    ######################################################################################
    colNum <- 2
    rowNum <- 2
    tolNum <- colNum*rowNum
    
    ## a matrix object specifying the location of the next N figures on the output device. Each value in the matrix must be 0 or a positive integer
    layout_vec <- c(0,0,1,2)
    layout_matrix <- matrix(layout_vec, rowNum, colNum, byrow=T)
    
    ## relative heights and widths
    frac <- colorbar.fraction
    ## relative heights for each row
    row_first <- frac*10
    row_rest <- (10-row_first)
    layout_heights <- c(row_first, row_rest)
    ## relative widths for each column
    col_last <- frac*10
    col_rest <- (10-col_last)
    layout_widths <- c(col_rest, col_last)

    ######################################################################################

    if (newpage){
        grDevices::dev.new(width=width, height=height)
    }
    graphics::par(mfrow=c(rowNum,colNum), mar=margin)
    graphics::layout(layout_matrix, widths=layout_widths, heights=layout_heights)
    
    Rgraphviz::plot(agDAG, "dot")
    #slotNames(agDAG)
    #agDAG@AgNode
    #agDAG@AgEdge
    
    ## colorbar
    if(colorbar & !is.null(data)){
        
        ## empty plot
        plot(c(0,1),c(0,1),xlab="", ylab="", axes=F, type="n")
            
        palette.name <- visColormap(colormap=colormap)
        colors <- palette.name(ncolors)
        lab.scale <- length(colors)/(zlim[2]-zlim[1])
        for (i in 1:length(colors)) {
            yValue <- (i-1)/ncolors * 0.15
            hValue <- 1/ncolors * 0.15
            xValue <- 0.2
            wValue <- 0.2
        
            ## for rect
            xleft <- xValue
            ybottom <- yValue
            xright <- xValue+wValue
            ytop <- yValue+hValue
            graphics::rect(xleft,ybottom,xright,ytop, col=colors[i], border="transparent")
        
            if(i == 1 | i == 1+length(colors)/2){
                tx <- (i-1)/lab.scale + zlim[1]
                graphics::text(x=xright+0.2, y=ybottom, labels=tx, cex=1)
            }else if(i==length(colors)){
                tx <- i/lab.scale + zlim[1]
                graphics::text(x=xright+0.2, y=ytop, labels=tx, cex=1)
            }
        }
    }
    
    # restore original settings
    #suppressWarnings(graphics::par(opar))
    
    invisible(agDAG)
}