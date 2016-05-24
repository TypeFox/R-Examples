#' Function to visualise enrichment analysis outputs in the context of the ontology hierarchy 
#'
#' \code{visEnrichment} is supposed to visualise enrichment analysis outputs (represented as an 'Eoutput' object) in the context of the ontology hierarchy (direct acyclic graph; DAG). Only part of DAG induced by those nodes/terms specified in query nodes (and the mode defining the paths to the root of DAG) will be visualised. Nodes in query are framed in black (by default), and all nodes (in query plus induced) will be color-coded according to a given data.type ('zscore'; otherwise taking the form of 10-based negative logarithm for 'adjp' or 'pvalue'). If no nodes in query, the top 5 significant terms (in terms of adjusted p-value) will be used for visualisation
#'
#' @param e an object of S4 class \code{\link{Eoutput}}
#' @param nodes_query a verctor containing a list of nodes/terms in query. These nodes are used to produce a subgraph of the ontology DAG induced by them. If NULL, the top significant terms (in terms of p-value) will be determined by the next 'num_top_nodes'
#' @param num_top_nodes a numeric value specifying the number of the top significant terms (in terms of p-value) will be used. This parameter does not work if the previous 'nodes_query' has been specified
#' @param path.mode the mode of paths induced by nodes in query. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param data.type a character telling which data type for nodes in query is used to color-code nodes. It can be one of 'adjp' for adjusted p-values (by default), 'pvalue' for p-values and 'zscore' for z-scores. When 'adjp' or 'pvalue' is used, 10-based negative logarithm is taken. For the style of how to color-code, please see the next arguments: colormap, ncolors, zlim and colorbar
#' @param height a numeric value specifying the height of device
#' @param width a numeric value specifying the width of device
#' @param margin margins as units of length 4 or 1
#' @param colormap short name for the colormap. It can be one of "yr" (yellow-red colormap; by default), "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/data values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If data is null, it always sets to false
#' @param colorbar.fraction the relative fraction of colorbar block against the device size
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param layout.orientation the orientation of the DAG layout. It can be one of "left_right" for the left-right layout (viewed from the DAG root point; by default), "top_bottom" for the top-bottom layout, "bottom_top" for the bottom-top layout, and "right_left" for the right-left layout
#' @param node.info tells the ontology term information used to label nodes. It can be one of "both" for using both of Term ID and Name (the first 15 characters; by default), "none" for no node labeling, "term_id" for using Term ID, "term_name" for using Term Name (the first 15 characters), and "full_term_name" for using the full Term Name
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
#' @importFrom dnet dDAGinduce visDAG
#' @seealso \code{\link{dcEnrichment}}, \code{\link{dcRDataLoader}}, \code{\link{dcConverter}}
#' @include visEnrichment.r
#' @examples
#' \dontrun{
#' # 1) load SCOP.sf (as 'InfoDataFrame' object)
#' SCOP.sf <- dcRDataLoader('SCOP.sf')
#' # randomly select 20 domains
#' data <- sample(rowNames(SCOP.sf), 20)
#' 
#' # 2) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")
#' eoutput
#'
#' # 3) visualise the top 10 significant terms
#' # color-coded according to 10-based negative logarithm of p-values
#' visEnrichment(eoutput)
#' # color-coded according to zscore
#' visEnrichment(eoutput, data.type='zscore')
#'
#' # 4) visualise the top 5 significant terms in the ontology hierarchy
#' nodes_query <- names(sort(adjp(eoutput))[1:5])
#' visEnrichment(eoutput, nodes_query=nodes_query)
#' # change the frame color: highlight (framed in blue) nodes/terms in query
#' nodes.highlight <- rep("blue", length(nodes_query))
#' names(nodes.highlight) <- nodes_query
#' visEnrichment(eoutput, nodes_query=nodes_query, node.attrs=list(color=nodes.highlight))
#' }

visEnrichment <- function (e, nodes_query=NULL, num_top_nodes=5, path.mode=c("all_shortest_paths","shortest_paths","all_paths"), data.type=c("adjp","pvalue","zscore"), height=7, width=7, margin=rep(0.1,4), colormap=c("yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.1, newpage=T, layout.orientation=c("left_right","top_bottom","bottom_top","right_left"), node.info=c("both", "none", "term_id", "term_name", "full_term_name"), graph.node.attrs=NULL, graph.edge.attrs=NULL, node.attrs=NULL)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    path.mode <- match.arg(path.mode)
    data.type <- match.arg(data.type)
    layout.orientation <- match.arg(layout.orientation)
    node.info<- match.arg(node.info)
    
    ## check input Eoutput
    if (class(e) != "Eoutput"){
        stop("The function must apply to 'Eoutput' object.\n")
    }
    
    ######################################################################################
    ## load ontology (as an 'igraph' object)
    
    if(class(suppressWarnings(try(dcRDataLoader(paste('onto.', e@ontology, sep=''), verbose=F), T)))=="try-error"){
        g <- ''
        eval(parse(text=paste("g <- get(load('", e@ontology,"'))", sep="")))
    }else{
        g <- dcRDataLoader(paste('onto.', e@ontology, sep=''), verbose=F)
    }
    if(class(g)=="Onto"){
        g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
    }
    
    num_top_nodes <- as.integer(num_top_nodes)
    num_all <- length(pvalue(e))
    
    if(is.null(nodes_query)){
        if(num_top_nodes<1 & num_top_nodes>num_all){
            num_top_nodes <- min(5, num_all)
        }
        nodes_query <- names(sort(pvalue(e))[1:num_top_nodes])
    }else{
        ind <- match(nodes_query, V(g)$name)
        nodes_query <- nodes_query[!is.na(ind)]
        if(length(nodes_query)==0){
            warnings("Nodes/terms in your query are not found in the ontology!\nInstead, the top 5 significant terms (in terms of p-value) will be used.\n")
            if(num_top_nodes<1 & num_top_nodes>num_all){
                num_top_nodes <- min(5, num_all)
            }
            nodes_query <- names(sort(pvalue(e))[1:num_top_nodes])
        }
    }
    
    ## induce DAG only including nodes/terms in query
    subg <- dnet::dDAGinduce(g, nodes_query, path.mode=path.mode)
    
    ## prepare data used for color-coding
    msg <- ''
    if(data.type=="adjp"){
        data <- -1*log10(adjp(e))
        msg <- '-1*log10(adjusted p-values)'
    }else if(data.type=="pvalue"){
        data <- -1*log10(pvalue(e))
        msg <- '-1*log10(p-values)'
    }else{
        data <- zscore(e)
        msg <- 'z-scores'
    }
    
    if(sum(is.infinite(data))>0){
        data[is.infinite(data)] <- max(data[is.finite(data)])
    }
    
    ## customise highlighting of (framed) nodes/terms in query
    if(is.null(node.attrs)){
        ## color for nodes
        nodes.highlight <- rep("black", length(nodes_query))
        names(nodes.highlight) <- nodes_query
        ## color for text
        nodes.fontcolor <- rep("blue", length(nodes_query))
        names(nodes.fontcolor) <- nodes_query
        ## fontsize for text
        nodes.fontsize <- rep(18, length(nodes_query))
        names(nodes.fontsize) <- nodes_query
        
        #node.attrs <- list(color=nodes.highlight, fontcolor=nodes.fontcolor, fontsize=nodes.fontsize)
        node.attrs <- list(color=nodes.highlight, fontsize=nodes.fontsize)
        #node.attrs <- list(color=nodes.highlight)        
    }

    ## do visualisation
    agDAG <- dnet::visDAG(g=subg, data=data, height=height, width=width, margin=margin, colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=colorbar, colorbar.fraction=colorbar.fraction, newpage=newpage, layout.orientation=layout.orientation, node.info=node.info, graph.node.attrs=graph.node.attrs, graph.edge.attrs=graph.edge.attrs, node.attrs=node.attrs)
    
    message(sprintf("Ontology '%s' containing %d nodes/terms (including %d in query; also highlighted in frame) has been shown in your screen, with colorbar indicating %s", e@ontology, vcount(subg), length(nodes_query), msg), appendLF=T)
    
    invisible(agDAG)
}