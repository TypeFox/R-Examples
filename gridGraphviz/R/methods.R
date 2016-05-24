## Methods for laying out and plotting Ragraph objects

## The following methods return and Ragraph's dimensions in inches.
## Useful for determining appropriate sizes for output to device
graphWidth <- function (graph) {
  return(getX(upRight(boundBox(graph)))/72)
}
graphHeight <- function (graph) {
    return(getY(upRight(boundBox(graph)))/72)
}

## groomAttrs() takes the 'attrs' list passed to agopenTrue() and alters
## certain values to preferred defaults, overriding Rgraphviz defaults.
## Arguments:
##   attrs - list
## Values
##   attrs - list
groomAttrs <- function(attrs) {
    ## Rgraphviz sets the size of the graph to that of the current device
    ## or 7x7 inches by default. We would prefer that graphviz determines
    ## the graph size.
    ## Assume that if the user has a desired graph size this will be passed to
    ## agopenTrue() in the attrs list by:
    ##   agopenTrue(..., attrs=list(graph=list(size="XX,YY")))
    ## If attrs$graph$size in NULL assume user has not set a size
    ## and set this "" to allow graphviz to
    ## determine size of graph
    if (is.null(attrs$graph$size)) attrs$graph$size <- ""
    
    ## Rgraphviz sets node fixedsize="TRUE", width="0.75", height="0.5" on
    ## all nodes by default. We would prefer that none of these values be
    ## specified unless the user does so explicitly.
    ## Assume if user wants to control node sizes this will be passed to
    ## agopenTrue() in the attrs list by:
    ##   agopenTrue(..., attrs=list(node=list(fixedsize="VALUE", width="XX",
    ##                                        height="YY")))
    ## If attrs$node$fixedsize attrs$node$width attrs$node$height are NULL
    ## assume user has not set a size and set to "" to allow graphviz to
    ## determine values.
    if (is.null(attrs$node$fixedsize)) {
        attrs$node$fixedsize <- ""
    }
    if (is.null(attrs$node$width)) {
        attrs$node$width <- ""
    }
    if (is.null(attrs$node$height)) {
        attrs$node$height <- ""
    }

    ## return modified attrs list
    attrs
}

## groomEdgeAttrs() takes the edgeAttrs list passed to agopenTrue() and
## certain values to preferred defaults, overriding Rgraphviz defaults.
## Arguments:
##   graph - "graphNEL" object
##   edgeAttrs - list
## Values
##   edgeAttrs - list
groomEdgeAttrs <- function(graph, edgeAttrs) {
    ## Rgraphviz does not pass edge weights from a graph through to graphviz.
    ## This function will extract edge weights from graph and add a named
    ## vector of weights to the 'edgeAttrs' list
    weights <- unlist(edgeWeights(graph))
    names(weights) <- edgeNames(graph)
    edgeAttrs$weight <- weights

    ## return modified edgeAttrs list
    edgeAttrs
}

## agopenTrue() is a wrapper for Rgraphviz's agopen()
## it replaces certain Rgraphviz defaults with preferred defaults
## and attempts to pass through graph features which Rgraphviz does not
## NB: agopenTrue() does not accept a 'laidout' argument as graph should
##     always be laid out
agopenTrue <- function(graph, name, nodes, edges, kind = NULL,
                        layoutType = "dot", attrs = list(), nodeAttrs = list(), 
                        edgeAttrs = list(), subGList = list(),
                        edgeMode = edgemode(graph),
                        recipEdges = c("combined", "distinct")) {
    ## modify 'attrs' list to replace Rgraphviz defaults with preferred default
    ## values
    attrs <- groomAttrs(attrs)
    ## set edge weights to those specified in graph as Rgraphviz 
    edgeAttrs <- groomEdgeAttrs(graph, edgeAttrs)

    agopen(graph=graph, name=name, nodes=nodes, edges=edges, kind=kind, 
           layout=TRUE, layoutType=layoutType, attrs=attrs, 
           nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs, subGList=subGList,
           edgeMode=edgeMode, recipEdges=recipEdges)
}
