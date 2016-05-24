
getObjectList <- function(grobs=TRUE, viewports=TRUE) {
	object <- grid.ls(grobs=grobs, viewports=viewports, print=FALSE)
    	object <- lapply(object, "[", !(object$type == "vpPopListing" |
                                        object$type == "vpUpListing" |
                                        object$type == "vpNameListing" |
                                        object$type == "vpNameTreeListing"))
	object
}



appendDupNames <- function(names) {
    for (i in names) {
        dups <- names == i
        ndups <- sum(dups)
        if (ndups > 1) {
            names[dups] <- paste(names[dups],"#", 
                                 1:ndups, sep="")
        }
    }
    names
}


generateNodeNames <- function(objs) {
    # Add prefix
    names <- paste(ifelse(objs$type == "vpListing" |
                          objs$type == "vpTreeListing",
                          "vp", "g"),
                   objs$name, sep=":")
    # Add suffixes
    names <- appendDupNames(names)
    names
}


renameGrobPaths <- function(path, objs) {
    for (i in 1:length(objs$name)) {
        # If I am a gTree
        if (objs$type[i] == "gTreeListing") {
            if (path[i] == "") {
                searchPattern <- paste("^", objs$name[i], sep="")
                replacement <- objs$nodes[i]
                path <- sub(searchPattern, replacement, path)
            } else {
                searchPattern <- paste(paste("^", path[i], sep=""),
                                       objs$name[i], sep="::")
                replacement <- paste(path[i], objs$nodes[i], sep="::")
                path <- sub(searchPattern, replacement, path)
            }
        }
    }
    path
}

renameVPPaths <- function(path, objs) {
    for (i in 2:length(objs$name)) {
        if (objs$type[i-1] == "vpListing" |
            objs$type[i-1] == "vpTreeListing") {
            if (objs$name[i-1] == "ROOT") {
                path = sub(objs$name[i-1], objs$nodes[i-1], path)
            } else {
                path = sub(paste0(path[i-1], "::", objs$name[i-1], "($|::)"),
                  paste0(path[i-1], "::", objs$nodes[i-1], "\\1"), path)
            }
        }
    }
    path
}

edgeEnds <- function(node, nodes, paths) {
    ends <- nodes[grep(paste(node, "$", sep=""), paths)]
    if (length(ends) > 0) {
        list(edges=ends)
    } else {
        list()
    }
}

generateEdgeList <- function(objs) { 
    
    objs$vpPath = renameVPPaths(objs$vpPath, objs)
    objs$gPath = renameGrobPaths(objs$gPath, objs)
    
    edgeVPList <- lapply(objs$nodes, edgeEnds, objs$nodes, objs$vpPath)
    edgeGrobList <- lapply(objs$nodes, edgeEnds, objs$nodes, objs$gPath)
    names(edgeVPList) <- objs$nodes
    names(edgeGrobList) <- objs$nodes
    
    edgeList <- mapply(c, edgeVPList, edgeGrobList)
    edgeList
}



makeGraph <- function(nodes, edges, mode) {
    new("graphNEL",
        nodes=nodes,
        edgeL=edges,
        edgemode=mode)
}






makeRaGraph <- function(graph, objs,
                        gattr, vpattr,
                        g2gattr, vp2vpattr, g2vpattr, vp2gattr,
                        split) {
##########
#node attributes
##########
    nodeAttrs <- list()
    for (i in unique(c(names(gattr), names(vpattr)))) {
        if (is.null(gattr[[i]]))
            gattr[[i]] <- ""
        if (is.null(vpattr[[i]]))
            vpattr[[i]] <- ""
        nodeAttrs[[i]] <- ifelse(grepl("grobListing|gTreeListing", objs$type),
                                 gattr[[i]], vpattr[[i]])
        names(nodeAttrs[[i]]) <- objs$nodes
    }
    if (is.null(nodeAttrs$label)) {
        if (split) {
            nodeAttrs$label <- sapply(strsplit(objs$name, "[.-]"),
                                      function(bits) {
                                          if (length(bits) == 1)
                                              bits
                                          else {
                                              label <- bits[1]
                                              length <- nchar(bits[1])
                                              sep <- "."
                                              for (i in 2:length(bits)) {
                                                  nc <- nchar(bits[i])
                                                  if (length + nc > 7) {
                                                      sep <- ".\\\n"
                                                      length <- nc
                                                  } else {
                                                      sep="."
                                                      length <- length + nc
                                                  }
                                                  label <- paste(label, bits[i],
                                                                 sep=sep)
                                              }
                                              label
                                          }
                                      })
        } else {
            nodeAttrs$label <- objs$name
        }
        names(nodeAttrs$label) <- objs$nodes
    }
    
##########
#edge attributes
##########
        
    begType <- ifelse(grepl("grobListing|gTreeListing", objs$type),
                      "g" ,"vp")
    endType <- lapply(edges(graph),
                      function(edges) {
                          begType[match(edges, objs$nodes)]
                      })
    edgeTypes <- unlist(mapply(function(begin, ends) {
                                   if (length(ends) == 0)
                                       NULL
                                   else
                                       paste(begin, ends, sep="~")
                               }, begType, endType))
    objs$edgeNames <- edgeNames(graph)
    edgeAttrs <- list()
    for (i in unique(c(names(g2gattr), names(vp2vpattr),
                       names(g2vpattr), names(vp2gattr)))) {
        if (is.null(g2gattr[[i]]))
            g2gattr[[i]] <- ""
        if (is.null(vp2vpattr[[i]]))
            vp2vpattr[[i]] <- ""
        if (is.null(g2vpattr[[i]]))
            g2vpattr[[i]] <- ""
        if (is.null(vp2vpattr[[i]]))
            vp2vpattr[[i]] <- ""
        edgeAttrs[[i]] <- ifelse(edgeTypes == "g~g",
                                 g2gattr[[i]],
                                 ifelse(edgeTypes == "vp~vp",
                                        vp2vpattr[[i]],
                                        ifelse(edgeTypes == "g~vp",
                                               g2vpattr[[i]],
                                               vp2gattr[[i]])))
        names(edgeAttrs[[i]]) <- objs$edgeNames
    }

    agopenTrue(graph, name="", nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
}


gridTree <- function(grobNodeAttrs=list(shape="circle", fillcolor="black",
                       fontcolor="white"),
                     vpNodeAttrs=list(shape="box", fillcolor="grey90",
                       fontcolor="black"),
                     grob2grobAttrs=list(color="black", lty="solid", lwd=1),
                     vp2vpAttrs=list(color="black", lty="solid", lwd=1),
                     grob2vpAttrs=list(color="black", lty="dotted", lwd=1),
                     vp2grobAttrs=list(color="grey", lty="solid", lwd=2),
                     split=TRUE, grid=TRUE, 
                     grobs=TRUE, viewports=TRUE, draw=TRUE) {

    objs <- getObjectList(grobs, viewports)
    objs$nodes <- generateNodeNames(objs)
    objs$edgeList <- generateEdgeList(objs)
    grid.ls.GNEL <- makeGraph(objs$nodes, objs$edgeList, "directed")	
    grid.ls.Ragraph <- makeRaGraph(grid.ls.GNEL, objs,
                                   grobNodeAttrs, vpNodeAttrs,
                                   grob2grobAttrs, vp2vpAttrs,
                                   grob2vpAttrs, vp2grobAttrs,
                                   split)
    if (draw) {
        if (grid) {
            grid.graph(grid.ls.Ragraph, newpage=TRUE)
        } else {
            plot(grid.ls.Ragraph)
        }
    }
    invisible(grid.ls.Ragraph)
}
