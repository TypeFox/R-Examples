## "smoothing2" uses "smooth2", and is used in "RNMF".
## "smoothing2" takes in a "graph" (data frame, or vector, or matrix). It returns either a matrix or a data frame of the smoothed data frame.
## NOTE: Examples of "smoothing2" can be found in "smoothing2.r"
smoothing2 = function(graph, outs, p = NULL, n = NULL, frame = FALSE){
    if(is.vector(graph)){
        if(missing(p) & missing(n)) stop("Missing p and n!")
        graph = data.frame(x = rep(1:p,n), y = rep(1:n,each = p), value = graph)
    }else if(is.matrix(graph)){
        p = nrow(graph); n = ncol(graph)
        graph = data.frame(x = rep(1:p,n), y = rep(1:n,each = p), value = as.vector(graph))
    }
    graph$outs = FALSE
    graph$outs[outs] = TRUE
    graph[outs,"value"] = unlist(lapply(outs, smooth2, graph, p, n))
    if(frame == TRUE){
        return(graph)
    }else{
        return(matrix(graph$value, nrow = p))
    }
}

