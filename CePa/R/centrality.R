# caculate the centrality of nodes in pathways
# graph: igrpah object
# method: string or function or name
centrality = function(graph, method="equal.weight") {
    
    if(length(method) > 1) {
        stop("Length of method must be equal to 1.\n") 
    }
    
    if(is.function(method)) {
        return(method(graph))
    } else if(mode(method) == "name") {
        method = eval(method)
        return(method(graph))
    } else if(method == "equal.weight") {
        return(rep(1, vcount(graph)))
    } else if(method == "in.degree") {
        return(degree(graph, mode="in"))
    } else if(method == "out.degree") {
        return(degree(graph, mode="out"))
    } else if(method == "degree") {
        return(degree(graph))
    } else if(method == "betweenness") {
        return(betweenness(graph))
    } else if(method == "in.reach") {
        return(reach(graph, mode="in"))
    } else if(method == "out.reach") {
        return(reach(graph, mode="out"))
    } else if(method == "reach") {
        return(reach(graph))
    } else if(method == "in.spread") {
        return(spread(graph, mode="in"))
    } else if(method == "out.spread") {
        return(spread(graph, mode="out"))
    } else if(method == "spread") {
        return(spread(graph))
    } else {
        stop("Wrong centrality method.\n")
    }
}


spread = function(graph, mode=c("all", "in", "out"),
                  weights=E(graph)$weight, f=function(x)1/x) {
    mode = mode[1]
    
    sp = shortest.paths(graph, mode=mode, weights=weights)
    s = apply(sp, 1, function(x) {
            return(sum(f(x[x>0])))
        })
    return(s)
}

# the longest shortest path from v to other nodes
reach = function(graph, weights=E(graph)$weight, mode=c("all", "in", "out")) {
    mode = mode[1]
    
    sp = shortest.paths(graph, weights=weights, mode=mode)
    s = apply(sp, 1, function(x) {
            if(all(x == Inf)) {
                return(0)
            }
            else {
                return(max(x[x != Inf]))
            }
        })
    return(s)
}

radiality = function(graph, mode=c("all", "in", "out")) {
    mode = mode[1]
    
    sp = shortest.paths(graph, mode=mode, weights=NA)
    n = vcount(graph)
    diam = diameter(graph, directed = ifelse(mode == "all", FALSE, TRUE))
    s = apply(sp, 1, function(x) {
            if(all(x == Inf)) {
                return(0)
            }
            else {
                return(sum(diam+1-x[x != Inf])/(n-1))
            }
        })
    return(s)
}

