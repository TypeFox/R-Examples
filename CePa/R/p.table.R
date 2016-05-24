# p values of all pathways
p.table = function(x, adj.method = NA, cutoff = ifelse(adj.method == "none", 0.01, 0.05)) {

    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    n.pathway = length(x)
    p.value = matrix(0, nrow=length(x), ncol= length(x[[1]]))
    for(i in 1:length(x)) {
        p.value[i, ] = sapply(x[[i]], function(x) x$p.value)
    }

    rownames(p.value) = names(x)
    colnames(p.value) = names(x[[1]])
    
    if(!is.na(adj.method)) {
        p.value = apply(p.value, 2, p.adjust, adj.method)
        l = apply(p.value, 1, function(x) { sum(x <= cutoff) > 0})
        p.value = p.value[l, , drop=FALSE]
    }
    
    return(p.value)
}
