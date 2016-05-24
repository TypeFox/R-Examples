cepa.univariate.all = function(mat, label, pc, cen = default.centralities,
           cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
           nlevel = "tvalue_abs", plevel = "mean", iter = 1000) {
    
    nlevelFun = find.nlevelFun(nlevel)  # if nlevelFun is binary, the adj.method and cutoff have been hard coded in the functions
    plevelFun = find.plevelFun(plevel)
    
    # if cen argument is a function, the function should be quoted or substituted
    # because we need the function name
    for(ce in cen) {
        if(is.function(ce)) {
            stop("Functions cannot be used directly because we need the function name, use quote or substitute.\n")
        }
    }
    
    if(class(pc) != "pathway.catalogue") {
        stop("pc argument should be a pathway.catalogue object.")
    }
    
    # if binary transformation
    # function need a whole expression matrix to calculate adjust p-values
    cat("  Calculate gene level values.\n")
    
    n.pathway = length(pc$pathList)
    
    # initialize pathway.result
    pathway.name = names(pc$pathList)    
    pathway.result = list()
    length(pathway.result) = n.pathway
    # pathway.result is like a two layer list
    pathway.result = lapply(pathway.result, function(x) {
                              y = list()
                              length(y) = length(cen.name)
                              names(y) = cen.name
                              return(y)
                            })
    names(pathway.result) = pathway.name
    
    cat("  Calculate pathway score...\n")
    
    for(i in 1:length(pc$pathList)) {
        
        cat("    ", i, "/", length(pc$pathList), ", ", pathway.name[i], "...\n", sep="")
        
        path = pc$pathList[[i]]
        inter = pc$interactionList[pc$interactionList[, 1] %in% path, 2:3]
        
        pathway = generate.pathway(as.matrix(inter))
        
        
        
        
        cat("      Calculate node level value and permutate sample labels...\n")
        # nodes in the pathway
        node = pathway.nodes(pathway)
        
         # only the mapping in the pathway
        mapping = pc$mapping[pc$mapping[, 1] %in% node, ]
    
        # gene expression just in the pathway
        l = rownames(mat) %in% mapping[, 2]
        cat("      ", sum(l), " genes measured in the pathway...\n", sep = "")
        # if no genes are measure in the pathway
        if(sum(l) == 0) {
            node.level.from.expr = rep(0, length(node))
            node.level.t.value = rep(0, length(node))
            r.node.level.from.expr = matrix(0, nrow = length(node), ncol = iter)
        } else {
            mat.gene = mat[l, ,drop = FALSE]
            mat.gene = t(apply(mat.gene, 1, scale))  # really need to be scaled?
            
            # generate node expression from its member genes
            # if there is no member gene, the expression is zero
            mat.node = matrix(0, nrow=length(node), ncol = dim(mat)[2])
            rownames(mat.node) = node
            for(k in 1:length(node)) {
                l = mapping[, 1] == node[k]
                gene.in.node = unique(mapping[l, 2])
                gene.in.node.in.mat = gene.in.node[gene.in.node %in% rownames(mat.gene)]
                if(length(gene.in.node.in.mat) == 1) {
                    mat.node[k, ] = mat.gene[gene.in.node.in.mat, ]
                } else if(length(gene.in.node.in.mat > 1)) {  # use the biggeset component of member gene expression matrix
                    mm = t(mat.gene[gene.in.node.in.mat, ])
                    pcar = prcomp(mm)
                    mat.node[k, ] = predict(pcar, mm)[, 1]
                }
            }
            
            node.level.from.expr = apply(mat.node, 1, function(x) nlevelFun(x[.treatment(label)], x[.control(label)]))
            node.level.t.value = apply(mat.node, 1, function(x) nodeLevelFun.tvalue(x[.treatment(label)], x[.control(label)]))
            
            node.level.from.expr[is.na(node.level.from.expr)] = 0
            node.level.t.value[is.na(node.level.t.value)] = 0
    
            r.node.level.from.expr = matrix(0, nrow = length(node), ncol = iter)
            for(k in 1:iter) {
                r.label = .permutate(label)
                r.node.level.from.expr[, k] = apply(mat.node, 1, function(x) nlevelFun(x[.treatment(r.label)], x[.control(r.label)]))
                r.node.level.from.expr[is.na(r.node.level.from.expr[, k]), k] = 0
            }
        }
        
        j = 0
        for(ce in cen) {
            j = j + 1
            pathway.result[[i]][[j]] = cepa.univariate(mat = mat, label = label, pc = pc, pathway = pathway,
                                                    cen = ce, iter = iter, nlevel = nlevel, plevel = plevel,
                                                    node.level.from.expr = node.level.from.expr, node.level.t.value = node.level.t.value,
                                                    r.node.level.from.expr = r.node.level.from.expr)
            
            cat("      - ", ce, ": ", round(pathway.result[[i]][[j]]$p.value, 3), "\n", sep = "")
        }
    }

    class(pathway.result) = "cepa.all"
    return(pathway.result)

}

# if gene.level and r.gene.level 
# mat, label, ...
cepa.univariate = function(mat, label, pc, pathway = NULL, id = NULL, cen = "equal.weight",
                   cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
                   iter = 1000, nlevel = "tvalue_abs", plevel = "mean",
                   node.level.from.expr = NULL, node.level.t.value = NULL,
                   r.node.level.from.expr = NULL) {
    
    # in this version, we do not allow user-defined nlevel, nlevel and plevel functions
    nlevelFun = find.nlevelFun(nlevel)  # if nlevelFun is binary, the adj.method and cutoff have been hard coded in the functions
    plevelFun = find.plevelFun(plevel)
    
    
    # you can specify a pathway object or a pathway id
    if(! is.null(pathway)) {   # a pathway is specified
        if(is.matrix(pathway) || is.data.frame(pathway)) {   # only edge list
            if(length(dim(pathway)) != 2 || dim(pathway)[2] != 2) {
                stop("if pathway is a matrix or data frame, it should be 2 dimension and the number of columns is 2.\n")
            }
            pathway = generate.pathway(pathway)    # generate an igraph object
        } else if(class(pathway) != "igraph") {    # it should be an igraph object
            stop("Since pathway is not formatted as edge list, it should be an igraph object.")
        }
    } else if(! is.null(id)) {  # if the pathway is not specified, but the pathway ID is available
        # get interactions in the pathway
        path = pc$pathList[[id]]
        inter = pc$interactionList[pc$interactionList[, 1] %in% path, 2:3]
        # generate graph from edge list
        pathway = generate.pathway(inter)
    } else {  # one of pathway and id should be set
        stop("You should specify pathway argument or id argument.")
    }

    if(iter < 100) {
        stop("Iterations should not be smaller than 100.\n")
    }
    
    # single centrality!
    if(length(cen) != 1) {
        stop("Length of cen must be equal to 1.\n") 
    }
    
    weight = centrality(pathway, cen)

    add = 0
    # if there are none-zero weight values
    if(sum(weight == 0) != length(weight)) {
        add = min(weight[weight > 0])/100
    }
    weight = weight + ifelse(sum(weight == weight[1]) == length(weight), 0, add)
    
    # nodes in the pathway
    node = pathway.nodes(pathway)
    
     # only the mapping in the pathway
    mapping = pc$mapping[pc$mapping[, 1] %in% node, ]
    
    # get node names formatted with genes
    node.name = node
    member = character(0)
    for(i in 1:length(node)) {
        # genes that exsit in the node
        l = mapping[, 1] == node[i]
        
        # if find nodes with genes mapped
        if(sum(l)) {
            member = sort(unique(mapping[l, 2]))
            node.name[i] = paste(member, collapse = "\n")
        }
    }
    
    if(is.null(node.level.from.expr)) {
        # gene expression just in the pathway
        l = rownames(mat) %in% mapping[, 2]
        if(sum(l) == 0) {
            node.level.from.expr = rep(0, length(node))
            node.level.t.value = rep(0, length(node))
        } else {
            mat.gene = mat[l, ,drop=FALSE]
            mat.gene = t(apply(mat.gene, 1, scale))  # really need to be scaled?
            
            # generate node expression from its member genes
            # if there is no member gene, the expression is zero
            mat.node = matrix(0, nrow=length(node), ncol = dim(mat)[2])
            rownames(mat.node) = node
            for(i in 1:length(node)) {
                l = mapping[, 1] == node[i]
                gene.in.node = unique(mapping[l, 2])
                gene.in.node.in.mat = gene.in.node[gene.in.node %in% rownames(mat.gene)]
                if(length(gene.in.node.in.mat) == 1) {
                    mat.node[i, ] = mat.gene[gene.in.node.in.mat, ]
                } else if(length(gene.in.node.in.mat > 1)) {  # use the biggeset component of member gene expression matrix
                    mm = t(mat.gene[gene.in.node.in.mat, ])
                    pcar = prcomp(mm)
                    mat.node[i, ] = predict(pcar, mm)[, 1]
                }
            }
            
            node.level.from.expr = apply(mat.node, 1, function(x) nlevelFun(x[.treatment(label)], x[.control(label)]))
            node.level.t.value = apply(mat.node, 1, function(x) nodeLevelFun.tvalue(x[.treatment(label)], x[.control(label)]))
        }
    }
    
    node.level.from.expr[is.na(node.level.from.expr)] = 0
    node.level.t.value[is.na(node.level.t.value)] = 0
    
    node.level = weight * node.level.from.expr
    ds = quantile(node.level, c(1, 0.75, 0.5, 0))
    names(ds) = c("max", "q75", "median", "min")
    
    s = plevelFun(node.level)    # calculate nodes
    
    # now the permutation
    s.random = numeric(iter)
    ds.random = matrix(0, iter, 4)   # descriptive statistic of the node
    colnames(ds.random) = c("max", "q75", "median", "min")
    if(sum(rownames(mat) %in% mapping[, 2]) == 0) {
        s.random = rep(0, iter)
    } else {
        for(i in 1:iter) {
            
            if(is.null(r.node.level.from.expr)) {
                r.label = .permutate(label)
                r.node.level.from.expr.current = apply(mat.node, 1, function(x) nlevelFun(x[.treatment(r.label)], x[.control(r.label)]))
            } else {
                r.node.level.from.expr.current = r.node.level.from.expr[, i]
            }
            
            r.node.level.from.expr.current[is.na(r.node.level.from.expr.current)] = 0
    
            r.node.level = weight * r.node.level.from.expr.current
            s.random[i] = plevelFun(r.node.level)

            ds.random[i, ] = quantile(r.node.level, c(1, 0.75, 0.5, 0))
        }
    }
    
    p.value = (sum(s.random >= s) + 1) / (iter + 1)
    
    res = list("score" = s,                                  # pathway score
               "score.distribution" = ds,                    # distribution of node value in the pathway
               "score.random" = s.random,                    # simulated pathway scores
               "score.distribution.random" = ds.random,      # distribution of node value in the pathway in each simulation
               "p.value" = p.value,                          # p value
               "centrality" = cen.name,                      # centrality name
               "weight" = weight,                            # weight for each node
               "node.level.t.value" = node.level.t.value,    # value for each node, exclude the centrality part
               "node.level" = node.level,                    # value for each node, exclude the centrality part
               "node.name" = node.name,                      # node names
               "pathway" = pathway,                          # pathway in igraph format
               "framework" = "gsa.univariate")                          
               
    class(res) = "cepa"
    
    return(invisible(res))

}

# calculate node.level.only.from.expr by expression of member genes
node.score = function(gene.level = NULL, gene.in.node = NULL, nlevelFun = NULL) {

    node.level.from.expr = numeric(length(gene.in.node))
    for(i in 1:length(gene.in.node)) {
        node.level.from.expr[i] = ifelse(length(gene.in.node[[i]]), nlevelFun(gene.level[names(gene.level) %in% gene.in.node[[i]]]), 0)
    }
    return(node.level.from.expr)
}

nodeLevelFun.tvalue = function(x1, x2, ...) {
    n1 = length(x1)
    n2 = length(x2)
    v1 = var(x1)
    v2 = var(x2)
    ifelse(v1 + v2 == 0, 0, (mean(x1) - mean(x2)) / sqrt(v1/n1 + v2/n2))
}
nodeLevelFun.tvalue_sq = function(x1, x2, ...) {
    nodeLevelFun.tvalue(x1, x2)^2
}
nodeLevelFun.tvalue_abs = function(x1, x2, ...) {
    abs(nodeLevelFun.tvalue(x1, x2))
}


pathwayLevelFun.max = function(x) {
    max(x, na.rm = TRUE)
}
pathwayLevelFun.min = function(x) {
    min(x, na.rm = TRUE)
}
pathwayLevelFun.median = function(x) {
    median(x, na.rm = TRUE)
}
pathwayLevelFun.sum = function(x) {
    sum(x, na.rm = TRUE)
}
pathwayLevelFun.mean = function(x) {
    mean(x, na.rm = TRUE)
}
pathwayLevelFun.rank = function(x) {
    wilcox.test(x, exact = FALSE)$statistic
}


# arguments are different according to different nlevel functions
find.nlevelFun = function(method) {
    if(is.character(method)) {
        f = switch(method,
                   tvalue     = nodeLevelFun.tvalue,
                   tvalue_sq  = nodeLevelFun.tvalue_sq,
                   tvalue_abs = nodeLevelFun.tvalue_abs,
                   stop("Default node-level functions in string format are only tvalue, tvalue_sq and tvalue_abs."))
    } else if ( is.function(method) ) {
        if(length(as.list(args(method))) != 3) {
            stop("Self-defined node-level function only can have two arguments.")
        }
        f = method
    }
    return(f)
}

find.plevelFun = function(method) {
    if(is.character(method)) {
        f = switch(method,
               max    = pathwayLevelFun.max,
               min    = pathwayLevelFun.min,
               median = pathwayLevelFun.median,
               sum    = pathwayLevelFun.sum,
               mean   = pathwayLevelFun.mean,
               rank   = pathwayLevelFun.rank,
               stop("Wrong pathway level method."))
    } else if ( is.function(method) ) {
        if(length(as.list(args(method))) != 2) {
            stop("Self-defined pathway-level function only can have one argument.")
        }
        f = method
    }
    return(f)
}
