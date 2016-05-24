
# centraliy extension on Over-representative analysis method
# on a list of pathways and under several centralities
cepa.ora.all = function(dif, pc, bk = NULL, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    iter = 1000) {
    
    # if no background gene list is specified, use whole human genome
    if(is.null(bk)) {
        dir = system.file(package = "CePa")
        bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        cat("Background gene list is not specified, use whole human genome instead.\n")
    }
    
    dif = dif[dif %in% bk]
    
    if(length(cen) < 1) {
        stop("cen argument must be specified.\n")
    }
    
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
    
    # generate simulated differential gene list
    # a list in which each item is a vector of differential genes
    #dif.random = list()
    #length(dif.random) = iter
    #for(i in 1:iter) {
    #    dif.random[[i]] = sample(bk, length(dif), replace = FALSE)
    #}
    
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
    
    cat("  Calculate pathway scores...\n")
    
    for(i in 1:n.pathway) {
        
        cat("    ", i, "/", n.pathway, ", ", pathway.name[i], "...\n", sep="")
        
        # interaction ID in this pathway
        path = pc$pathList[[i]]
        # interactions in this pathway
        inter = pc$interactionList[pc$interactionList[, 1] %in% path, 2:3]
        # generate graph from edge list
        pathway = generate.pathway(as.matrix(inter))
        
        j = 0
        # to this pathway, applying various centralities
        for(ce in cen) {
            j = j + 1
            pathway.result[[i]][[j]] = cepa.ora(dif = dif, bk = bk, pathway = pathway, pc = pc, cen = ce, iter = iter)
            cat("      - ", ce, ": ", round(pathway.result[[i]][[j]]$p.value, 3), "\n", sep = "")
        }
    }

    class(pathway.result) = "cepa.all"
    return(pathway.result)
}

# apply centrality-extension of Over-representative analysis method
# on a single pathway and a single centrality
cepa.ora = function(dif, pc, bk = NULL, pathway = NULL, id = NULL, cen = "equal.weight",
                cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
                iter = 1000) {
    
    # if no background gene list is specified, use whole human genome
    if(is.null(bk)) {
        dir = system.file(package = "CePa")
        bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        cat("Background gene list is not specified, use whole human genome instead.\n")
    }
    
    # some checking of the arguments
    if(length(dif) > length(bk)) {
        stop("Length of differential genes should not be larger than the length of background genes.\n")
    }
    if(sum(dif %in% bk) != length(dif)) {
        stop("Differential genes must be all in background list.\n")
    }
    
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
    
    # is our algorithm, only non-negative centrality is allowed
    if(any(weight < 0)) {
        stop("Weight should not be negative.")
    }

    add = 0
    # if there are none-zero weight values
    if(any(weight > 0)) {
        add = min(weight[weight > 0])/100
    }
    weight = weight + add
    
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
            # mark the diff genes
            member[member %in% dif] = paste("[", member[member %in% dif], "]", sep="")
            node.name[i] = paste(member, collapse = "\n")
        }
    }
    
    # map dif genes to node id
    dif.node = unique(mapping[mapping[, 2] %in% dif, 1])
    
    is.dif.node = as.numeric(node %in% dif.node)
    node.level = is.dif.node * weight
    s = sum(node.level)
    
    # distribution of node level value (combined with weight)
    if(sum(is.dif.node > 0) == 0) {   # if there is no differential genes
        ds = c(0, 0, 0, 0)
    } else {
        ds = quantile(node.level, c(1, 0.75, 0.5, 0))
    }
    names(ds) = c("max", "q75", "median", "min")
    
    # sampling
    p.dif = length(dif) / length(bk)  # probability to be a differential gene
    s.random = numeric(iter)
    ds.random = matrix(0, iter, 4)   # descriptive statistic of the node
    colnames(ds.random) = c("max", "q75", "median", "min")
    
    # genes in the pathway
    gene = unique(mapping[mapping[, 1] %in% node, 2])
    
    for(i in 1:iter) {
        # simulated random genes
        dif.gene.random = gene[as.logical(rbinom(length(gene), 1, p.dif))]
        
        # then map to node id
        dif.node.random = unique(mapping[mapping[, 2] %in% dif.gene.random, 1])
        # find which node is differentially affected
        is.dif.node.random = as.numeric(node %in% dif.node.random)
        # calculate the score
        node.level.random = is.dif.node.random * weight
        s.random[i] = sum(node.level.random)
        if(sum(is.dif.node.random > 0) == 0) {
            ds.random[i, ] = c(0, 0, 0, 0)
        }
        else {
            ds.random[i, ] = quantile(node.level.random, c(1, 0.75, 0.5, 0))
        }
    }
    
    p.value = (sum(s.random >= s) + 1) / (iter + 1)
    
    dif.gene = intersect(dif, gene)
    n.dif.node = length(dif.node)
    n.dif.gene = length(dif.gene)
    n.node = length(node)
    n.gene = length(gene)
    
    count = c(n.dif.node, n.node, n.dif.gene, n.gene)
    names(count) = c("n.dif.node", "n.node", "n.dif.gene", "n.gene")
    
    
    res = list("score" = s,                                  # pathway score
               "score.distribution" = ds,                    # distribution of node value in the pathway
               "score.random" = s.random,                    # simulated pathway scores
               "score.distribution.random" = ds.random,      # distribution of node value in the pathway in each simulation
               "p.value" = p.value,                          # p value
               "centrality" = cen.name,                      # centrality name
               "weight" = weight,                            # weight for each node
               "node.level.t.value" = as.integer(is.dif.node),    # value for each node, exclude the centrality part
               "node.level" = node.level,                    # value for each node, exclude the centrality part
               "node.name" = node.name,                      # node names
               "pathway" = pathway,                          # pathway in igraph format
               "count" = count,
               "framework" = "ora")                          
    
    class(res) = "cepa"
    
    return(invisible(res))

}

