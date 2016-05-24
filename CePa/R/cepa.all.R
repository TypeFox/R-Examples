
# wrapper of cepa.ora.all, cepa.univariate.all
# choose corresponding fucntions according to the arguments
# arguments for cepa.ora.all:
#   dif, bk, pc, cen, cen.name, iter
# arguments for cepa.univariate.all:
#   mat, label, pc, cen, cen.name, nlevel, plevel, iter
cepa.all = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000 ) {
    
    # for those who are lazy to specify argument names
    # if the first argument is a vector, then it is ora method
    if(is.vector(dif)) {
        # do nothing
        if(is.null(bk)) {
            dir = system.file(package = "CePa")
            bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        }
    } else if(is.matrix(dif)) {
        mat = dif
        dif = NULL
    } else if(is.data.frame(dif)) {
        mat = as.matrix(dif)
        dif = NULL
    }
    
    if(! is.null(dif)) {     # if dif is specified
        res = cepa.ora.all(dif = dif, bk = bk, pc = pc, cen = cen, cen.name = cen.name, iter = iter)
    } else {
        res = cepa.univariate.all(mat = mat, label = label, pc = pc, cen = cen, cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter)
    }
    
    return(res)
}


# wrapper of cepa.ora, cepa.univariate, cepa.multivariate
# choose corresponding fucntions according to the arguments
# arguments for cepa.ora:
#   dif, bk, pc, pathway, id, cen, cen.name, iter
# arguments for cepa.univariate:
#   mat, label, pc, pathway, id, cen, cen.name, nlevel, nlevel, plevel, iter, gene.level, r.gene.level
cepa = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, pathway = NULL, id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000) {
    
    # if the first argument is a vector, then it is ora method
    if(is.vector(dif)) {
        if(is.null(bk)) {
            dir = system.file(package = "CePa")
            bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        }
    } else if(is.matrix(dif)) {
        mat = dif
        dif = NULL
    } else if(is.data.frame(dif)) {
        mat = as.matrix(dif)
        dif = NULL
    }
    
    if(! is.null(dif)) {     # if dif is specified
        cat("  Applying Over-representative analysis.\n")
        res = cepa.ora(dif = dif, bk = bk, pc = pc, pathway = pathway, id = id, cen = cen, cen.name = cen.name, iter = iter)
    } else {
        cat("  Applying Gene-set analysis (univariate procedure).\n")
        res = cepa.univariate(mat = mat, label = label, pathway = pathway, pc = pc, id = id, cen = cen,
                   cen.name = cen.name, iter = iter, nlevel = nlevel, plevel = plevel)
    }
    
    return(res)
    
}

# A pathway data contains a list of pathways in which pathways are represented as a list of interactions.
# An interactions ID is assigned to each interaction, so each pathway contains a list of interaction IDs.
# Also, a list describing interactions should be provided. The interaction data contains three columns where
# the first column is the interaction ID, the second column is the input node ID and the third column is
# the output node ID. At last, a mapping list from node ID to gene ID is provided.
# Since pathway data is not changed in the analysis, we integrate all the three kinds of pathway data into 
# a whole pathway.catalogue object.
# There are some simple pre-process of the pathway data, such as retrict pathways from the number of nodes
# and the number of genes.
set.pathway.catalogue = function(pathList, interactionList, mapping,
    min.node = 5, max.node = 500, min.gene = min.node, max.gene = max.node, ...) {
    
    # pathway should be list
    if(!is.list(pathList)) {
        stop("pathList should be a list.\n")
    }
    
    # pathway should have name
    if(is.null(names(pathList))) {
        stop("pathList should have names.\n")
    }
    
    # interactionList is a matrix
    if(length(dim(interactionList)) != 2) {
        stop("interactionList should be two dimension matrix.\n")
    }
    
    # interactionList contains three columns
    if(dim(interactionList)[2] != 3) {
        stop("interactinList should contain 3 columns.\n")
    }
    
    l = sapply(pathList, function(x) {
                   it = interactionList[interactionList[, 1] %in% x, 2:3]   # interaction list for the pathway
                   node = unique(c(it[, 1], it[, 2]))                       # nodes in the pathway
                   l.node = length(node)                                    # number of nodes
                   gene = unique(mapping[mapping[,1] %in% node, 2])         # genes in the pathway
                   l.gene = length(gene)                                    # number of genes
                   return(c(l.node, l.gene))
               })
    pathList = pathList[l[1, ] >= min.node & l[1, ] <= max.node
                        & l[2, ] >= min.gene & l[2, ] <= max.gene]
    
    if(length(pathList) == 0) {
        warning("Your pathway catalogue is empty!")
    }
    
    pc = list(pathList = pathList, 
              interactionList = interactionList,
              mapping = mapping,
              ...)        # name of the catalogue, such as KEGG, PID ...
    class(pc) = "pathway.catalogue"
    return(pc)
}

# simply print the summary of pathway.catalogue object.
# simply print the summary of pathway.catalogue object.
print.pathway.catalogue = function(x, ...) {
    cat("\n  The catalogue contains", length(x$pathList), "pathways.\n\n")
}

# plot the distribution of resident of genes and number genes for nodes
# scatter plot of number of genes and number of nodes in the catalogue
plot.pathway.catalogue = function(x, ...) {

    r1 = numeric(0)
    r2 = numeric(0)
    i = 0
    
    # first column is the number of genes in each pathway
    # second column is the number of nodes in each pathway
    gn = matrix(0, nrow=length(x$pathList), ncol=2)
    colnames(gn) = c("node", "gene")
    for(pa in x$pathList) {
        # interaction ID list
        int = pa
        # node IDs in the pathway
        l  = x$interactionList[, 1] %in% int
        node = unique(c(x$interactionList[l, 2], x$interactionList[l, 3]))
        # mapping in this pathway
        m = x$mapping[x$mapping[, 1] %in% node, ]
        # number of nodes contain k genes
        ta = table(table(m[, 1]))
        # how many genes that a node contains
        tname = names(ta)
        tname = as.integer(tname)
        r1 = c(r1, rep(tname, ta))
        
        ta = table(table(m[, 2]))
        tname = names(ta)
        tname = as.integer(tname)
        r2 = c(r2, rep(tname, ta))
                
        i = i+1
        
        # number of nodes
        gn[i, 1] = length(node)
        # number of genes
        gn[i, 2] = length(unique(m[, 2]))
        
    }
    
    op = par(no.readonly = TRUE)
    par(mfrow=c(1,3))
    t1 = table(r1)
    plot((as.integer(names(t1))), (as.vector(t1)), pch=16, cex=0.8, xlab="Number of member genes in each node", ylab="Frequency", log="y", main="(A) Distribution of the number\nof member genes in each node")
    t2 = table(r2)
    plot((as.integer(names(t2))), (as.vector(t2)), pch=16, cex=0.8, xlab="Number of nodes in which a single gene resides", ylab="Frequency", log="y", main="(B) Distribution of the number\nof nodes in which a single gene resides")
    plot(gn[,2], gn[,1], pch=16, cex=0.8, xlim=range(gn), ylim=range(gn), xlab="Count of genes", ylab="Count of nodes", main="(C) Relationship between node count\nand gene count in pathways")
    abline(a=0,b=1)
    
    par(op)
    
}

is.ora = function(x) {
    return(class(x) == "cepa.all" && x[[1]][[1]]$framework == "ora")
}

is.gsa = function(x) {
    return(class(x) == "cepa.all" && x[[1]][[1]]$framework == "gsa.univariate")
}

