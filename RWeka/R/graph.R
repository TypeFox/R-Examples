write_to_dot <-
function(x, con = stdout(), ...)
    UseMethod("write_to_dot")

write_to_dot.Weka_classifier <-
function(x, con = stdout(), ...)
{
    ## Most Weka classifiers do not implement the 'Drawable' interface
    ## and hence have no graph() method.
    if(!.has_method(x$classifier, "graph"))
        stop("Cannot create dot description from 'x'.")
    writeLines(.jcall(x$classifier, "S", "graph"), con)
}

parse_Weka_digraph <-
function(x, plainleaf = TRUE)
{
    ## Create a simple node/edge representation of digraphs obtained
    ## from Weka's graph() methods.  Note that this could easily be
    ## turned into a full-featured dot parser possibly providing a
    ## graphNEL representation, but this is already done by package
    ## 'graph' ...
    
    ## Use the individual lines apart from the first and last ones.
    x <- strsplit(x, "\n")[[1L]]
    x <- x[-c(1L, length(x))]
    
    ind <- regexpr("->", x, fixed = TRUE)
    nodes <- x[ind == -1L]
    edges <- x[ind != -1L]
    
    nval <- matrix(rep("", 2L * length(nodes)), ncol = 2L)
    colnames(nval) <- c("name", "splitvar")
    nval[, 1L] <- sapply(strsplit(nodes, " "), "[", 1L)
    nval[, 2L] <- sapply(strsplit(nodes, "\""), "[", 2L)
    if(plainleaf)
        nval[grep("(", nval[, 2L], fixed = TRUE), 2L] <- ""
    
    eval <- matrix(rep("", 3L * length(edges)), ncol = 3L)
    colnames(eval) <- c("from", "to", "label")
    eval[, 1L] <- sapply(strsplit(edges, "->"), "[", 1L)
    eval[, 2L] <-
        sapply(strsplit(as.character(sapply(strsplit(edges, "->"),
                                            "[", 2L)),
                        " "),
               "[", 1L)
    eval[, 3L] <- sapply(strsplit(edges, "\""), "[", 2L)
    
    return(list(nodes = nval, edges = eval))
}

make_Weka_classifier_tree <-
function(x)
{
    ## For a fitted Weka classifier tree from a class which implements
    ## the Drawable interface and hence has a graph() method creating a
    ## dot representation, create an intermediate representation of the
    ## graph which can then be coerced to e.g. a dendrogram object (note
    ## that the plot method for dendrograms really is not good enough
    ## for our purposes), or a BinaryTree object (provided that the tree
    ## is binary, of course).

    x <- .jcall(x$classifier, "S", "graph")
    nodes_and_edges <- parse_Weka_digraph(x, FALSE)
    nodes <- nodes_and_edges$nodes
    edges <- nodes_and_edges$edges

    max_n_of_children <- if(NROW(edges) > 0L) max(table(edges[, "from"])) else 0
    max_depth <- 0
            
    get_subtree <- function(node, depth = 0) {
        if (depth > max_depth)
            max_depth <<- depth

        ind <- which(nodes[, "name"] == node)
        ## message(ind, "\n")
        ## Sanitize ...
        label <- nodes[ind, "splitvar"]
        ## message(node, label, "\n")
        ind <- which(edges[, "to"] == node)
        edgetext <- if(any(ind))
            edges[ind, "label"]
        else
            ""
        ind <- which(edges[, "from"] == node)
        if (!length(ind)) {
            out <- 1
            attributes(out) <-
                list(members = 1, leaf = TRUE, depth = depth,
                     label = label, edgetext = edgetext, nodeID = node)
            return(out)
        }
        out <- vector("list", length = length(ind))
        for (i in seq_along(out))
            out[[i]] <- Recall(edges[ind[i], "to"], depth + 1)
        attributes(out) <-
            list(members = sum(sapply(out, attr, "members")),
                 leaf = FALSE, depth = depth, label = label,
                 edgetext = edgetext, nodeID = node)
        out
    }

    out <- get_subtree("N0")
    attr(out, "max_depth") <- max_depth
    attr(out, "max_n_of_children") <- max_n_of_children
    
    out
}

as.dendrogram.Rweka_classifier_tree <-
function(object, ...)
{
    max_depth <- attr(object, "max_depth")
    attr(object, "max_depth") <- NULL

    convert <- function(x) {
        y <- x
        class(y) <- "dendrogram"
        attr(y, "height") <- max_depth - attr(x, "depth")
        if(is.leaf(x))
            return(y)
        for(i in seq_along(x))
            y[[i]] <- Recall(x[[i]])
        y
    }

    convert(object)
}

plot.Weka_tree <-
function(x, ...)
{
    if(system.file(package = "partykit") == "")
        stop("Plotting Weka trees requires package 'partykit'.")
    plot(partykit::as.party(x), ...)
}
