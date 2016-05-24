## read.tree.R (2015-01-12)

##   Read Tree Files in Parenthetic Format

## Copyright 2002-2012 Emmanuel Paradis, Daniel Lawson and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

tree.build <- function(tp)
{
    add.internal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- current.node <<- node <<- node + 1L
        index[node] <<- j # set index
        j <<- j + 1L
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j # set index
        X <- unlist(strsplit(tpc[k], ":"))
        tip.label[tip] <<- X[1]
        edge.length[j] <<- as.numeric(X[2])
        k <<- k + 1L
        tip <<- tip + 1L
        j <<- j + 1L
    }
    go.down <- function() {
        l <- index[current.node]
        X <- unlist(strsplit(tpc[k], ":"))
        node.label[current.node - nb.tip] <<- X[1]
        edge.length[l] <<- as.numeric(X[2])
        k <<- k + 1L
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2L, 1L), 1, 2))
        tp <- unlist(strsplit(tp, "[\\(\\):;]"))
        obj$edge.length <- as.numeric(tp[3])
        obj$Nnode <- 1L
        obj$tip.label <- tp[2]
        if (tp[4] != "") obj$node.label <- tp[4]
        class(obj) <- "phylo"
        return(obj)
    }

    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[nzchar(tpc)]
    ## the following 2 lines are (slightly) faster than using gsub()
    tsp <- unlist(strsplit(tp, NULL))
    skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
    nsk <- length(skeleton)
    nb.node <- sum(skeleton == ")")
    nb.tip <- sum(skeleton == ",") + 1
    ## We will assume there is an edge at the root;
    ## if so, it will be removed and put into a vector
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)

    edge.length <- numeric(nb.edge)
    edge <- matrix(0L, nb.edge, 2)
    current.node <- node <- as.integer(nb.tip + 1) # node number
    edge[nb.edge, 2] <- node
    index <- numeric(nb.edge + 1) # hash index to avoid which
    index[node] <- nb.edge

    ## j: index of the line number of edge
    ## k: index of the line number of tpc
    ## tip: tip number
    j <- k <- tip <- 1L

    for (i in 2:nsk) {
        if (skeleton[i] == "(") add.internal() # add an internal branch (on top)
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") add.terminal() # add a terminal branch
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] == ",") { # add a terminal branch and go down one level
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") go.down() # go down one level
        }
    }

    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
    root.edge <- edge.length[nb.edge]
    edge.length <- edge.length[-nb.edge]
    if (!all(is.na(edge.length))) # added 2005-08-18
        obj$edge.length <- edge.length
    if (is.na(node.label[1])) node.label[1] <- ""
    if (any(nzchar(node.label))) obj$node.label <- node.label
    if (!is.na(root.edge)) obj$root.edge <- root.edge
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"
    obj
}

read.tree <- function(file = "", text = NULL, tree.names = NULL, skip = 0,
    comment.char = "#", keep.multi = FALSE, ...)
{
    unname <- function(treetext) {
        nc <- nchar(treetext)
	tstart <- 1
	while (substr(treetext, tstart, tstart) != "(" && tstart <= nc)
            tstart <- tstart + 1
	if (tstart > 1)
            return(c(substr(treetext, 1, tstart - 1),
                     substr(treetext, tstart, nc)))
	return(c("", treetext))
    }
    if (!is.null(text)) {
        if (!is.character(text))
          stop("argument `text' must be of mode character")
        tree <- text
    } else {
        tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }
    ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }
    tree <- gsub("[ \t]", "", tree)
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    ## Suggestion from Olivier Francois (added 2006-07-15):
    if (is.na(y[1])) return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree) {
        tmp <- paste(tree[x[i]:y[i]], sep = "", collapse = "")
        STRING[i] <- gsub("\\[[^]]*\\]", "", tmp) # delete comments (fix 2015-01-12)
    }

    tmp <- unlist(lapply(STRING, unname))
    tmpnames <- tmp[c(TRUE, FALSE)]
    STRING <- tmp[c(FALSE, TRUE)]
    if (is.null(tree.names) && any(nzchar(tmpnames)))
        tree.names <- tmpnames

    colon <- grep(":", STRING)
    if (!length(colon)) {
        obj <- lapply(STRING, clado.build)
    } else if (length(colon) == Ntree) {
        obj <- lapply(STRING, tree.build)
    } else {
        obj <- vector("list", Ntree)
        obj[colon] <- lapply(STRING[colon], tree.build)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        obj[nocolon] <- lapply(STRING[nocolon], clado.build)
    }
    for (i in 1:Ntree) {
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root
        ROOT <- length(obj[[i]]$tip.label) + 1
        if(sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 1)
            stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading Newick file aborted at tree no.", i))
    }
    if (Ntree == 1 && !keep.multi) obj <- obj[[1]] else {
        if (!is.null(tree.names)) names(obj) <- tree.names
        class(obj) <- "multiPhylo"
    }
    obj
}
