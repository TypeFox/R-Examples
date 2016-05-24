##
## @author Marc A. Suchard
##
## A class for reading Newick formatted trees with BEAST-style annotations
##

.strip.annotations <- function(text) {
    annotations <- list()
    end <- 1

    # Merge node and branch annotations
    text <- gsub("\\[&(.*?)\\]:\\[&(.*?)\\]", ":\\[&\\1,\\2\\]", text)
    text <- gsub("\\[&(.*?)\\]:", ":\\[&\\1]", text)

    pattern = "\\[&.*?\\]"

    repeat {
        match = regexpr(pattern=pattern,text=text)
        if (!(match[1] > 0)) {
            break
        }
        annotations[[end]] = regmatches(text, match)
        text = sub(pattern,paste("[",end,"]",sep=""), text)
        end = end + 1
    }
    return(list(annotations=annotations,tree=text))
}

.split.tree.names <- function(text) {
    text = gsub(pattern="\\[.*?\\]=",x=text,replacement="")
    text = gsub(pattern="^tree",x=text,replacement="")
    return(text)
}

.split.tree.traits <- function(text) {

    ## Pull out annotation
    text = regmatches(text,regexpr(pattern="\\[.*?\\]",text))
    ## Remove leading and trailing delimitors
    text = substring(text,3,nchar(text)-1)
    return(text)
}

.parse.value <- function(text) {
    value = text
    if (length(grep("^\\{",value))) { ## starts with {
        save = value
        value = substring(value, 2, nchar(value)-1)

        depth = 0
        r = regexpr(pattern="\\{+",value,perl=TRUE)
        match.length = attr(r, "match.length")

        if (match.length > 0) {
            depth = match.length
        }

        if (depth == 0) {
            split = ","
        } else {
            split = paste(
            "(?<=",rep("\\}",depth),")",
            ",",
            "(?=" ,rep("\\{",depth),")",
            sep="")
        }

        if (depth >= 1) {
            return(save) # TODO Still error in recursion
        }

        part = strsplit(value, split, perl=TRUE)[[1]]
        value = list()
        for (i in 1:length(part)) {
            value[[i]] = .parse.value(part[i])
        }
        ## TODO Unlist when simple array?
    } else {
        if (!is.na(suppressWarnings(as.numeric(value)))) { # is a number
            value = as.numeric(value)
        } else { # is a string
        	 value <- gsub("\\\"","", value)
        }
    }
    return(value)
}

.parse.traits <- function(text, header=FALSE) {

    if (header == TRUE) {
        text = substring(text,3,nchar(text)-1)
    }

    pattern <- "(\"[^\"]*\"+|[^,=\\s]+)\\s*(=\\s*(\\{[^=]*\\}|\"[^\"]*\"+|[^,]+))?"

    rgx <- gregexpr(pattern,text,perl=TRUE)
    n <- length(attr(rgx[[1]],"match.length"))
    traits <- list()
    start <- attr(rgx[[1]],"capture.start")
    names <- attr(rgx[[1]],"capture.names")
    length <- attr(rgx[[1]],"capture.length")
    names <- attr(rgx[[1]],"capture.names")
    for (i in 1:n) {
        s <- start[i,3]
        e <- s + length[i,3] - 1
        value <- substring(text,s,e)

        s <- start[i,1]
        e <- s + length[i,1] - 1
        key <- substring(text,s,e)

        traits[[key]] <- .parse.value(value)
    }

    return(traits)
}

.annotated.clado.build <- function(tp) {
    stop(paste("Annotated clado.build is not yet implemented.\n"))
}

## THE CODE BELOW COMES FROM 'ape'. MY GOAL IS TO DERIVE FROM THIS TO READ IN BEAST-STYLE ANNOTATIONS

.annotated.tree.build <- function(tp){

    add.internal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- current.node <<- node <<- node + 1L
        index[node] <<- j
        j <<- j + 1L
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j
        X <- unlist(strsplit(new.tpc[k], ":"))
        tip.label[tip] <<- X[1]
        index <- length(X)
        edge.length[j] <<- as.numeric(X[index])

        if (length(annotations) > 0) {
        	permute[[j]] <<- annotations[[as.numeric(X[2])]] ## permute traits
        }

        k <<- k + 1L
        tip <<- tip + 1L
        j <<- j + 1L
    }
    go.down <- function() {
        l <- index[current.node]
        X <- unlist(strsplit(new.tpc[k], ":"))
        node.label[current.node - nb.tip] <<- X[1]
        index <- length(X)
        edge.length[l] <<- as.numeric(X[index])

        if (length(annotations) >  0) {
        	permute[[l]] <<- annotations[[as.numeric(X[2])]] ## permute traits
        }

        k <<- k + 1L
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2L, 1L), 1, 2))
        tp <- unlist(strsplit(tp, "[\\(\\):;]"))
        obj$edge.length <- as.numeric(tp[3])
        obj$Nnode <- 1L
        obj$tip.label <- tp[2]
        if (tp[4] != "")
            obj$node.label <- tp[4]
        class(obj) <- "phylo"
        return(obj)
    }

    result = .strip.annotations(tp)
    annotations = result$annotations
    new.tp.stripped = result$tree

    # patched for 0.0 root branch length from BEAST2 (not confirmed)
    new.tp.stripped <- gsub("\\]0.0;", "\\];", new.tp.stripped)

    root.annotation.number <- NULL
    m <- regexpr("\\[\\d+\\];", new.tp.stripped)
    if (m != -1) {
    	root.annotation.number <- as.numeric(
    		gsub("\\[(\\d+)\\];", "\\1", regmatches(new.tp.stripped, m)))
    }

    annotations = lapply(annotations, .parse.traits, header=TRUE)

    tp.stripped = gsub("\\[.*?\\]","",tp)
    tpc <- unlist(strsplit(tp.stripped, "[\\(\\),;]"))
    tpc <- tpc[nzchar(tpc)]

    new.tp.stripped <- gsub("\\[\\d+\\];", ";", new.tp.stripped)
    new.tp.stripped <- gsub("\\[(\\d+)\\]","\\1:", new.tp.stripped)
    new.tpc <- unlist(strsplit(new.tp.stripped, "[\\(\\),;]"))
    new.tpc <- new.tpc[nzchar(new.tpc)]

    tsp <- unlist(strsplit(tp.stripped, NULL))
    skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
    nsk <- length(skeleton)
    nb.node <- sum(skeleton == ")")
    nb.tip <- sum(skeleton == ",") + 1
    nb.edge <- nb.node + nb.tip

    node.label <- character(nb.node)
    tip.label <- character(nb.tip)
    edge.length <- numeric(nb.edge)
    edge <- matrix(0L, nb.edge, 2)
    current.node <- node <- as.integer(nb.tip + 1)
    edge[nb.edge, 2] <- node
    index <- numeric(nb.edge + 1)
    index[node] <- nb.edge
    j <- k <- tip <- 1L

    permute = list()

    for (i in 2:nsk) {
        if (skeleton[i] == "(") {
            add.internal()
        }
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") {
                add.terminal()
            }
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] == ",") {
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") {
                go.down()
            }
        }
    }

    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
    root.edge <- edge.length[nb.edge]
    edge.length <- edge.length[-nb.edge]
    if (!all(is.na(edge.length)))
        obj$edge.length <- edge.length
    if (is.na(node.label[1]))
        node.label[1] <- ""
    if (any(nzchar(node.label)))
        obj$node.label <- node.label
    if (!is.na(root.edge))
        obj$root.edge <- root.edge
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"

    if (!is.null(root.annotation.number)) {
    	obj$root.annotation <- annotations[[root.annotation.number]]
    }

    obj$annotations = permute
    obj
}

read.annotated.tree <- function (file = "", text = NULL, tree.names = NULL, skip = 0,
                                  comment.char = "#", keep.multi = FALSE, ...)
{
    unname <- function(treetext) {
        nc <- nchar(treetext)
        tstart <- 1
        while (substr(treetext, tstart, tstart) != "(" && tstart <=
               nc) tstart <- tstart + 1
        if (tstart > 1)
            return(c(substr(treetext, 1, tstart - 1), substr(treetext,
                                                             tstart, nc)))
        return(c("", treetext))
    }

    if (!is.null(text)) {
        if (!is.character(text))
            stop("argument `text' must be of mode character")
        tree <- text
    }
    else {
        tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }

    tree <- gsub("[ \n\t]", "", tree)
    tree <- gsub("\\[&R\\]", "", tree)
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    if (is.na(y[1]))
        return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree) STRING[i] <- paste(tree[x[i]:y[i]], sep = "",
                                          collapse = "")

    tmp <- unlist(lapply(STRING, unname))
    tmpnames <- tmp[c(TRUE, FALSE)]
    STRING <- tmp[c(FALSE, TRUE)]
    if (is.null(tree.names) && any(nzchar(tmpnames)))
        tree.names <- tmpnames
    colon <- grep(":", STRING)

    if (!is.null(tree.names)) {
        traits.text = lapply(tree.names, .split.tree.traits)
        tree.names = lapply(tree.names, .split.tree.names)
        tree.traits = lapply(traits.text, .parse.traits)
    }

    if (!length(colon)) {
        stop(paste("Annotated clado.build is not yet implemented.\n"))
        obj <- lapply(STRING, .annotated.clado.build)
    }
    else if (length(colon) == Ntree) {
        obj <- lapply(STRING, .annotated.tree.build)
    }
    else {
        obj <- vector("list", Ntree)
        obj[colon] <- lapply(STRING[colon], .annotated.tree.build)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        obj[nocolon] <- lapply(STRING[nocolon], clado.build)
    }
    for (i in 1:Ntree) {
        ROOT <- length(obj[[i]]$tip.label) + 1
        if (sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] >
            1)
            stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading Newick file aborted at tree no.",
                       i))
    }
    if (Ntree == 1 && !keep.multi)
        obj <- obj[[1]]
    else {
        if (!is.null(tree.names)) {
            names(obj) <- tree.names
        }
        class(obj) <- "multiPhylo"
    }
    obj
}


read.annotated.nexus <- function (file, tree.names = NULL) {
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)

    ##    browser()
    ##
    ##    if (length(LEFT)) {
    ##        w <- LEFT == RIGHT
    ##        if (any(w)) {
    ##            s <- LEFT[w]
    ##            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
    ##        }
    ##        w <- !w
    ##        if (any(w)) {
    ##            s <- LEFT[w]
    ##            X[s] <- gsub("\\[.*", "", X[s])
    ##            sb <- RIGHT[w]
    ##            X[sb] <- gsub(".*\\]", "", X[sb])
    ##            if (any(s < sb - 1))
    ##                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
    ##        }
    ##    }

    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1)
        TRUE
    else FALSE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end]
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    start <- if (translation)
        semico[semico > i2][1] + 1
    else semico[semico > i1][1]
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]

    ##    browser()

    rm(X)
    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico)
    if (Ntree == 1 && length(tree) > 1)
        STRING <- paste(tree, collapse = "")
    else {
        if (any(diff(semico) != 1)) {
            STRING <- character(Ntree)
            s <- c(1, semico[-Ntree] + 1)
            j <- mapply(":", s, semico)
            if (is.list(j)) {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]],collapse = "")
            }
            else {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[,i]], collapse = "")
            }
        }
        else STRING <- tree
    }
    rm(tree)

    ##    browser()

    STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
    Ntree <- length(STRING)

    STRING <- gsub("\\[&R\\]", "", STRING)

    ## TODO Parse out tree-level traits
    nms.annontations.trees <- sub(" * = *.*", "", STRING)
    nms.annontations.trees <- sub("^ *tree *", "", nms.annontations.trees, ignore.case = TRUE)

    nms.trees <- sub("\\s+\\[&.*?\\]", "", nms.annontations.trees)

    if (any(nms.trees != nms.annontations.trees)) { # There are tree-level annontations
    	annotations.trees <- sub(".*\\[&", "\\[&", nms.annontations.trees)
    	annotations.trees = lapply(annotations.trees, .parse.traits, header=TRUE)
    } else {
    	annotations.trees <- NULL
    }

    STRING <- sub("^.*? = *", "", STRING)
    STRING <- gsub("\\s", "", STRING)

    ##    browser()

    colon <- grep(":", STRING)
    if (!length(colon)) {
        stop(".annotated.clado.build is not yet implemented.\n")
        trees <- lapply(STRING, .annotated.clado.build)
    } else if (length(colon) == Ntree) {
        ##        trees <- if (translation) {
        ##            browser()
        ##            stop("treeBuildWithTokens is not yet implemented.\n")
        ##            lapply(STRING, .treeBuildWithTokens)
        ##        }
        ##        else lapply(STRING, .annotated.tree.build)
        trees <- lapply(STRING, .annotated.tree.build)
        ##        browser()
    } else {
        ##        trees <- vector("list", Ntree)
        ##        trees[colon] <- lapply(STRING[colon], .annotated.tree.build)
        ##        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        ##        trees[nocolon] <- lapply(STRING[nocolon], .annotated.clado.build)
        ##        if (translation) {
        ##            for (i in 1:Ntree) {
        ##                tr <- trees[[i]]
        ##                for (j in 1:n) {
        ##                  ind <- which(tr$tip.label[j] == TRANS[, 1])
        ##                  tr$tip.label[j] <- TRANS[ind, 2]
        ##                }
        ##                if (!is.null(tr$node.label)) {
        ##                  for (j in 1:length(tr$node.label)) {
        ##                    ind <- which(tr$node.label[j] == TRANS[,
        ##                      1])
        ##                    tr$node.label[j] <- TRANS[ind, 2]
        ##                  }
        ##                }
        ##                trees[[i]] <- tr
        ##            }
        ##            translation <- FALSE
        ##        }
        stop("Unknown error in read.annotated.nexus.\n")
    }
    for (i in 1:Ntree) {
        tr <- trees[[i]]
        if (!translation)
            n <- length(tr$tip.label)
        ROOT <- n + 1
        if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] >
            1) {
            stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading NEXUS file aborted at tree no.",
                       i, sep = ""))
        }
    }
    if (Ntree == 1) {
        trees <- trees[[1]]
        if (translation) {
            trees$tip.label <- # if (length(colon)) {
                # TRANS[, 2]
            #} else {
            	TRANS[, 2][as.numeric(trees$tip.label)]
            #}
        }
    }
    else {
        if (!is.null(tree.names))
            names(trees) <- tree.names
        if (translation) {
            # if (length(colon) == Ntree)
              #  attr(trees, "TipLabel") <- TRANS[, 2]
            #else {
                for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
                trees <- .compressTipLabel(trees)
            #}
        }
        class(trees) <- "multiPhylo"
        if (!all(nms.trees == ""))
            names(trees) <- nms.trees
    }

    # Add tree-level annotations back on
    if (!is.null(annotations.trees)) {
    	if (Ntree == 1) {
    		trees$tree.annotations <- annotations.trees[[1]]
    	} else {
    		for (i in 1:Ntree) {
    			trees[[i]]$tree.annotations <- annotations.trees[[i]]
    		}
    	}
    }

    trees
} # end read.annotated.nexus

