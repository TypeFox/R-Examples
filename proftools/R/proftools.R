.EmptyEnv <- if (exists("emptyenv")) emptyenv() else NULL
mkHash <- function() new.env(hash = TRUE, parent = .EmptyEnv)

n2d <- function(name, color = NULL, fontSize = 14, shape = "ellipse") {
    if (is.null(color) || is.na(color))
        paste0("\"", name, "\"[shape=", shape, ",fontsize=", fontSize,"];\n")
    else
        ## If you don't specify fillcolor, graphviz defaults both background &
        ## border to color. Original dot graph only specified color, which I 
        ## believe is the reason you made the dot file colors range from light 
        ## blue to red (to prevent borders from disappearing).
        ## See http://stackoverflow.com/questions/9106079/graphviz-how-to-change-border-color
        paste0("\"", name, "\"[shape=", shape,
               ",style=filled,color=black,fillcolor=\"",
               color, "\"fontsize=", fontSize, "];\n")
}

e2d <- function(from, to, color = NULL, edgeLabel="", edgeWidth = 1) {
    e <- paste0("\"", from, "\" -> \"", to, "\"")
    if (is.null(color))
        paste0(e, "[label=\"", edgeLabel, "\", penwidth=", edgeWidth, "];\n")
    else
        paste0(e, "[color=\"", color, "\", label=\"", edgeLabel,"\",
               penwidth=", edgeWidth, "];\n")
}

# **** A plausible size is 10,7.5
g2d <- function(g, filename = "g.dot", landscape = TRUE,
                nodeColors = NULL, nodeDetails = TRUE,
                edgeColors = NULL, edgeDetails = TRUE,
                size, center = FALSE, rankdir = c("TB","LR"),
                score = NULL, shape = "ellipse",
                nodeSizeScore = c("none", "total", "self"),
                edgeSizeScore = c("none", "total")) {
    rankdir <- match.arg(rankdir)
    nodeSizeScore <- match.arg(nodeSizeScore)
    edgeSizeScore <- match.arg(edgeSizeScore)

    eL <- g$edgeL
    
    con <- file(filename, open = "w")
    on.exit(close(con))

    cat("digraph xyz {\n", file = con)
    if (! missing(size))
        cat(paste0("size=\"", size, "\";\n"), file = con)
    if (landscape)
        cat("rotate=90;\n", file = con)
    if (center)
        cat("center=1;\n", file = con)
    cat(paste0("rankdir=", rankdir, ";\n"), file = con)

    ## The quadruple backslashes needed for dot file
    g$gNodes <- gsub("\n", '\\\\n', g$gNodes)
    edgeCounts <- g$callCounts
    if (edgeSizeScore != "none")
        edgeWidths <- lapply(edgeCounts, function(x) round(log10(x+10)))
    else
        edgeWidths <- lapply(edgeCounts, function(x) rep(1, length(x)))
    if (edgeDetails)
        edgeCounts <- lapply(edgeCounts, function(x) paste(" ", x))
    else
        edgeCounts <- lapply(edgeCounts, function(x) rep("", length(x)))
    ## This assigns the sizes of fonts, which determines the sizes of nodes
    if (nodeDetails)
        labels <- g$gNodes
    else
        labels <- g$nodes
    if (nodeSizeScore != "none") {
        if (nodeSizeScore == "total") val <- g$totalPercent
        else val <- g$selfPercent
        fontSizes <- pmax(7, sqrt(val) * 4) * 2
    }
    else
        fontSizes <- rep(14, length(labels))

    for (i in seq(along = labels)) {
        from <- labels[i]
        cat(n2d(from, nodeColors[[i]], fontSizes[i], shape), file = con)
        ## Have to use eL to know toList because g$edges contains plain
        ## node names
        toList <- labels[eL[[i]]$edges]
        toColors <- edgeColors[[i]]
        for (j in seq(along = toList))
            cat(e2d(from, toList[[j]], toColors[[j]], edgeCounts[[i]][j],
                    edgeWidths[[i]][j]), file = con)
    }

    cat("}", file = con)
}

g2g <- function(g, nodeDetails) {
    eL <- g$edgeL

    if (nodeDetails) {
        names(eL) <- g$gNodes
        nodes <- g$gNodes
    }
    else {
        nodes <- g$nodes
        names(eL) <- nodes
    }

    graph::graphNEL(nodes = nodes, edgeL = eL, edgemode = "directed")
}

getProfCallGraphNodeEntry <- function(name, env)
    get(name, envir = env)

incProfCallGraphNodeEntry <- function(name, what, env, count = 1) {
    if (exists(name, envir = env, inherits = FALSE))
        entry <- get(name, envir = env)
    else
        entry <- list(self = 0, total = 0, edges = mkHash())
    entry[[what]] <- entry[[what]] + count
    assign(name, entry, envir = env)
}

getProfCallGraphEdgeEntry <- function(from, to, env) {
    fromEntry <- getProfCallGraphNodeEntry(from, env)
    get(to, envir = fromEntry$edges)
}

incProfCallGraphEdgeEntry <- function(from, to, what, env, count = 1) {
    ## To allow for node trimming, only increment if both endpoints exist.
    if (exists(from, env) && exists(to, env)) {
        fromEntry <- getProfCallGraphNodeEntry(from, env)
        if (exists(to, envir = fromEntry$edges, inherits = FALSE))
            entry <- get(to, envir = fromEntry$edges)
        else entry <- list(self = 0, total = 0)
        entry[[what]] <- entry[[what]] + count
        assign(to, entry, envir = fromEntry$edges)
    }
}

rawProfCallGraph <- function(pd) {
    data <- mkHash()
    rvStacks <- lapply(pd$stacks, rev)
    fun <- function(line, count = 1) {
        incProfCallGraphNodeEntry(line[1], "self", data, count)
        for (n in unique(line))
            incProfCallGraphNodeEntry(n, "total", data, count)
        if (length(line) > 1) {
            incProfCallGraphEdgeEntry(line[2], line[1], "self",
                                      data, count)
            le <- lineEdges(line)
            for (i in seq(along = le$nodes)) {
                from <- le$nodes[i]
                for (to in le$edges[[i]])
                    incProfCallGraphEdgeEntry(from, to, "total",
                                              data, count)
            }
        }
    }
    mapply(fun, rvStacks, pd$counts, SIMPLIFY = FALSE)
    data
}

charMatch <- function(x, table, nomatch = NA)
    match(x, table, nomatch)

isIn <- function(x, table)
    match(x, table, 0)

lineEdges <- function(line) {
    if (length(line) > 1) {
        from <- unique(line[-1])
        edges <- rep(list(character(0)), length(from))
        for (i in 2 : length(line)) {
            j <- charMatch(line[i], from)
            if (! isIn(line[i - 1], edges[[j]]))
                edges[[j]] <- c(edges[[j]], line[i - 1])
        }
        list(nodes = from, edges = edges)
    }
}

lsEnv <- function(env)
    ls(env, all.names = TRUE)

profCallGraphEdges <- function(data) {
    edges <- mkHash()
    for (from in lsEnv(data)) {
        entry <- getProfCallGraphNodeEntry(from, data)
        assign(from, lsEnv(entry$edges), envir = edges)
    }
    edges
}

makeCycleMap <- function(cycles) {
    cycleMap <- mkHash()
    for (i in seq(along = cycles))
        for (n in cycles[[i]])
            assign(n, paste0("<cycle ", i, ">"), envir = cycleMap)
    cycleMap
}

## New findCycles
## Uses an incidence matrix and matrix multiplication.
## Could use a sparse representation but proably not worth the trouble.
edges2funs <- function(edges) {
    ln <- as.list(edges, all = TRUE)
    sort(unique(c(names(ln), unlist(ln))))
}

edges2mat <- function(funs, edges) {
    m <- diag(length(funs))
    for (i in seq_along(funs))
        m[i, match(edges[[funs[i]]], funs)] <- 1
    m
}

findReachableMat <- function(m) {
    nr <- nrow(m)
    n <- 2
    m <- m %*% m
    while (n <= nr) {
        n <- 2 * n
        m <- sign(m %*% m)
    }
    m
}

findMatCycles <- function(m) {
    nr <- nrow(m)
    cycles <- NULL
    visited <- rep(FALSE, nr)

    for (i in seq_len(nr - 1)) {
        if (! visited[i]) {
            visited[i] <- TRUE ## not really needed
            v <- i
            for (j in (i + 1) : nr)
                if (m[i, j] > 0 && m[j, i] > 0) {
                    v <- c(v, j)
                    visited[j] <- TRUE
                }
            if (length(v) > 1)
                cycles <- c(cycles, list(v))
        }
    }
    cycles
}

findCycles <- function(data) {
    edges <- profCallGraphEdges(data)
    funs <- edges2funs(edges)
    m1 <- edges2mat(funs, edges)
    mr <- findReachableMat(m1)
    lapply(findMatCycles(mr), function(v) funs[v])
}

compressLineRuns <- function(line) {
    if (length(line) > 1) {
        keep <- rep(TRUE, length(line))
        last <- line[1]
        for (i in 2 : length(line)) {
            val <- line[i]
            if (val == last)
                keep[i] <- FALSE
            else last <- val
        }
        line[keep]
    }
    else line
}

addCycleInfo <- function(pd, data, cycles) {
    map <- makeCycleMap(cycles)
    rvStacks <- lapply(pd$stacks, rev)
    inCycle <- function(name) exists(name, envir = map, inherits = FALSE)
    cycleName <- function(name) get(name, envir = map, inherits = FALSE)
    renameCycles <- function(line)
        unlist(lapply(line,
                      function(n) if (inCycle(n)) cycleName(n) else n))
    # **** speed up by inlining loop and calls to 'exists', 'get'
    renameCycles <- function(line) {
        len <- length(line)
        if (len > 0)
            for (i in 1 : len) {
                n <- line[i]
                ## if (.Internal(exists(n, map, "any", FALSE)))
                ##     line[i] <- .Internal(get(n, map, "any", FALSE))
                if (exists(n, envir = map, inherits = FALSE))
                    line[i] <- get(n, envir = map, inherits = FALSE)
            }
        line
    }
    cnames <- unique(unlist(lapply(lsEnv(map), get, map)))
    fun <- function(line, count = 1) {
        line <- compressLineRuns(renameCycles(line))
        if (isIn(line[1], cnames))
            incProfCallGraphNodeEntry(line[1], "self", data, count)
        for (n in unique(line))
            if (isIn(n, cnames))
                incProfCallGraphNodeEntry(n, "total", data, count)
        if (length(line) > 1) {
            if (isIn(line[1], cnames) || isIn(line[2], cnames))
                incProfCallGraphEdgeEntry(line[2], line[1], "self",
                                          data, count)
            le <- lineEdges(line)
            for (i in seq(along = le$nodes)) {
                from <- le$nodes[i]
                for (to in le$edges[[i]])
                    if (isIn(from, cnames) || isIn(to, cnames))
                        incProfCallGraphEdgeEntry(from, to, "total",
                                                  data, count)
            }
        }
    }
    mapply(fun, rvStacks, pd$counts, SIMPLIFY = FALSE)
}

trimProfCallGraph <- function(data, maxnodes, total.pct, total) {
    v <- unlist(eapply(data, function(x) x$total, all.names = TRUE))

    if (! is.na(maxnodes) && length(v) > maxnodes)
        drop <- names(sort(v)[1 : (length(v) - maxnodes)])
    else
        drop <- character(0)
    if (total.pct > 0)
        drop <- union(drop, names(v)[v < total * (total.pct / 100)])

    if (length(drop) > 0) {
        rm(list = drop, envir = data)
        for (fun in lsEnv(data)) {
            x <- get(fun, data, inherits = FALSE)
            e <- x$edges
            rm(list = intersect(lsEnv(e), drop), envir = e)
        }
    }

    data
}

cvtProfileData <- function(pd, GC, maxnodes = NA, total.pct = 0) {
    if (inherits(pd, "proftools_profData")) {
        if (GC && pd$haveGC)
            pd <- mergeGC(pd)
        data <- rawProfCallGraph(pd)
        if (! is.na(maxnodes) || total.pct > 0)
            data <- trimProfCallGraph(data, maxnodes, total.pct, pd$total)
        cycles <- findCycles(data)
        if (! is.null(cycles))
            addCycleInfo(pd, data, cycles)
        ## **** add a class -- proftools_callgraph?
        list(interval = pd$interval, count = pd$total,
             data = data, cycles = cycles)
    }
    else pd
}

profileDataCycles <- function(pd, GC)
    cvtProfileData(pd, GC)$cycles

revProfCallGraphMap <- function(data) {
    rg <- mkHash()
    for (from in lsEnv(data)) {
        entry <- getProfCallGraphNodeEntry(from, data)
        for (to in lsEnv(entry$edges)) {
            if (exists(to, envir = rg, inherits = FALSE))
                edges <- get(to, envir = rg)
            else edges <- character(0)
            if (! from %in% edges)
                assign(to, c(from, edges), envir = rg)
        }
    }
    rg
}

flatProfile <- function(pd, byTotal = TRUE, GC = TRUE) {
    pd <- cvtProfileData(pd, GC)
    nodes <- lsEnv(pd$data)
    if (! is.null(pd$cycles)) {
        map <- makeCycleMap(pd$cycles)
        cnames <- unique(unlist(lapply(lsEnv(map), get, map)))
        nodes <- nodes[! nodes %in% cnames]
    }

    if (byTotal) {
        total <- sapply(nodes, function(n) get(n, envir = pd$data)$total)
        ord <- order(-total)
    }
    else {
        self <- sapply(nodes, function(n) get(n, envir = pd$data)$self)
        ord <- order(-self)
    }

    nodes <- nodes[ord]
    self <- sapply(nodes, function(n) get(n, envir = pd$data)$self)
    selftime <- self * pd$interval/1e+06
    selfpct <- 100 * self / pd$count
    total <- sapply(nodes, function(n) get(n, envir = pd$data)$total)
    totaltime <- total * pd$interval/1e+06
    totalpct <- 100 * total / pd$count
    
    if (byTotal) {
        val <- cbind(round(totalpct, 2), round(totaltime, 2),
                     round(selfpct, 2), round(selftime, 2))
        colnames(val) <- c("total.pct", "total.time",
                           "self.pct", "self.time")
    }
    else {
        cumselftime <- cumsum(selftime)
        val <- cbind(round(selfpct, 2), round(cumselftime, 2),
                     round(selftime, 2), round(totalpct, 2),
                     round(totaltime, 2))
        colnames(val) <- c("self.pct", "cum.self.time", "self.time",
                           "total.pct", "total.time")
    }
    rownames(val) <- nodes
    val
}

makePrimaryLine<- function(node, i, pg) {
    idx <- sprintf("%-6s", paste0("[", i, "]"))
    if (pg$percent) {
        self <- pg$selfpct[i]
        child <- pg$childpct[i]
    }
    else {
        self <- pg$selftime[i]
        child <- pg$childtime[i]
    }
    stats <- sprintf("%8.2f   %8.2f   %8.2f", pg$totalpct[i], self, child)
    if (node %in% pg$cnames)
        name <- paste(substr(node, 1, nchar(node) - 1), "as a whole>")
    else if (node == "<Anonymous>")
        name <- "[Anonymous]"
    else name <- node
    if (pg$inCycle(node))
        extra <- paste(pg$cycleName(node), idx)
    else extra <- idx
    paste(idx, stats, "   ", name, extra, "\n")
}

makeCallerLine <- function(n, node, pg) {
    idx <- paste0("[", match(n, pg$nodes), "]")
    if (pg$inCycle(n) && pg$inCycle(node) &&
        pg$cycleName(n) == pg$cycleName(node))
        stats <- "                                     "
    else {
        entry <- getProfCallGraphEdgeEntry(n, node, pg$data)
        if (pg$percent) {
            self <- 100 * entry$self / pg$count
            total <- 100 * entry$total / pg$count
        }
        else {
            self <- entry$self * pg$interval/1e+06
            total <- entry$total * pg$interval/1e+06
        }
        child <- total - self
        stats <- sprintf("                  %8.2f   %8.2f",self, child)
    }
    if (pg$inCycle(n))
        extra <- paste(pg$cycleName(n), idx)
    else extra <- idx
    if (n == "<Anonymous>")
        name <- "[Anonymous]"
    else name <- n
    paste(stats, "       ", name, extra, "\n")
}

# **** most of this is the same as for callers--extract the common part.
makeCalleeLine <- function(n, node, pg) {
    idx <- paste0("[", match(n, pg$nodes), "]")
    if (pg$inCycle(n) && pg$inCycle(node) &&
        pg$cycleName(n) == pg$cycleName(node))
        stats <- "                                     "
    else {
        entry <- getProfCallGraphEdgeEntry(node, n, pg$data)
        if (pg$percent) {
            self <- 100 * entry$self / pg$count
            total <- 100 * entry$total / pg$count
        }
        else {
            self <- entry$self * pg$interval/1e+06
            total <- entry$total * pg$interval/1e+06
        }
        child <- total - self
        stats <- sprintf("                  %8.2f   %8.2f",self, child)
    }
    if (pg$inCycle(n))
        extra <- paste(pg$cycleName(n), idx)
    else extra <- idx
    if (n == "<Anonymous>")
        name <- "[Anonymous]"
    else name <- n
    paste(stats, "       ", name, extra, "\n")
}

makeCycleMemberLine <- function(n, cycle, pg) {
    i <- match(n, pg$nodes)
    idx <- paste0("[", i, "]")
    extra <- paste(cycle, idx)
    if (pg$percent) self <- pg$selfpct[i]
    else self <- pg$selftime[i]
    stats <- sprintf("                  %8.2f           ", self)
    if (n == "<Anonymous>")
        name <- "[Anonymous]"
    else name <- n
    paste(stats, "       ", name, extra, "\n")
}

printProfileCallGraph <- function(pd, file = stdout(), percent = TRUE,
                                  GC = TRUE, maxnodes = NA, total.pct = 0) {
    pd <- cvtProfileData(pd, GC, maxnodes, total.pct)
    if (is.character(file)) {
        if (file == "")
            stop("'file' must be non-empty string")
        con <- file(file, "wb")
        on.exit(close(con))
    }
    else if (inherits(file, "connection"))
        con <- file
    else stop("bad file argument")

    map <- makeCycleMap(pd$cycles)
    if (is.null(pd$cycles))
        cnames <- character(0)
    else cnames <- unique(unlist(lapply(lsEnv(map), get, map)))
    inCycle <- function(name) exists(name, envir = map, inherits = FALSE)
    cycleName <- function(name) get(name, envir = map, inherits = FALSE)

    nodes <- lsEnv(pd$data)
    total <- sapply(nodes, function(n) get(n, envir = pd$data)$total)
    ord <- order(-total)
    nodes <- nodes[ord]
    total <- total[ord]
    totalpct <- 100 * total / pd$count
    totaltime <- total * pd$interval/1e+06
    self <- sapply(nodes, function(n) get(n, envir = pd$data)$self)
    selfpct <- 100 * self / pd$count
    selftime <- self * pd$interval/1e+06
    pge <- profCallGraphEdges(pd$data)
    rpge <- revProfCallGraphMap(pd$data)

    pd$cnames <- cnames
    pd$nodes <- nodes
    pd$totalpct <- totalpct
    pd$selftime <- selftime
    pd$childtime <- totaltime - selftime
    pd$selfpct <- selfpct
    pd$childpct <- totalpct - selfpct
    pd$inCycle <- inCycle
    pd$cycleName <- cycleName
    pd$percent = percent

    cat("Call graph\n\n", file = con)
    if (percent)
        cat("index    % time     % self   % children     name\n\n", file = con)
    else
        cat("index    % time       self   children       name\n\n", file = con)
    for (i in seq(along = nodes)) {
        node <- nodes[i]
        if (exists(node, envir = rpge, inherits = FALSE))
            for (n in get(node, envir = rpge))
                if (! n %in% cnames)
                    cat(makeCallerLine(n, node, pd), file = con)
        cat(makePrimaryLine(node, i, pd), file = con)
        if (node %in% cnames)
            for (n in lsEnv(map))
                if (cycleName(n) == node)
                    cat(makeCycleMemberLine(n, node, pd), file = con)
        if (exists(node, envir = pge, inherits = FALSE))
            for (n in get(nodes[i], envir = pge))
                if (! n %in% cnames)
                    cat(makeCalleeLine(n, node, pd), file = con)
        cat("-----------------------------------------------\n", file = con)
    }
}

getOmittedNodes <- function(pd, mergeCycles) {
    map <- makeCycleMap(pd$cycles)
    cnodes <- lsEnv(map)
    if (mergeCycles)
        cnodes
    else if (is.null(pd$cycles))
        character(0)
    else unique(unlist(lapply(cnodes, get, map)))
}

extractProfileNodes <- function(pd, score = c("self", "total", "none"),
                                mergeCycles = TRUE) {
    if (score == "none")
        score <- "total"
    else match.arg(score)
    nodes <- lsEnv(pd$data)
    omitted <- getOmittedNodes(pd, mergeCycles)
    nodes <- nodes[! nodes %in% omitted]
    getScore <- function(n, type) get(n, envir = pd$data)[[type]]
    ## totalCost & selfCost needed for Google-style node Labels
    totalCost <- unlist(lapply(nodes, getScore, "total"))
    selfCost <- unlist(lapply(nodes, getScore, "self"))
    if(score == "total")
        sval <- totalCost / pd$count
    else sval <- selfCost / pd$count
    list(nodes = nodes, scores = sval,
         totalCost = totalCost, selfCost = selfCost)
}

extractProfileEdges <- function(pd, score = c("self", "total", "none"),
                                mergeCycles = TRUE) {
    if (score == "none")
        score <- "total"
    else match.arg(score)
    nodes <- lsEnv(pd$data)
    omitted <- getOmittedNodes(pd, mergeCycles)
    nodes <- nodes[! nodes %in% omitted]
    getToNodes <- function(n) {
        to <- lsEnv(get(n, envir = pd$data)$edges)
        to[! to %in% omitted]
    }
    edges <- lapply(nodes, getToNodes)
    ## getScores no longer divides by pd$count so can we avoid calling it 
    ## twice (if possible) to get callCounts below
    getScores <- function(n, type) {
        env <- get(n, envir = pd$data)$edges
        to <- lsEnv(env)
        to <- to[! to %in% omitted]
        unlist(lapply(to, function(v) get(v, envir = env)[[type]]))
    }
    sval <- lapply(nodes, getScores, score)
    if (score == "total") callCounts <- sval
    else callCounts <- lapply(nodes, getScores, "total")
    list(edges = edges, scores = sval, callCounts = callCounts)
}

np2x <- function(pd, score = c("total", "self", "none"),
                 transfer = function(x) x, nodeColorMap = NULL,
                 edgeColorMap = NULL, mergeCycles = FALSE,
                 edgesColored = FALSE) {
    match.arg(score)
    nodes <- extractProfileNodes(pd, score, mergeCycles = mergeCycles)
    edges <- extractProfileEdges(pd, score, mergeCycles = mergeCycles)
    totalPercent <- round(nodes$totalCost*100/pd$count, 2)
    selfPercent <- round(nodes$selfCost*100/pd$count, 2)
    gNodes <- paste0(nodes$nodes, "\n", nodes$selfCost, " (", selfPercent, 
                     "%) \n of ", nodes$totalCost, " (", totalPercent, "%)")
    p <- list(nodes = nodes$nodes, edges = edges$edges, 
              callCounts = edges$callCounts, gNodes = gNodes, 
              totalPercent = totalPercent, selfPercent = selfPercent)
    if (score == "none")
        color <- ecolor <- NULL
    else {
        ## Scale by maxScore to always create a red node
        maxScore <- max(nodes$scores)
        nodes$scores <- nodes$scores/maxScore
        color <- lapply(transfer(nodes$scores), colorScore, nodeColorMap)
        if (edgesColored) {
            ecolor <- vector("list", length(p$nodes))
            for (i in seq(along = ecolor)) {
                ## edges$scores is actually counts (see extractProfileEdges),
                ## divide by total count to convert to scores 
                ## and scale by maxScore to always create a red edge
                escore <- transfer(edges$scores[[i]] / pd$count)
                escore <- escore/maxScore
                ecolor[[i]] <- lapply(escore, colorScore, edgeColorMap)
            }
        }
        else ecolor <- NULL
    }
    p$nodeColors <- color
    p$edgeColors <- ecolor

    mke <- function(e) list(edges = match(e, p$nodes))
    p$edgeL <- lapply(p$edges, mke)

    p
}

colorScore <- function(score, colorMap) {
    if (is.null(score) || is.na(score))
        NULL
    else if (! is.null(colorMap)) {
        nc <- length(colorMap)
        colorMap[min(nc, max(ceiling(nc * (1 - score)), 1))]
    }
    else {
        score = min(max(score, 0), 1)
        # from cgprof
        maxhue = 0.6    # from red (.0) to magenta (.6), cf rainbow
        minsat = 0.1    # low saturation
        bri = 1.0       # brightness, always 100%

        # following formulas are totally empirical
        hue <- maxhue * (1.0 - score)
        sat <- minsat + (3.0 - minsat) * score
        paste0(hue, ",", sat, ",", bri)
    }
}

profileCallGraph2Dot <- function(pd, score = c("none", "total", "self"),
                                 transfer = function(x) x, nodeColorMap = NULL,
                                 edgeColorMap = NULL, filename = "Rprof.dot",
                                 landscape = FALSE, mergeCycles = FALSE,
                                 edgesColored = FALSE,
                                 rankDir = c("TB", "LR"),
                                 nodeDetails = FALSE, edgeDetails = FALSE,
                                 nodeSizeScore = c("none", "total", "self"),
                                 edgeSizeScore = c("none", "total"),
                                 center = FALSE, size, shape = "ellipse",
                                 layout = "dot", style, GC = TRUE,
                                 maxnodes = NA, total.pct = 0) {
    pd <- cvtProfileData(pd, GC, maxnodes, total.pct)

    if (missing(style)) style <- default.style

    if (missing(layout)) layout <- style$layout
    if (missing(score)) score <- style$score
    if (missing(transfer)) transfer <- style$transfer
    if (missing(nodeColorMap)) nodeColorMap <- style$nodeColorMap
    if (missing(edgeColorMap)) edgeColorMap <- style$edgeColorMap
    if (missing(mergeCycles)) mergeCycles <- style$mergeCycles
    if (missing(edgesColored)) edgesColored <- style$edgesColored
    if (missing(rankDir)) rankDir <- style$rankDir
    if (missing(nodeDetails)) nodeDetails <- style$nodeDetails
    if (missing(edgeDetails)) edgeDetails <- style$edgeDetails
    if (missing(nodeSizeScore)) nodeSizeScore <- style$nodeSizeScore
    if (missing(edgeSizeScore)) edgeSizeScore <- style$edgeSizeScore
    if (missing(shape)) shape <- style$shape
    if (missing(maxnodes)) maxnodes <- style$maxnodes
    if (missing(total.pct)) total.pct <- style$total.pct

    score <- match.arg(score)
    rankDir <- match.arg(rankDir)
    nodeSizeScore <- match.arg(nodeSizeScore)
    edgeSizeScore <- match.arg(edgeSizeScore)

    if (score != "none") {
        if (is.null(nodeColorMap))
            nodeColorMap <- heat.colors(100)
        if (is.null(edgeColorMap))
            edgeColorMap <- hsv(1,1,seq(1,0,length.out=50))
    }

    p <- np2x(pd, score, transfer, nodeColorMap, edgeColorMap, mergeCycles,
              edgesColored)
    g2d(p, filename, nodeColors = p$nodeColors, edgeColors = p$edgeColors,
        landscape = landscape, rankdir = rankDir, size = size, center = center,
        score = score, nodeDetails = nodeDetails, edgeDetails = edgeDetails,
        shape = shape, nodeSizeScore = nodeSizeScore,
        edgeSizeScore = edgeSizeScore)
}

plotProfileCallGraph <- function(pd, layout = "dot",
                                 score = c("none", "total", "self"),
                                 transfer = function(x) x, nodeColorMap = NULL,
                                 edgeColorMap = NULL, mergeCycles = FALSE,
                                 edgesColored = FALSE, rankDir = c("TB", "LR"),
                                 nodeDetails = FALSE, edgeDetails = FALSE,
                                 nodeSizeScore = c("none", "total", "self"),
                                 edgeSizeScore = c("none", "total"),
                                 shape = "ellipse", style, GC = TRUE,
                                 maxnodes = NA, total.pct = 0, ...) {
    if (missing(style)) style <- default.style

    if (missing(layout)) layout <- style$layout
    if (missing(score)) score <- style$score
    if (missing(transfer)) transfer <- style$transfer
    if (missing(nodeColorMap)) nodeColorMap <- style$nodeColorMap
    if (missing(edgeColorMap)) edgeColorMap <- style$edgeColorMap
    if (missing(mergeCycles)) mergeCycles <- style$mergeCycles
    if (missing(edgesColored)) edgesColored <- style$edgesColored
    if (missing(rankDir)) rankDir <- style$rankDir
    if (missing(nodeDetails)) nodeDetails <- style$nodeDetails
    if (missing(edgeDetails)) edgeDetails <- style$edgeDetails
    if (missing(nodeSizeScore)) nodeSizeScore <- style$nodeSizeScore
    if (missing(edgeSizeScore)) edgeSizeScore <- style$edgeSizeScore
    if (missing(shape)) shape <- style$shape
    if (missing(maxnodes)) maxnodes <- style$maxnodes
    if (missing(total.pct)) total.pct <- style$total.pct

    score <- match.arg(score)
    rankDir <- match.arg(rankDir)
    nodeSizeScore <- match.arg(nodeSizeScore)
    edgeSizeScore <- match.arg(edgeSizeScore)

    if (score != "none") {
        if (is.null(nodeColorMap))
            nodeColorMap <- heat.colors(100)
        if (is.null(edgeColorMap))
            edgeColorMap <- hsv(1,1,seq(1,0,length.out=50))
    }

    pd <- cvtProfileData(pd, GC, maxnodes, total.pct)

    p <- np2x(pd, score, transfer, nodeColorMap,
              edgeColorMap, mergeCycles, edgesColored)

    g <- g2g(p, nodeDetails)
    names(labels) <- labels <- graph::nodes(g)
    if (! is.null(p$nodeColors)) {
        p$nodeColors <- unlist(p$nodeColors)
        names(p$nodeColors) <- labels
    }
    ## Create the edgeNames list based on the edgeL(g) which only has
    ## numbers matching the from and to.
    edges <- graph::edgeL(g); edgeNames <- list()
    for (i in seq(along = edges))
        if(length(edges[[i]]$edges)>0)
            edgeNames[[i]] <- paste(labels[i], labels[edges[[i]]$edges],
                                    sep = "~")
    edgeNames <- unlist(edgeNames)
    if (! is.null(p$edgeColors)) {
        p$edgeColors <- unlist(p$edgeColors)
        names(p$edgeColors) <- edgeNames
    }

    if (nodeSizeScore != "none") {
        ## This uses the same font size calculation as the dot file
        ## version
        if (nodeSizeScore == "total") val <- p$totalPercent
        else val <- p$selfPercent
        fontSizes <- pmax(7, sqrt(val) * 4) * 2
        names(fontSizes) <- labels
        graph::nodeRenderInfo(g) <- list(label = labels, fontsize = fontSizes)
    }

    edgeCounts <- unlist(p$callCounts)
    anyEdges <- (length(edgeCounts) > 0)
    if (anyEdges) {
        if (edgeSizeScore != "none") {
            edgeWeights <- log10(edgeCounts+10)
            names(edgeWeights) <- edgeNames
            graph::edgeRenderInfo(g) <- list(lwd = edgeWeights)
        }

        if (edgeDetails) {
            spacing <- sapply(edgeCounts, function(x){rep("   ", length(x))})
            edgeCounts <- paste0(spacing, edgeCounts)
            names(edgeCounts) <- edgeNames
            graph::edgeRenderInfo(g) <- list(label = edgeCounts, fontsize = 8)
        }
    }

    ## The order of layout and rendering info calls below
    ## matters. Certain attributes don't work well if we change the
    ## order.  See
    ## https://stat.ethz.ch/pipermail/bioconductor/2012-November/049253.html
    attrs <- list(node = list(shape = shape, fixedsize = FALSE))
    if (layout == "dot")
        attrs$graph <- list(rankdir = rankDir)
    g <- Rgraphviz::layoutGraph(g, layoutType = layout, attrs = attrs)
    if (anyEdges)
        graph::edgeRenderInfo(g) <- list(col = p$edgeColors)
    graph::nodeRenderInfo(g) <- list(fill = p$nodeColors)
    tryCatch(Rgraphviz::renderGraph(g, ...),
             error = function(e) NULL) ## catch and ignore errors
}

plain.style <- list(layout = "dot", score = "none",
                    transfer = function(x) x, nodeColorMap = NULL,
                    edgeColorMap = NULL, mergeCycles = FALSE,
                    edgesColored = FALSE, rankDir = "TB",
                    nodeDetails = FALSE, edgeDetails = FALSE,
                    nodeSizeScore = "none", edgeSizeScore = "none",
                    shape = "ellipse", maxnodes = NA, total.pct = 0)

google.style <- list(layout = "dot", score = "none",
                     transfer = function(x) x, nodeColorMap = NULL,
                     edgeColorMap = NULL, mergeCycles = FALSE,
                     edgesColored = FALSE, rankDir = "TB",
                     nodeDetails = TRUE, edgeDetails = TRUE,
                     nodeSizeScore = "self", edgeSizeScore = "total",
                     shape = "box", maxnodes = NA, total.pct = 0)

default.style <- google.style
default.style$score <- "total"
default.style$maxnodes <- 30
