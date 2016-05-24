"simpleGraphToGraph" <-
function (sdg = NULL, frameModels = NULL, dg = NULL, 
          control = dg.control(...), ...) 
{
    if ((is.null(sdg))) 
        sdg <- .newDgSimpleGraph(...)
    else if ((class(sdg) != "dg.simple.graph")) 
        sdg <- .newDgSimpleGraph(vertex.names = sdg, ...)
    if (!is.null(sdg@types) && (length(sdg@types) > 0) && (sdg@types == "")) 
        sdg@types <- character()
    if (length(sdg@from) == 0) 
        sdg@from <- vector()
    if (length(sdg@to) == 0) 
        sdg@to <- vector()
    if (!is.null(sdg@edge.types) && 
        (length(sdg@edge.types) > 0) && (sdg@edge.types == "")) 
        sdg@edge.types <- character()
    if (.IsEmpty(sdg@edge.list)) 
        sdg@edge.list <- list(NULL)
    if (!is.null(sdg@labels) && 
        (length(sdg@labels) > 0) && (sdg@labels == "")) 
        sdg@labels <- sdg@vertex.names
    if (.IsEmpty(sdg@blocks)) 
        sdg@blocks <- list(NULL)
    if (.IsEmpty(sdg@block.tree)) 
        sdg@block.tree <- list(NULL)
    if (.IsEmpty(sdg@factors)) 
        sdg@factors <- list(NULL)
    if (!is.null(sdg@texts) && (length(sdg@texts) > 0) && (sdg@texts == "")) 
        sdg@texts <- character()
    if (length(sdg@extra.from) == 0) 
        sdg@extra.from <- vector()
    if (length(sdg@extra.to) == 0) 
        sdg@extra.to <- vector()
    if (.IsEmpty(sdg@extra.edge.list)) 
        sdg@extra.edge.list <- list(NULL)
    "two.to.pairs" <- function(from, to) {
        edge.list <- vector("list", length(to))
        for (j in seq(along = to)) edge.list[[j]] <- c(from[j], to[j])
        return(edge.list)
    }
    "check.names" <- function(x, label) if (length(x) > 0) 
        if (!is.numeric(x)) {
            x <- unique(sort(x))
            x.indices <- match(x, sdg@vertex.names)
            if (any(is.na(x.indices))) {
                new.names <- x[is.na(x.indices)]
                warning(paste("Invalid names in '", label, "' : ", 
                  paste(new.names, collapse = ", ")))
                sdg@vertex.names <<- c(sdg@vertex.names, new.names)
                sdg@types <<- c(sdg@types, rep(sdg@types[1], 
                  length(new.names)))
                vertex.colors <<- c(vertex.colors, rep("magenta", 
                  length(new.names)))
            }
        }
        else if (max(x) > length(sdg@vertex.names)) {
            warning(paste("Invalid index", max(x), "in '", label, "'"))
            new.names <- paste("v", (length(sdg@vertex.names) + 1):max(x), 
                sep = "")
            sdg@vertex.names <<- c(sdg@vertex.names, new.names)
            sdg@types <<- c(sdg@types, rep(sdg@types[1], length(new.names)))
            vertex.colors <<- c(vertex.colors, rep("magenta", 
                length(new.names)))
        }
    Arguments <- list(...)
    Vertices <- NULL
    if (!is.null(frameModels)) {
        Vertices <- frameModels@vertices
        sdg@vertex.names <- Names(Vertices)
    }
    else if (!is.null(dg) && !is.null(dg@vertexList)) {
        Vertices <- dg@vertexList
        sdg@vertex.names <- Names(Vertices)
    }
    else if (!is.null(Arguments$Vertices)) {
        Vertices <- Arguments$Vertices
        sdg@vertex.names <- Names(Vertices)
    }
    else if (!is.null(Arguments$vertexList)) {
        Vertices <- Arguments$vertexList
        sdg@vertex.names <- Names(Vertices)
    }
    n <- length(sdg@vertex.names)
    if (is.numeric(sdg@vertex.names)) 
        sdg@vertex.names <- paste(sdg@vertex.names)
    if (.IsEmpty(sdg@types)) 
        sdg@types <- rep("Discrete", n)
    if (!.IsEmpty(sdg@types) && (length(sdg@types) == 1)) 
        sdg@types <- rep(sdg@types, n)
    if (n != length(sdg@types)) 
        stop("Invalid data: Length of types")
    if (length(sdg@from) != length(sdg@to)) 
        stop("Invalid data: Length of from and to")
    vertex.colors <- rep(control$vertexColor, length(sdg@vertex.names))
    check.names(sdg@from, "from")
    check.names(sdg@to, "to")
    check.names(unlist(sdg@edge.list), "edge.list")
    check.names(unlist(sdg@factors), "factors")
    if (!.IsEmpty(sdg@blocks)) {
        names.blocks <- unlist(sdg@blocks)
        names.blocks <- names.blocks[grep("Vertices", names(names.blocks))]
        suppressWarnings(x <- as.numeric(names.blocks))
        if (!any(is.na(x))) 
            names.blocks <- x
        if (length(names.blocks) > 0) 
            warning("Did you give a block.tree for blocks?")
        check.names(unlist(sdg@blocks), "blocks")
    }
    if (!.IsEmpty(sdg@block.tree)) {
        names.block.tree <- unlist(sdg@block.tree)
        names.block.tree <- names.block.tree[grep("Vertices", 
            names(names.block.tree))]
        suppressWarnings(x <- as.numeric(names.block.tree))
        if (!any(is.na(x))) 
            names.block.tree <- x
        check.names(names.block.tree, "block.tree")
    }
    if (length(sdg@labels) == 0) 
        sdg@labels <- sdg@vertex.names
    if (is.null(Vertices)) 
        Vertices <- returnVertexList(sdg@vertex.names, labels = sdg@labels, 
                                     types = sdg@types, N = control$N, 
                                     colors = vertex.colors, 
                                     vertexClasses = control$vertexClasses)
    else control$N <- dim(Positions(Vertices))[2]
    X <- c("edgeList", "blockEdgeList", "factorVertexList", "factorEdgeList", 
        "extraList", "extraEdgeList")
    if (any(X %in% names(list(...)))) {
        if (!is.null(sdg)) {
            if (!.IsEmpty(sdg@from)) 
                warning("Argument 'from      ' ignored")
            if (!.IsEmpty(sdg@to)) 
                warning("Argument 'to        ' ignored")
            if (!.IsEmpty(sdg@edge.types)) 
                warning("Argument 'edge.types' ignored")
            if (!.IsEmpty(sdg@edge.list)) 
                warning("Argument 'edge.list ' ignored")
            if (!.IsEmpty(sdg@blocks)) 
                warning("Argument 'blocks    ' ignored")
            if (!.IsEmpty(sdg@block.tree)) 
                warning("Argument 'block.tree' ignored")
            if (!.IsEmpty(sdg@factors)) 
                warning("Argument 'factors   ' ignored")
            if (!.IsEmpty(sdg@texts)) 
                warning("Argument 'texts     ' ignored")
        }
        dg <- .newDgGraph(viewType = Arguments$viewType, 
                          vertexList = Arguments$vertexList, 
                          edgeList = Arguments$vertexEdges, 
                          oriented = Arguments$oriented, 
                          blockList = Arguments$blockList, 
                          blockEdgeList = Arguments$blockEdges, 
                          factorVertexList = Arguments$factorVertices, 
                          factorEdgeList = Arguments$factorEdges, 
                          extraList = Arguments$extraVertices, 
                          extraEdgeList = Arguments$extraEdges)
    }
    else {
        BlockList <- NULL
        BlockTree <- NULL
        if (!is.null(frameModels)) {
            BlockList <- frameModels@blocks
            if (!is.null(BlockList) && ((length(BlockList) == 
                0) || is.null(BlockList[[1]]))) 
                BlockList <- NULL
        }
        else if (!is.null(Arguments$BlockList)) {
            BlockList <- Arguments$BlockList
            if (is.null(Arguments$Vertices) && is.null(Arguments$vertexList)) 
                warning(
    "Argument 'vertexList' should also be given to put verices in blocks!")
        }
        else if (!is.null(Arguments$blockList)) {
            BlockList <- Arguments$blockList
            if (is.null(Arguments$Vertices) && is.null(Arguments$vertexList)) 
                warning(
    "Argument 'vertexList' should also be given to put verices in blocks!")
        }
        else if (!(.IsEmpty(sdg@blocks))) {
            result <- setBlocks(sdg@blocks, Vertices, 
                                right.to.left = control$right.to.left, 
                                nested.blocks = control$nested.blocks, 
                                blockColors = control$blockColors, 
                                N = control$N)
            if (is.null(Arguments$Vertices)) 
                Vertices <- result$Vertices
            if (control$drawblocks) 
                BlockList <- result$Blocks
        }
        else if (!(.IsEmpty(sdg@block.tree))) {
            result <- setTreeBlocks(sdg@block.tree, Vertices, 
                                    root.label = "ROOT", N = control$N, 
                                    blockColors = control$blockColors, 
                                    overlaying = control$overlaying)
            if (is.null(Arguments$Vertices)) 
                Vertices <- result$Vertices
            if (control$drawblocks) 
                BlockTree <- result$BlockTree
        }
        if (!is.null(dg)) {
            ExtraVertices = dg@extraList
            Edges = dg@edgeList
            BlockEdges = dg@blockEdgeList
            FactorVertices = dg@factorVertexList
            FactorEdges = dg@factorEdgeList
            ExtraEdges = dg@extraEdgeList
            if (length(ExtraVertices) > 0) 
                if (is.null(ExtraVertices[[1]])) 
                  ExtraVertices <- NULL
            if (length(Edges) > 0) 
                if (is.null(Edges[[1]])) 
                  Edges <- NULL
            if (length(BlockEdges) > 0) 
                if (is.null(BlockEdges[[1]])) 
                  BlockEdges <- NULL
            if (length(FactorVertices) > 0) 
                if (is.null(FactorVertices[[1]])) 
                  FactorVertices <- NULL
            if (length(FactorEdges) > 0) 
                if (is.null(FactorEdges[[1]])) 
                  FactorEdges <- NULL
            if (length(ExtraEdges) > 0) 
                if (is.null(ExtraEdges[[1]])) 
                  ExtraEdges <- NULL
        }
        else {
            if (!is.null(Arguments$ExtraVertices)) 
                ExtraVertices <- Arguments$ExtraVertices
            else if (.IsEmpty(sdg@texts)) 
                ExtraVertices <- NULL
            else ExtraVertices <- returnVertexList(
                paste("T", 1:length(sdg@texts), sep = ""), 
                labels = sdg@texts, 
                types = rep("TextVertex", length(sdg@texts)), 
                line = TRUE, N = control$N, 
                colors = rep(control$extraVertexColor, length(sdg@texts)), 
                vertexClasses = control$vertexClasses)
                # 'vertexClasses' = control$vertexClasses ???
            if (!is.null(ExtraVertices) && !is.null(control$diagonal) && 
                control$diagonal) {
                i <- 1:length(ExtraVertices) - 2
                Positions(ExtraVertices) <- (matrix(c(rep(i * 50, 2), 
                  rep(200, (control$N - 2) * length(ExtraVertices))), 
                  ncol = control$N) - 200)/4
            }
            if (!is.null(Arguments$ExtraEdges)) 
                ExtraEdges <- Arguments$ExtraEdges
            else {
                ExtraEdges <- NULL
                if (.IsEmpty(sdg@extra.edge.list)) 
                  sdg@extra.edge.list <- two.to.pairs(sdg@extra.from, 
                    sdg@extra.to)
                if (!(.IsEmpty(sdg@extra.edge.list))) 
                  ExtraEdges <- returnExtraEdgeList(sdg@extra.edge.list, 
                    Vertices, ExtraVertices, color = control$extraEdgeColor)
            }
            FactorVertices <- NULL
            FactorEdges <- NULL
            if (!is.null(Arguments$FactorVertices)) {
                FactorVertices <- Arguments$FactorVertices
                if (is.null(Arguments$FactorEdges)) 
                  warning(
    "Argument FactorEdges should also be given with FactorVertices!")
            }
            if (!(.IsEmpty(sdg@factors))) {
                result <- returnFactorVerticesAndEdges(Vertices, 
                          sdg@factors, 
                          factorVertexColor = control$factorVertexColor, 
                          factorEdgeColor = control$factorEdgeColor, 
                          fixedFactorPositions = control$fixedFactorPositions, 
                          factorClasses = control$factorClasses)
                if (!is.null(Arguments$FactorVertices)) {
                  FactorVertices <- Arguments$FactorVertices
                }
                else FactorVertices <- result$FactorVertices
                if (!is.null(Arguments$FactorEdges)) 
                  FactorEdges <- Arguments$FactorEdges
                else FactorEdges <- result$FactorEdges
                if ((.IsEmpty(sdg@from))) {
                  sdg@from <- result$PairEdges[, 1]
                  sdg@to <- result$PairEdges[, 2]
                }
            }
            Oriented <- sdg@oriented
            if (is.list(Oriented)) 
                Oriented <- unlist(Oriented)
            if (.IsEmpty(Oriented)) 
                Oriented <- NA
            if (length(Oriented) == 0) 
                Oriented <- NA
            if (.IsEmpty(sdg@edge.list)) 
                sdg@edge.list <- two.to.pairs(sdg@from, sdg@to)
            if (!is.null(Arguments$Edges)) 
                Edges <- Arguments$Edges
            else Edges <- returnEdgeList(sdg@edge.list, Vertices, 
                color = control$edgeColor, oriented = Oriented, 
                N = control$N, types = sdg@edge.types, 
                edgeClasses = control$edgeClasses)
            if (length(Oriented) > 1) 
                Oriented <- TRUE
            BlockEdges <- NULL
            if (!is.null(Arguments$BlockEdges)) 
                BlockEdges <- Arguments$BlockEdges
            else if (TRUE && (!is.null(BlockList) || !is.null(BlockTree))) {
                if (!(.IsEmpty(sdg@factors))) 
                  message("Edges between blocks and factors not implemented!")
                if (is.null(BlockList) && !is.null(BlockTree)) 
                  BlockList <- blockTreeToList(BlockTree)
                if (.IsEmpty(BlockList)) 
                  BlockEdges <- NULL
                else BlockEdges <- returnBlockEdgeList(sdg@edge.list, 
                  Vertices, BlockList, color = control$blockEdgeColor, 
                  oriented = Oriented)
            }
        }
        dg <- .newDgGraph(viewType = sdg@viewType, vertexList = Vertices, 
                          edgeList = Edges, oriented = Oriented, 
                          blockList = BlockList, blockEdgeList = BlockEdges, 
                          factorVertexList = FactorVertices, 
                          factorEdgeList = FactorEdges, 
                          extraList = ExtraVertices, 
                          extraEdgeList = ExtraEdges)
    }
    return(dg)
}
