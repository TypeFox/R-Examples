    setClass("dg.list", contains = "list", 
        representation("list"), prototype = list())

    setClass("dg.NodeList", 
        contains = c("dg.list", "list"))
    setClass("dg.VertexList", 
        contains = c("dg.NodeList", "dg.list", "list"))
    setClass("dg.BlockList", 
        contains = c("dg.NodeList", "dg.list", "list"))
    setClass("dg.FactorVertexList", 
        contains = c("dg.NodeList", "dg.list", "list"))

    setClass("dg.EdgeList", 
        contains = c("dg.NodeList", "dg.list", "list"))
    setClass("dg.VertexEdgeList", 
        contains = c("dg.EdgeList", "dg.NodeList", "dg.list", "list"))
    setClass("dg.BlockEdgeList", 
        contains = c("dg.EdgeList", "dg.NodeList", "dg.list", "list"))
    setClass("dg.FactorEdgeList", 
        contains = c("dg.EdgeList", "dg.NodeList", "dg.list", "list"))
    setClass("dg.ExtraEdgeList", 
        contains = c("dg.EdgeList", "dg.NodeList", "dg.list", "list"))

    setMethod("initialize", "dg.VertexList", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("names", names(Args)) || is.element("types", 
                names(Args))) {
                names <- NULL
                labels <- NULL
                types <- NULL
                strata <- NULL
                line <- FALSE
                N <- 3
                colors <- ifelse(types == "TextVertex", 
                                 "FloralWhite", "DarkRed")
                vertexClasses <- validVertexClasses()
                if (is.element("names", names(Args))) 
                  names <- Args$names
                if (is.element("labels", names(Args))) 
                  labels <- Args$labels
                if (is.element("types", names(Args))) 
                  types <- Args$types
                if (is.element("strata", names(Args))) 
                  strata <- Args$strata
                if (is.element("line", names(Args))) 
                  line <- Args$line
                if (is.element("N", names(Args))) 
                  N <- Args$N
                if (is.element("colors", names(Args))) 
                  colors <- Args$colors
                if (is.element("vertexClasses", names(Args))) 
                  vertexClasses <- Args$vertexClasses
                result <- returnVertexList(names = names, labels = labels, 
                  types = types, strata = strata, line = line, 
                  N = N, colors = colors, vertexClasses = vertexClasses)
                .Object@.Data <- result
            }
            else .Object@.Data <- checkList("dg.VertexList", 
                "dg.Vertex", Args)
        }
        return(.Object)
    })
    setMethod("initialize", "dg.VertexEdgeList", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("edge.list", names(Args))) {
                edge.list <- NULL
                vertices <- NULL
                width <- 2
                color <- "DarkSlateGrey"
                N <- 3
                oriented <- NA
                types <- NULL
                edgeClasses <- validEdgeClasses()
                if (is.element("edge.list", names(Args))) 
                  edge.list <- Args$edge.list
                if (is.element("vertices", names(Args))) 
                  vertices <- Args$vertices
                if (is.element("width", names(Args))) 
                  width <- Args$width
                if (is.element("colors", names(Args))) 
                  color <- Args$color
                if (is.element("N", names(Args))) 
                  N <- Args$N
                if (is.element("oriented", names(Args))) 
                  oriented <- Args$oriented
                if (is.element("types", names(Args))) 
                  types <- Args$types
                if (is.element("edgeClasses", names(Args))) 
                  edgeClasses <- Args$edgeClasses
                result <- returnEdgeList(edge.list = edge.list, 
                  vertices = vertices, width = width, color = color, 
                  N = N, oriented = oriented, types = types, 
                  edgeClasses = edgeClasses)
                .Object@.Data <- result
            }
            else .Object@.Data <- checkList("dg.VertexEdgeList", 
                "dg.VertexEdge", Args)
        }
        return(.Object)
    })
    setMethod("initialize", "dg.ExtraEdgeList", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("edge.list", names(Args))) {
                edge.list <- NULL
                vertices <- NULL
                extravertices <- NULL
                width <- 2
                color <- "DarkSlateGrey"
                N <- 3
                type <- NULL
                if (is.element("edge.list", names(Args))) 
                  edge.list <- Args$edge.list
                if (is.element("vertices", names(Args))) 
                  vertices <- Args$vertices
                if (is.element("extravertices", names(Args))) 
                  extravertices <- Args$extravertices
                if (is.element("width", names(Args))) 
                  width <- Args$width
                if (is.element("colors", names(Args))) 
                  color <- Args$color
                if (is.element("N", names(Args))) 
                  N <- Args$N
                if (is.element("type", names(Args))) 
                  type <- Args$type
                result <- returnExtraEdgeList(edge.list = edge.list, 
                  vertices = vertices, extravertices = extravertices, 
                  width = width, color = color, N = N, type = type)
                .Object@.Data <- result
            }
            else .Object@.Data <- checkList("dg.ExtraEdgeList", 
                "dg.ExtraEdge", Args)
        }
        return(.Object)
    })
    setMethod("initialize", "dg.BlockEdgeList", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("edge.list", names(Args))) {
                edge.list <- NULL
                vertices <- NULL
                blocks <- NULL
                visibleBlocks <- 1:length(blocks)
                width <- 2
                color <- "default"
                N <- 3
                oriented <- NA
                type <- NULL
                if (is.element("edge.list", names(Args))) 
                  edge.list <- Args$edge.list
                if (is.element("vertices", names(Args))) 
                  vertices <- Args$vertices
                if (is.element("blocks", names(Args))) 
                  blocks <- Args$blocks
                if (is.element("visibleBlocks", names(Args))) 
                  visibleBlocks <- Args$visibleBlocks
                if (is.element("width", names(Args))) 
                  width <- Args$width
                if (is.element("colors", names(Args))) 
                  color <- Args$color
                if (is.element("N", names(Args))) 
                  N <- Args$N
                if (is.element("oriented", names(Args))) 
                  oriented <- Args$oriented
                if (is.element("type", names(Args))) 
                  type <- Args$type
                result <- returnBlockEdgeList(edge.list = edge.list, 
                  vertices = vertices, blocks = blocks, 
                  visibleBlocks = visibleBlocks, 
                  width = width, color = color, N = N, 
                  oriented = oriented, 
                  type = type)
                .Object@.Data <- result
            }
            else .Object@.Data <- checkList("dg.BlockEdgeList", 
                "dg.BlockEdge", Args)
        }
        return(.Object)
    })
    setMethod("initialize", "dg.FactorEdgeList", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("edge.list", names(Args))) {
                edge.list <- NULL
                vertices <- NULL
                factorvertices <- NULL
                width <- 2
                color <- "DarkSlateGrey"
                N <- 3
                type <- NULL
                if (is.element("edge.list", names(Args))) 
                  edge.list <- Args$edge.list
                if (is.element("vertices", names(Args))) 
                  vertices <- Args$vertices
                if (is.element("factorvertices", names(Args))) 
                  factorvertices <- Args$factorvertices
                if (is.element("width", names(Args))) 
                  width <- Args$width
                if (is.element("colors", names(Args))) 
                  color <- Args$color
                if (is.element("N", names(Args))) 
                  N <- Args$N
                if (is.element("type", names(Args))) 
                  type <- Args$type
                result <- returnFactorEdgeList(edge.list = edge.list, 
                  vertices = vertices, factorvertices = factorvertices, 
                  width = width, color = color, N = N, type = type)
                .Object@.Data <- result
            }
            else .Object@.Data <- checkList("dg.FactorEdgeList", 
                "dg.FactorEdge", Args)
        }
        return(.Object)
    })
    if (!isGeneric("descendantsBlockList")) {
        if (is.function("descendantsBlockList")) 
            fun <- descendantsBlockList
        else fun <- function(blockList, index = NULL) 
                                    standardGeneric("descendantsBlockList")
        setGeneric("descendantsBlockList", fun)
    }
    setMethod("descendantsBlockList", "dg.BlockList", 
        function(blockList, 
        index = NULL) {
        if (is.null(index)) {
            result <- lapply(0:length(blockList), 
                              function(i) descendantsBlockList(blockList, i))
            names(result) <- c("root", Names(blockList))
        }
        else {
            result <- NULL
            if (!.IsEmpty(blockList)) 
                if (index > 0) {
                  for (j in children(blockList[[index]])) if (j > 0) 
                    result <- unique(sort(c(result, j, 
                                  descendantsBlockList(blockList, index = j))))
                }
                else for (j in 1:length(blockList)) if (index == 
                  parent(blockList[[j]])) 
                  result <- unique(sort(c(result, j, 
                                descendantsBlockList(blockList, index = j))))
        }
        return(result)
    })
    if (!isGeneric("ancestorsBlockList")) {
        if (is.function("ancestorsBlockList")) 
            fun <- ancestorsBlockList
        else fun <- function(blockList, index = NULL) 
                              standardGeneric("ancestorsBlockList")
        setGeneric("ancestorsBlockList", fun)
    }
    setMethod("ancestorsBlockList", "dg.BlockList", 
        function(blockList, index = NULL) {
        if (is.null(index)) {
            result <- lapply(1:length(blockList), function(i) 
                                             ancestorsBlockList(blockList, i))
            names(result) <- c(Names(blockList))
        }
        else {
            result <- NULL
            if ((index > 0) && !.IsEmpty(blockList)) {
                parent <- parent(blockList[[index]])
                if (parent > 0) 
                  result <- unique(sort(c(result, parent, 
                    ancestorsBlockList(blockList, index = parent))))
            }
        }
        return(result)
    })
    if (!isGeneric("checkBlockList")) {
        if (is.function("checkBlockList")) 
            fun <- checkBlockList
        else fun <- function(blockList) standardGeneric("checkBlockList")
        setGeneric("checkBlockList", fun)
    }
    setMethod("checkBlockList", "dg.BlockList", 
        function(blockList) {
        for (j in seq(along = blockList)) {
            parent <- parent(blockList[[j]])
            if (parent > 0) {
                children <- children(blockList[[parent]])
                if (!is.element(j, children)) {
                  message(paste("Block ", j, ", '", name(blockList[[j]]), 
                    "', inserted as child of block ", parent, 
                    ", '", name(blockList[[parent]]), "'.", sep = ""))
                  z <- unique(c(children, parent))
                  children(blockList[[parent]]) <- z[z != 0]
                }
            }
        }
        for (j in seq(along = blockList)) {
            children <- children(blockList[[j]])
            for (k in children) if (k > 0) {
                parent <- parent(blockList[[k]])
                if (parent != j) {
                  message(paste("Block ", k, ", '", name(blockList[[k]]), 
                    "', deleted as child of block ", j, ", '", 
                    name(blockList[[j]]), "' .", sep = ""))
                  children <- children[children != k]
                  if (is.null(children)) 
                    children <- numeric(0)
                  children(blockList[[j]]) <- children
                }
            }
        }
        for (j in seq(along = blockList)) {
            x <- descendantsBlockList(blockList, j)
            y <- blockList[[j]]@descendants
            if (!setequal(x, y) && !(is.null(x) && (length(y) == 1) 
                && (y == 0))) {
                message(paste("Descendants of block", j, ",", 
                  name(blockList[[j]]), ",", "changed from", 
                  paste(y, collapse = ","), "to", paste(x, collapse = ","), 
                  "."))
                if (is.null(x)) 
                  x <- numeric(0)
                blockList[[j]]@descendants <- x
            }
            x <- ancestorsBlockList(blockList, j)
            y <- blockList[[j]]@ancestors
            if (!setequal(x, y) && !(is.null(x) && 
                (length(y) == 1) && (y == 0))) {
                message(paste("Ancestors of block", j, ",", 
                  name(blockList[[j]]), 
                  ",", "changed from", paste(y, collapse = ","), 
                  "to", paste(x, collapse = ","), "."))
                if (is.null(x)) 
                  x <- numeric(0)
                blockList[[j]]@ancestors <- x
            }
        }
        return(blockList)
    })
    setClass("dg.graphedges", 
        representation(viewType = "character", 
                       visibleVertices = "numeric", visibleBlocks = "numeric", 
                       oriented = "logical", edgeList = "dg.VertexEdgeList", 
                       blockEdgeList = "dg.BlockEdgeList", 
                       factorVertexList = "dg.FactorVertexList", 
                       factorEdgeList = "dg.FactorEdgeList", 
                       extraList = "dg.VertexList", 
                       extraEdgeList = "dg.ExtraEdgeList"), 
        prototype = list(viewType = "Simple", 
                         visibleVertices = numeric(), 
                         visibleBlocks = numeric(), 
                         edgeList = new("dg.VertexEdgeList"), 
                         oriented = NA, 
                         blockEdgeList = new("dg.BlockEdgeList"), 
                         factorVertexList = new("dg.FactorVertexList"), 
                         factorEdgeList = new("dg.FactorEdgeList"), 
                         extraList = new("dg.VertexList"), 
                         extraEdgeList = new("dg.ExtraEdgeList")))
    if (!isGeneric("visibleVertices")) {
        if (is.function("visibleVertices")) 
            fun <- visibleVertices
        else fun <- function(object) standardGeneric("visibleVertices")
        setGeneric("visibleVertices", fun)
    }
    setMethod("visibleVertices", "dg.graphedges", 
        function(object) object@visibleVertices)
    setGeneric("visibleVertices<-", 
        function(x, value) standardGeneric("visibleVertices<-"))
    setReplaceMethod("visibleVertices", "dg.graphedges", 
        function(x, value) {
        x@visibleVertices <- value
        x
    })
    if (!isGeneric("visibleBlocks")) {
        if (is.function("visibleBlocks")) 
            fun <- visibleBlocks
        else fun <- function(object) standardGeneric("visibleBlocks")
        setGeneric("visibleBlocks", fun)
    }
    setMethod("visibleBlocks", "dg.graphedges", 
        function(object) object@visibleBlocks)
    setGeneric("visibleBlocks<-", 
        function(x, value) standardGeneric("visibleBlocks<-"))
    setReplaceMethod("visibleBlocks", "dg.graphedges", 
        function(x, value) {
        x@visibleBlocks <- value
        x
    })
    if (!isGeneric("edgeList")) {
        if (is.function("edgeList")) 
            fun <- edgeList
        else fun <- function(object) standardGeneric("edgeList")
        setGeneric("edgeList", fun)
    }
    setMethod("edgeList", "dg.graphedges", 
        function(object) object@edgeList)
    setGeneric("edgeList<-", 
        function(x, value) standardGeneric("edgeList<-"))
    setReplaceMethod("edgeList", "dg.graphedges", 
        function(x, value) {
        x@edgeList <- value
        x
    })
    if (!isGeneric("viewType")) {
        if (is.function("viewType")) 
            fun <- viewType
        else fun <- function(object) standardGeneric("viewType")
        setGeneric("viewType", fun)
    }
    setMethod("viewType", "dg.graphedges", 
        function(object) object@viewType)
    setGeneric("viewType<-", 
        function(x, value) standardGeneric("viewType<-"))
    setReplaceMethod("viewType", "dg.graphedges", 
        function(x, value) {
        x@viewType <- value
        x
    })
    if (!isGeneric("blockList")) {
        if (is.function("blockList")) 
            fun <- blockList
        else fun <- function(object) standardGeneric("blockList")
        setGeneric("blockList", fun)
    }
    setMethod("blockList", "dg.graphedges", 
        function(object) object@blockList)
    setGeneric("blockList<-", 
        function(x, value) standardGeneric("blockList<-"))
    setReplaceMethod("blockList", "dg.graphedges", 
        function(x, value) {
        x@blockList <- value
        x
    })
    if (!isGeneric("blockEdgeList")) {
        if (is.function("blockEdgeList")) 
            fun <- blockEdgeList
        else fun <- function(object) standardGeneric("blockEdgeList")
        setGeneric("blockEdgeList", fun)
    }
    setMethod("blockEdgeList", "dg.graphedges", 
        function(object) object@blockEdgeList)
    setGeneric("blockEdgeList<-", 
        function(x, value) standardGeneric("blockEdgeList<-"))
    setReplaceMethod("blockEdgeList", "dg.graphedges", 
        function(x, value) {
        x@blockEdgeList <- value
        x
    })
    if (!isGeneric("factorVertexList")) {
        if (is.function("factorVertexList")) 
            fun <- factorVertexList
        else fun <- function(object) standardGeneric("factorVertexList")
        setGeneric("factorVertexList", fun)
    }
    setMethod("factorVertexList", "dg.graphedges", 
        function(object) object@factorVertexList)
    setGeneric("factorVertexList<-", 
        function(x, value) standardGeneric("factorVertexList<-"))
    setReplaceMethod("factorVertexList", "dg.graphedges", 
        function(x, value) {
        x@factorVertexList <- value
        x
    })
    if (!isGeneric("factorEdgeList")) {
        if (is.function("factorEdgeList")) 
            fun <- factorEdgeList
        else fun <- function(object) standardGeneric("factorEdgeList")
        setGeneric("factorEdgeList", fun)
    }
    setMethod("factorEdgeList", "dg.graphedges", 
        function(object) object@factorEdgeList)
    setGeneric("factorEdgeList<-", 
        function(x, value) standardGeneric("factorEdgeList<-"))
    setReplaceMethod("factorEdgeList", "dg.graphedges", 
        function(x, value) {
        x@factorEdgeList <- value
        x
    })
    if (!isGeneric("extraList")) {
        if (is.function("extraList")) 
            fun <- extraList
        else fun <- function(object) standardGeneric("extraList")
        setGeneric("extraList", fun)
    }
    setMethod("extraList", "dg.graphedges", 
        function(object) object@extraList)
    setGeneric("extraList<-", 
        function(x, value) standardGeneric("extraList<-"))
    setReplaceMethod("extraList", "dg.graphedges", 
        function(x, value) {
        x@extraList <- value
        x
    })
    if (!isGeneric("extraEdgeList")) {
        if (is.function("extraEdgeList")) 
            fun <- extraEdgeList
        else fun <- function(object) standardGeneric("extraEdgeList")
        setGeneric("extraEdgeList", fun)
    }
    setMethod("extraEdgeList", "dg.graphedges", 
        function(object) object@extraEdgeList)
    setGeneric("extraEdgeList<-", 
        function(x, value) standardGeneric("extraEdgeList<-"))
    setReplaceMethod("extraEdgeList", "dg.graphedges", 
        function(x, value) {
        x@extraEdgeList <- value
        x
    })
    setClass("dg.graph", contains = "dg.graphedges", 
        representation(vertexList = "dg.VertexList", 
                       blockList = "dg.BlockList"), 
        prototype = list(viewType = "Simple", 
                         vertexList = new("dg.VertexList"), 
                         visibleVertices = numeric(), 
                         visibleBlocks = numeric(), 
                         edgeList = new("dg.VertexEdgeList"), 
                         oriented = NA, blockList = new("dg.BlockList"), 
                         blockEdgeList = new("dg.BlockEdgeList"), 
                         factorVertexList = new("dg.FactorVertexList"), 
                         factorEdgeList = new("dg.FactorEdgeList"), 
                         extraList = new("dg.VertexList"), 
                         extraEdgeList = new("dg.ExtraEdgeList")))
    setClass("dg.simple.graph", 
        representation(viewType = "character", 
                       vertex.names = "vector", 
                       types = "character", labels = "vector", 
                       from = "vector", to = "vector", 
                       edge.list = "list", edge.types = "character", 
                       blocks = "list", block.tree = "list", 
                       oriented = "vector", factors = "list", 
                       texts = "character", extra.from = "vector", 
                       extra.to = "vector", extra.edge.list = "list"), 
        prototype = list(viewType = "Simple", vertex.names = vector(), 
                         labels = vector(), types = "", from = vector(), 
                         to = vector(), edge.list = list(NULL), 
                         edge.types = "", blocks = list(NULL), 
                         block.tree = list(NULL), 
                         oriented = NA, factors = list(NULL), 
                         texts = "", extra.from = vector(), 
                         extra.to = vector(), extra.edge.list = list(NULL)))
    setAs("dg.simple.graph", "dg.graph", 
        function(from, to) simpleGraphToGraph(from))
    if (!isGeneric("dg")) {
        if (is.function("dg")) 
            fun <- position
        else fun <- function(object, ...) standardGeneric("dg")
        setGeneric("dg", fun)
    }
    setMethod("dg", "dg.graph", 
        function(object, ...) {
        dots <- list(...)
        modelObject <- dots$modelObject
        modelObjectName <- dots$modelObjectName
        dynamicGraphMain(vertexList = object@vertexList, 
            blockList = object@blockList, dg = object, 
            object = modelObject, objectName = modelObjectName, ...)
    })
    setMethod("dg", "dg.simple.graph", 
        function(object, ...) {
        dg <- simpleGraphToGraph(sdg = object, ...)
        dg(dg, ...)
    })
    if (!isGeneric("addModel")) {
        if (is.function("addModel")) 
            fun <- position
        else fun <- function(object, frameModels, ...) 
            standardGeneric("addModel")
        setGeneric("addModel", fun)
    }
    setMethod("addModel", "dg.graph", 
        function(object, frameModels, ...) 
        .addModel(dg = object, frameModels = frameModels, 
                  overwrite = FALSE, ...))
    setMethod("addModel", "dg.graphedges", 
        function(object, frameModels, ...) 
        .addModel(dg = object, frameModels = frameModels, 
                  overwrite = FALSE, ...))
    setMethod("addModel", "dg.simple.graph", 
        function(object, frameModels, ...) 
        .addModel(dg = object, frameModels = frameModels, 
                  overwrite = FALSE, ...))
    if (!isGeneric("addView")) {
        if (is.function("addView")) 
            fun <- position
        else fun <- function(object, frameModels, 
                             frameViews = frameModels@models[[modelIndex]], 
                             modelIndex = 1, ...) 
            standardGeneric("addView")
        setGeneric("addView", fun)
    }
    setMethod("addView", "dg.graph", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, ...) 
        .addView(dg = object, frameModels = frameModels, 
                 frameViews = frameViews, graphWindow = NULL, 
                 overwrite = FALSE, ...))
    setMethod("addView", "dg.graphedges", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, ...) 
        .addView(dg = object, frameModels = frameModels, 
                 frameViews = frameViews, graphWindow = NULL, 
                 overwrite = FALSE, ...))
    setMethod("addView", "dg.simple.graph", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, ...) 
        .addView(dg = object, frameModels = frameModels, 
                 frameViews = frameViews, graphWindow = NULL, 
                 overwrite = FALSE, ...))
    if (!isGeneric("replaceModel")) {
        if (is.function("replaceModel")) 
            fun <- position
        else fun <- function(object, frameModels, 
                             frameViews = frameModels@models[[modelIndex]], 
                             modelIndex = 1, 
                             graphWindow = frameViews@graphs[[viewIndex]], 
                             viewIndex = 1, ...) 
                                 standardGeneric("replaceModel")
        setGeneric("replaceModel", fun)
    }
    setMethod("replaceModel", "dg.graph", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, 
                 graphWindow = frameViews@graphs[[viewIndex]],
                 viewIndex = 1, ...) 
        .addModel(dg = object, frameModels = frameModels, 
                  frameViews = frameViews, graphWindow = graphWindow, 
                  overwrite = TRUE, ...))
    setMethod("replaceModel", "dg.graphedges", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
                 viewIndex = 1, ...) 
        .addModel(dg = object, frameModels = frameModels, 
                  frameViews = frameViews, graphWindow = graphWindow, 
                  overwrite = TRUE, ...))
    setMethod("replaceModel", "dg.simple.graph", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
                 viewIndex = 1, control = dg.control(...), ...) 
        .addModel(dg = object, 
                  frameModels = frameModels, frameViews = frameViews, 
                  graphWindow = graphWindow, overwrite = TRUE, ...))
    if (!isGeneric("replaceView")) {
        if (is.function("replaceView")) 
            fun <- position
        else fun <- function(object, frameModels, 
                             frameViews = frameModels@models[[modelIndex]], 
                             modelIndex = 1, 
                             graphWindow = frameViews@graphs[[viewIndex]], 
                             viewIndex = 1, ...) standardGeneric("replaceView")
        setGeneric("replaceView", fun)
    }
    setMethod("replaceView", "dg.graph", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, 
                 graphWindow = frameViews@graphs[[viewIndex]], 
                 viewIndex = 1, ...) 
        .addView(dg = object, frameModels = frameModels, 
                 frameViews = frameViews, graphWindow = graphWindow, 
                 overwrite = TRUE, ...))
    setMethod("replaceView", "dg.graphedges", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
                 viewIndex = 1, ...) 
        .addView(dg = object, frameModels = frameModels, 
                 frameViews = frameViews, graphWindow = graphWindow, 
                 overwrite = TRUE, ...))
    setMethod("replaceView", "dg.simple.graph", 
        function(object, frameModels, 
                 frameViews = frameModels@models[[modelIndex]], 
                 modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
                 viewIndex = 1, control = dg.control(...), ...) 
        .addView(dg = object, frameModels = frameModels, 
                 frameViews = frameViews, 
                 graphWindow = graphWindow, overwrite = TRUE, ...))
    setClass("DynamicGraph", 
             representation(id.env = "character", label = "character", 
                            vertices = "dg.VertexList", 
                            blocks = "dg.BlockList", 
                            control = "list", models = "list"))
    if (!isGeneric("vertices")) {
        if (is.function("vertices")) 
            fun <- vertices
        else fun <- function(object) standardGeneric("vertices")
        setGeneric("vertices", fun)
    }
    setMethod("vertices", "DynamicGraph", 
        function(object) object@vertices)
    setGeneric("vertices<-", 
        function(x, value) standardGeneric("vertices<-"))
    setReplaceMethod("vertices", "DynamicGraph", 
        function(x, value) {
        x@vertices <- value
        replaceVertexList(value, x)
        x
    })
    if (!isGeneric("blocks")) {
        if (is.function("blocks")) 
            fun <- blocks
        else fun <- function(object) standardGeneric("blocks")
        setGeneric("blocks", fun)
    }
    setMethod("blocks", "DynamicGraph", 
        function(object) object@blocks)
    setGeneric("blocks<-", 
        function(x, value) standardGeneric("blocks<-"))
    setReplaceMethod("blocks", "DynamicGraph", 
        function(x, value) {
        x@blocks <- value
        replaceBlockList(value, x)
        x
    })
    if (!isGeneric("control")) {
        if (is.function("control")) 
            fun <- control
        else fun <- function(object) standardGeneric("control")
        setGeneric("control", fun)
    }
    setMethod("control", "DynamicGraph", 
        function(object) object@control)
    setGeneric("control<-", 
        function(x, value) standardGeneric("control<-"))
    setReplaceMethod("control", "DynamicGraph", 
        function(x, value) {
        x@control <- value
        replaceControls(value, x)
        x
    })
    if (!isGeneric("models")) {
        if (is.function("models")) 
            fun <- models
        else fun <- function(object) standardGeneric("models")
        setGeneric("models", fun)
    }
    setMethod("models", "DynamicGraph", 
        function(object) object@models)
    setGeneric("models<-", 
        function(x, value) standardGeneric("models<-"))
    setReplaceMethod("models", "DynamicGraph", 
        function(x, value) {
        x@models <- value
        x
    })
    setMethod("dg", "DynamicGraph", 
        function(object, ...) {
        dots <- list(...)
        modelObject <- dots$modelObject
        modelObjectName <- dots$modelObjectName
        control <- dots$control
        if (is.null(control)) {
            control <- object@control
            if (is.null(control)) 
                control <- dg.control()
            dynamicGraphMain(vertexList = object@vertices, 
                             blockList = object@blocks, 
                             object = modelObject, 
                             objectName = modelObjectName, 
                             frameModels = object, redraw = TRUE, 
                             control = control, ...)
        }
        else dynamicGraphMain(vertexList = object@vertices, 
                              blockList = object@blocks, 
                              object = modelObject, 
                              objectName = modelObjectName, 
                              frameModels = object, redraw = TRUE, ...)
    })
    setClass("DynamicGraphModel", representation(id.env = "character", 
        label = "character", index = "numeric", model = "list", 
        graphs = "list"))
    if (!isGeneric("model")) {
        if (is.function("model")) 
            fun <- model
        else fun <- function(object) standardGeneric("model")
        setGeneric("model", fun)
    }
    setMethod("model", "DynamicGraphModel", 
        function(object) object@model)
    setGeneric("model<-", 
        function(x, value) standardGeneric("model<-"))
    setReplaceMethod("model", "DynamicGraphModel", 
        function(x, value) {
        x@model <- value
        x
    })
    if (!isGeneric("graphs")) {
        if (is.function("graphs")) 
            fun <- graphs
        else fun <- function(object) standardGeneric("graphs")
        setGeneric("graphs", fun)
    }
    setMethod("graphs", "DynamicGraphModel", 
        function(object) object@graphs)
    setGeneric("graphs<-", 
        function(x, value) standardGeneric("graphs<-"))
    setReplaceMethod("graphs", "DynamicGraphModel", 
        function(x, value) {
        x@graphs <- value
        x
    })
    setReplaceMethod("control", "DynamicGraphModel", 
        function(x, value) {
        replaceControls(value, frameViews = x)
        x
    })
    setClass("DynamicGraphView", representation(id.env = "character", 
        id = "numeric", label = "character", index = "numeric", 
        dg = "dg.graphedges"))
    if (!isGeneric("label")) {
        if (is.function("label")) 
            fun <- label
        else fun <- function(object) standardGeneric("label")
        setGeneric("label", fun)
    }
    setMethod("label", "DynamicGraph", 
        function(object) object@label)
    setMethod("label", "DynamicGraphModel", 
        function(object) object@label)
    setMethod("label", "DynamicGraphView", 
        function(object) object@label)
    setGeneric("label<-", 
        function(x, value) standardGeneric("label<-"))
    setReplaceMethod("label", "DynamicGraph", 
        function(x, value) {
        x@label <- value
        x
    })
    setReplaceMethod("label", "DynamicGraphModel", 
        function(x, value) {
        x@label <- value
        x
    })
    setReplaceMethod("label", "DynamicGraphView", 
        function(x, value) {
        x@label <- value
        x
    })
    setMethod("dg", "DynamicGraphView", 
        function(object, ...) object@dg)
    setGeneric("dg<-", 
        function(x, value) standardGeneric("dg<-"))
    setReplaceMethod("dg", "DynamicGraphView", 
        function(x, value) {
        x@dg <- value
        x
    })
    setReplaceMethod("control", "DynamicGraphView", 
        function(x, value) {
        replaceControls(value, graphWindow = x)
        x
    })
    if (!isGeneric("top")) {
        if (is.function("top")) 
            fun <- top
        else fun <- function(object) standardGeneric("top")
        setGeneric("top", fun)
    }
    setMethod("top", "DynamicGraphView", 
        function(object) .get.env.graphWindow(graphWindow = object)$env$top)
    if (!isGeneric("vbox")) {
        if (is.function("vbox")) 
            fun <- vbox
        else fun <- function(object) standardGeneric("vbox")
        setGeneric("vbox", fun)
    }
    setMethod("vbox", "DynamicGraphView", 
        function(object) 
            .get.env.graphWindow(graphWindow = object)$env$top$env$box)
    if (!isGeneric("canvas")) {
        if (is.function("canvas")) 
            fun <- canvas
        else fun <- function(object) standardGeneric("canvas")
        setGeneric("canvas", fun)
    }
    setMethod("canvas", "DynamicGraphView", 
        function(object) 
            .get.env.graphWindow(graphWindow = object)$env$top$env$canvas)
    if (!isGeneric("viewLabel")) {
        if (is.function("viewLabel")) 
            fun <- viewLabel
        else fun <- function(object) standardGeneric("viewLabel")
        setGeneric("viewLabel", fun)
    }
    setMethod("viewLabel", "DynamicGraphView", 
        function(object) 
            .get.env.graphWindow(graphWindow = object)$env$top$env$viewLabel)
    if (!isGeneric("tags")) {
        if (is.function("tags")) 
            fun <- tags
        else fun <- function(object) standardGeneric("tags")
        setGeneric("tags", fun)
    }
    setMethod("tags", "DynamicGraphView", 
        function(object) .get.env.graphWindow(graphWindow = object)$env$tags)
    for (prototype in paste(validViewClasses()[, 2])) 
        setClass(prototype, contains = "DynamicGraphView")
    setClass("dg.Node", 
             representation(color = "character", label = "character", 
                            label.position = "numeric"), 
             prototype = list(color = "black", label = "Label",
                              label.position = c(0, 0, 0)))
    setClass("dg.Vertex", contains = "dg.Node", 
             representation(name = "character", index = "numeric", 
                            position = "numeric", constrained = "logical", 
                            blockindex = "numeric", stratum = "numeric"), 
             prototype = list(color = "black", 
                              label = "Label", 
                              label.position = c(0, 0, 0), 
                              name = "Name", 
                              index = 0, 
                              position = c(0, 0, 0), 
                              constrained = FALSE, 
                              blockindex = 0, stratum = 0))
    if (!isGeneric("setSlots")) {
        if (is.function("setSlots")) 
            fun <- setSlots
        else fun <- function(object, arguments) standardGeneric("setSlots")
        setGeneric("setSlots", fun)
    }
    setMethod("setSlots", "dg.Node", 
        function(object, arguments) {
        for (i in seq(along = arguments)) {
            name <- names(arguments)[i]
            if (is.element(name, slotNames(object))) 
                slot(object, name) <- arguments[[i]]
            else message(paste("Argument '", name, "' not valid slot of '", 
                class(object), "', thus ignored.", sep = ""))
        }
        return(object)
    })
    setMethod("initialize", "dg.Vertex", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("N", names(Args))) {
                N <- Args$N
                .Object@position <- rep(0, N)
                .Object@label.position <- rep(0, N)
                Args <- (Args[!names(Args) == "N"])
            }
            if (class(.Object) == "dg.Node") 
                .Object@color <- "red"
            else if (class(.Object) == "dg.Vertex") 
                .Object@color <- "SaddleBrown"
            else if (class(.Object) == "dg.TextVertex") 
                .Object@color <- "GhostWhite"
            else if (class(.Object) == "dg.DiscreteVertex") 
                .Object@color <- "red"
            else if (class(.Object) == "dg.OrdinalVertex") 
                .Object@color <- "red"
            else if (class(.Object) == "dg.ContinuousVertex") 
                .Object@color <- "red"
            .Object <- setSlots(.Object, Args)
            if (!(length(.Object@position) == length(.Object@label.position))) 
                .Object@label.position <- rep(0, length(.Object@position))
        }
        return(.Object)
    })
    for (prototype in paste(validVertexClasses()[, 2])) 
        setClass(prototype, 
            contains = c("dg.Vertex", "dg.Node"), 
            prototype = list(color = "black", label = "Label",
                             label.position = c(0, 0, 0), name = "Name", 
                             index = 0, position = c(0, 0, 0), 
                             constrained = FALSE, blockindex = 0, stratum = 0))
    setClass("dg.TextVertex", contains = c("dg.Vertex", "dg.Node"), 
        prototype = list(color = "black", label = "Label", 
            label.position = c(0, 0, 0), name = "Name", 
            index = 0, position = c(0, 0, 0), 
            constrained = FALSE, blockindex = 0, stratum = 0))
    setClass("dg.Block", contains = "dg.Node", 
        representation(stratum = "numeric", index = "numeric", 
                       parent = "numeric", children = "numeric", 
                       ancestors = "numeric", descendants = "numeric", 
                       position = "matrix", 
                       closed = "logical", visible = "logical"), 
        prototype = list(color = "black", label = "Label", 
                         label.position = c(0, 0, 0), stratum = 0, 
                         index = 0, parent = 0, children = 0, ancestors = 0, 
                         descendants = 0, 
                         position = matrix(rep(0, 6), ncol = 3), 
                         closed = TRUE, visible = TRUE))
    setMethod("initialize", "dg.Block", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("N", names(Args))) {
                N <- Args$N
                .Object@position <- matrix(rep(0, 2 * N), ncol = N)
                .Object@label.position <- rep(0, N)
                Args <- (Args[!names(Args) == "N"])
            }
            .Object <- setSlots(.Object, Args)
            if ((.Object@parent == 0) && !(.Object@ancestors == 0)) 
                .Object@parent <- .Object@ancestors[length(.Object@ancestors)]
        }
        return(.Object)
    })
    setClass("dg.FactorVertex", contains = c("dg.Vertex", "dg.Node"), 
        representation(vertex.indices = "numeric", 
                       fixed.positions = "logical"), 
        prototype = list(color = "black", label = "Label", 
                         label.position = c(0, 0, 0), name = "Name", 
                         index = 0, position = c(0, 0, 0), 
                         constrained = FALSE, fixed.positions = FALSE, 
                         blockindex = 0, stratum = 0, 
                         vertex.indices = c(0, 0)))
    setMethod("initialize", "dg.FactorVertex", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("vertices", names(Args)) || 
                is.element("vertexList", names(Args))) {
                if (is.element("vertices", names(Args))) {
                  vertices <- Args$vertices
                  if (!(length(vertices) == length(Args$vertex.indices))) 
                    message(paste("Different lengths of 'vertex.indices' ", 
                      "and 'vertices'. Use argument ", 
                      "'vertexList' for subsetting.", sep = ""))
                }
                else vertices <- new("dg.VertexList", 
                                         Args$vertexList[Args$vertex.indices])
                if (length(vertices) > 0) {
                  position <- apply(Positions(vertices), 2, mean)
                  name <- paste(Labels(vertices), collapse = ":")
                  .Object@position <- position
                  if ((length(.Object@name) > 0) && 
                    (.Object@name == "Name")) 
                    .Object@name <- name
                  if ((length(.Object@label) > 0) && 
                    (.Object@label == "Label")) 
                    .Object@label <- name
                }
                Args <- (Args[!names(Args) == "vertices"])
                Args <- (Args[!names(Args) == "vertexList"])
            }
            if (class(.Object) == "dg.Generator") 
                .Object@color <- "yellow"
            else if (class(.Object) == "dg.DiscreteGenerator") 
                .Object@color <- "cyan"
            else if (class(.Object) == "dg.LinearGenerator") 
                .Object@color <- "magenta"
            else if (class(.Object) == "dg.QuadraticGenerator") 
                .Object@color <- "blue"
            .Object <- setSlots(.Object, Args)
        }
        return(.Object)
    })
    for (prototype in paste(validFactorClasses()[, 2])) 
        setClass(prototype, 
            contains = c("dg.FactorVertex", "dg.Vertex", "dg.Node"), 
            prototype = list(color = "black", label = "Label", 
                             label.position = c(0, 0, 0), name = "Name", 
                             index = 0,  position = c(0, 0, 0), 
                             constrained = FALSE, blockindex = 0, 
                             stratum = 0, vertex.indices = c(0, 0)))
    setClass("dg.Edge", contains = "dg.Node", 
             representation(dash = "character", 
                            vertex.indices = "numeric", width = "numeric"), 
             prototype = list(color = "black", label = "Label", 
                              label.position = c(0, 0, 0), dash = "", 
                              vertex.indices = c(0, 0), width = 2))
    setMethod("initialize", "dg.Edge", 
        function(.Object, ...) {
        Args <- list(...)
        if (length(Args) > 0) {
            if (is.element("vertices", names(Args)) || 
                is.element("vertexList", names(Args))) {
                if (is.element("vertices", names(Args))) {
                  vertices <- Args$vertices
                  if (!(length(vertices) == length(Args$vertex.indices))) 
                    message(paste("Different lengths of 'vertex.indices' ", 
                      "and 'vertices'. Use argument ", 
                      "'vertexList' for subsetting.", sep = ""))
                }
                else vertices <- new("dg.VertexList", 
                                        Args$vertexList[Args$vertex.indices])
                if ((length(.Object@label) > 0) && 
                  (.Object@label == "Label") && ((length(vertices) > 0))) 
                  .Object@label <- paste(Labels(vertices), collapse = "~")
                Args <- (Args[!names(Args) == "vertices"])
                Args <- (Args[!names(Args) == "vertexList"])
            }
            if (class(.Object) == "dg.BlockEdge") 
                .Object@color <- "DarkOliveGreen"
            else if (class(.Object) == "dg.FactorEdge") 
                .Object@color <- "DarkOliveGreen"
            else if (class(.Object) == "dg.ExtraEdge") 
                .Object@color <- "DarkOliveGreen"
            else if (class(.Object) == "dg.VertexEdge") 
                .Object@color <- "DarkOliveGreen"
            if (is.element("N", names(Args))) {
                .Object@label.position <- rep(0, Args$N)
                Args <- (Args[!names(Args) == "N"])
            }
            .Object <- setSlots(.Object, Args)
        }
        return(.Object)
    })
    setClass("dg.VertexEdge", contains = c("dg.Edge", "dg.Node"), 
        representation(oriented = "logical"), 
        prototype = list(color = "black", label = "Label", 
                         label.position = c(0, 0, 0), dash = "", 
                         vertex.indices = c(0, 0), width = 2, 
                         oriented = FALSE))
    for (prototype in paste(validEdgeClasses()[-1, 2])) 
        setClass(prototype, 
            contains = c("dg.VertexEdge", "dg.Edge", "dg.Node"))
    setClass("dg.FactorEdge", contains = c("dg.Edge", "dg.Node"), 
        representation(), 
        prototype = list(color = "black", label = "Label", 
                         label.position = c(0, 0, 0), dash = "", 
                         vertex.indices = c(0, 0), width = 2))
    setClass("dg.ExtraEdge", contains = c("dg.Edge", "dg.Node"), 
        representation(), 
        prototype = list(color = "black", label = "Label", 
                         label.position = c(0, 0, 0), dash = "", 
                         vertex.indices = c(0, 0), width = 2))
    setClass("dg.BlockEdge", contains = c("dg.Edge", "dg.Node"), 
        representation(oriented = "logical"), 
        prototype = list(color = "black", label = "Label", 
                         label.position = c(0, 0, 0), dash = "", 
                         vertex.indices = c(0, 0), width = 2, oriented = NA))
    if (!isGeneric("color")) {
        if (is.function("color")) 
            fun <- color
        else fun <- function(object) standardGeneric("color")
        setGeneric("color", fun)
    }
    setMethod("color", "dg.Node", 
        function(object) object@color)
    setGeneric("color<-", 
        function(x, value) standardGeneric("color<-"))
    setReplaceMethod("color", "dg.Node", 
        function(x, value) {
        if (!(value %in% colors())) 
            message(paste("Invalid color: ", value))
        else x@color <- value
        x
    })
    if (!isGeneric("label")) {
        if (is.function("label")) 
            fun <- label
        else fun <- function(object) standardGeneric("label")
        setGeneric("label", fun)
    }
    setMethod("label", "dg.Node", 
        function(object) object@label)
    setGeneric("label<-", 
        function(x, value) standardGeneric("label<-"))
    setReplaceMethod("label", "dg.Node", 
        function(x, value) {
        x@label <- value
        x
    })
    if (!isGeneric("labelPosition")) {
        if (is.function("labelPosition")) 
            fun <- labelPosition
        else fun <- function(object) standardGeneric("labelPosition")
        setGeneric("labelPosition", fun)
    }
    setMethod("labelPosition", "dg.Node", 
        function(object) object@label.position)
    setGeneric("labelPosition<-", 
        function(x, value) standardGeneric("labelPosition<-"))
    setReplaceMethod("labelPosition", "dg.Node", 
        function(x, value) {
        x@label.position <- value
        x
    })
    if (!isGeneric("name")) {
        if (is.function("name")) 
            fun <- name
        else fun <- function(object) standardGeneric("name")
        setGeneric("name", fun)
    }
    setMethod("name", "dg.Vertex", 
        function(object) object@name)
    setMethod("name", "dg.Edge", 
        function(object) object@label)
    setMethod("name", "dg.Block", 
        function(object) object@label)
    setGeneric("name<-", 
        function(x, value) standardGeneric("name<-"))
    setReplaceMethod("name", "dg.Vertex", 
        function(x, value) {
        x@name <- value
        x
    })
    if (!isGeneric("index")) {
        if (is.function("index")) 
            fun <- index
        else fun <- function(object) standardGeneric("index")
        setGeneric("index", fun)
    }
    setMethod("index", "dg.Vertex", 
        function(object) object@index)
    setMethod("index", "dg.Block", 
        function(object) abs(object@index))
    setMethod("index", "dg.FactorVertex", 
        function(object) abs(object@index))
    setGeneric("index<-", 
        function(x, value) standardGeneric("index<-"))
    setReplaceMethod("index", "dg.Vertex", 
        function(x, value) {
        x@index <- value
        x
    })
    setReplaceMethod("index", "dg.Block", 
        function(x, value) {
        x@index <- -abs(value)
        x
    })
    setReplaceMethod("index", "dg.FactorVertex", 
        function(x, value) {
        x@index <- -abs(value)
        x
    })
    setMethod("index", "DynamicGraphModel", 
        function(object) object@index)
    setMethod("index", "DynamicGraphView", 
        function(object) object@index)
    if (!isGeneric("nodeIndices")) {
        if (is.function("nodeIndices")) 
            fun <- nodeIndices
        else fun <- function(object) standardGeneric("nodeIndices")
        setGeneric("nodeIndices", fun)
    }
    setMethod("nodeIndices", "dg.FactorVertex", 
        function(object) object@vertex.indices)
    setMethod("nodeIndices", "dg.Edge", 
        function(object) object@vertex.indices)
    setGeneric("nodeIndices<-", 
        function(x, value) standardGeneric("nodeIndices<-"))
    setReplaceMethod("nodeIndices", "dg.FactorVertex", 
        function(x, value) {
        x@vertex.indices <- value
        x
    })
    if (!isGeneric("constrained")) {
        if (is.function("constrained")) 
            fun <- constrained
        else fun <- function(object) standardGeneric("constrained")
        setGeneric("constrained", fun)
    }
    setMethod("constrained", "dg.Vertex", 
        function(object) if (("constrained" %in% 
        slotNames(object))) 
        object@constrained
    else FALSE)
    setGeneric("constrained<-", 
        function(x, value) standardGeneric("constrained<-"))
    setReplaceMethod("constrained", "dg.Vertex", 
        function(x, value) if (("constrained" %in% slotNames(x))) {
        x@constrained <- value
        x
    })
    if (!isGeneric("fixed.positions")) {
        if (is.function("fixed.positions")) 
            fun <- fixed.positions
        else fun <- function(object) standardGeneric("fixed.positions")
        setGeneric("fixed.positions", fun)
    }
    setMethod("fixed.positions", "dg.FactorVertex", 
        function(object) if (("fixed.positions" %in% 
        slotNames(object))) 
        object@fixed.positions
    else FALSE)
    setGeneric("fixed.positions<-", 
        function(x, value) standardGeneric("fixed.positions<-"))
    setReplaceMethod("fixed.positions", "dg.FactorVertex", 
        function(x, value) 
        if (("fixed.positions" %in% slotNames(x))) {
        x@fixed.positions <- value
        x
    })
    if (!isGeneric("position")) {
        if (is.function("position")) 
            fun <- position
        else fun <- function(object) standardGeneric("position")
        setGeneric("position", fun)
    }
    setMethod("position", "dg.Vertex", 
        function(object) object@position)
    setGeneric("position<-", 
        function(x, value) standardGeneric("position<-"))
    setReplaceMethod("position", "dg.Vertex", 
        function(x, value) {
        x@position <- value
        x
    })
    setMethod("position", "dg.Block", 
        function(object) t(object@position))
    setReplaceMethod("position", "dg.Block", 
        function(x, value) {
        x@position <- t(value)
        x
    })
    if (!isGeneric("stratum")) {
        if (is.function("stratum")) 
            fun <- stratum
        else fun <- function(object) standardGeneric("stratum")
        setGeneric("stratum", fun)
    }
    setMethod("stratum", "dg.Vertex", 
        function(object) object@stratum)
    setGeneric("stratum<-", 
        function(x, value) standardGeneric("stratum<-"))
    setReplaceMethod("stratum", "dg.Vertex", 
        function(x, value) {
        x@stratum <- value
        x
    })
    setMethod("stratum", "dg.Block", 
        function(object) object@stratum)
    setReplaceMethod("stratum", "dg.Block", 
        function(x, value) {
        x@stratum <- value
        x
    })
    if (!isGeneric("blockindex")) {
        if (is.function("blockindex")) 
            fun <- blockindex
        else fun <- function(object) standardGeneric("blockindex")
        setGeneric("blockindex", fun)
    }
    setMethod("blockindex", "dg.Vertex", 
        function(object) object@blockindex)
    setGeneric("blockindex<-", 
        function(x, value) standardGeneric("blockindex<-"))
    setReplaceMethod("blockindex", "dg.Vertex", 
        function(x, value) {
        x@blockindex <- value
        x
    })
    if (!isGeneric("parent")) {
        if (is.function("parent")) 
            fun <- parent
        else fun <- function(object) standardGeneric("parent")
        setGeneric("parent", fun)
    }
    setMethod("parent", "dg.Block", 
        function(object) abs(object@parent))
    setGeneric("parent<-", 
        function(x, value) standardGeneric("parent<-"))
    setReplaceMethod("parent", "dg.Block", 
        function(x, value) {
        x@parent <- value
        x
    })
    if (!isGeneric("children")) {
        if (is.function("children")) 
            fun <- children
        else fun <- function(object, ...) standardGeneric("children")
        setGeneric("children", fun)
    }
    setMethod("children", "dg.Block", 
        function(object, ...) object@children)
    setGeneric("children<-", 
        function(x, value) standardGeneric("children<-"))
    setReplaceMethod("children", "dg.Block", 
        function(x, value) {
        x@children <- value
        x
    })
    if (!isGeneric("ancestors")) {
        if (is.function("ancestors")) 
            fun <- ancestors
        else fun <- function(object, blockList = NULL, vertexList = NULL, 
            ...) standardGeneric("ancestors")
        setGeneric("ancestors", fun)
    }
    setMethod("ancestors", "dg.Vertex", 
        function(object, blockList = NULL, 
        vertexList = NULL, ...) warning("Not implemented"))
    setGeneric("ancestors<-", 
        function(x, value) standardGeneric("ancestors<-"))
    setReplaceMethod("ancestors", "dg.Vertex", 
        function(x, value) {
        warning("Not implemented")
    })
    setMethod("ancestors", "dg.Block", 
        function(object, blockList = NULL, 
        vertexList = NULL, ...) object@ancestors)
    setReplaceMethod("ancestors", "dg.Block", 
        function(x, value) {
        x@ancestors <- value
        x
    })
    if (!isGeneric("descendants")) {
        if (is.function("descendants")) 
            fun <- descendants
        else fun <- function(object, blockList = NULL, vertexList = NULL, ...) 
            standardGeneric("descendants")
        setGeneric("descendants", fun)
    }
    setMethod("descendants", "dg.Vertex", 
        function(object, blockList = NULL, 
        vertexList = NULL, ...) warning("Not implemented"))
    setGeneric("descendants<-", 
        function(x, value) standardGeneric("descendants<-"))
    setReplaceMethod("descendants", "dg.Vertex", 
        function(x, 
        value) {
        warning("Not implemented")
    })
    setMethod("descendants", "dg.Block", 
        function(object, blockList = NULL, 
        vertexList = NULL, ...) object@descendants)
    setReplaceMethod("descendants", "dg.Block", 
        function(x, value) {
        x@descendants <- value
        x
    })
    if (!isGeneric("closed")) {
        if (is.function("closed")) 
            fun <- closed
        else fun <- function(object) standardGeneric("closed")
        setGeneric("closed", fun)
    }
    setMethod("closed", "dg.Block", 
        function(object) object@closed)
    setGeneric("closed<-", 
        function(x, value) standardGeneric("closed<-"))
    setReplaceMethod("closed", "dg.Block", 
        function(x, value) {
        x@closed <- value
        x
    })
    if (!isGeneric("visible")) {
        if (is.function("visible")) 
            fun <- visible
        else fun <- function(object) standardGeneric("visible")
        setGeneric("visible", fun)
    }
    setMethod("visible", "dg.Vertex", 
        function(object) warning("Not implemented"))
    setGeneric("visible<-", 
        function(x, value) standardGeneric("visible<-"))
    setReplaceMethod("visible", "dg.Vertex", 
        function(x, value) {
        warning("Not implemented")
    })
    setMethod("visible", "dg.Block", 
        function(object) object@visible)
    setReplaceMethod("visible", "dg.Block", 
        function(x, value) {
        x@visible <- value
        x
    })
    if (!isGeneric("draw")) {
        if (is.function("draw")) 
            fun <- draw
        else fun <- function(object, canvas, position, x = position[1], 
            y = position[2], stratum = 0, w = 2, color = "green", 
            background = "white", ...) standardGeneric("draw")
        setGeneric("draw", fun)
    }
    setMethod("draw", "dg.TextVertex", 
        function(object, canvas, 
        position, x = position[1], y = position[2], stratum = 0, 
        w = 2, color = "black", background = "white") {
        s <- w * sqrt(4/pi)/10
        p <- tkcreate(canvas, "oval", x - s, y - s, x + s, y + s, 
            fill = color(object), activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.DiscreteVertex", 
        function(object, canvas, 
        position, x = position[1], y = position[2], stratum = 0, 
        w = 2, color = "green", background = "white") {
        s <- w * sqrt(4/pi)
        p <- tkcreate(canvas, "oval", x - s, y - s, x + s, y + s, 
            fill = color(object), activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.OrdinalVertex", 
        function(object, canvas, 
        position, x = position[1], y = position[2], stratum = 0, 
        w = 2, color = "green", background = "white") {
        p <- tkcreate(canvas, "rectangle", x - w, y - w, x + w, y + w, 
            fill = color(object), activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.ContinuousVertex", 
        function(object, 
        canvas, position, x = position[1], y = position[2], stratum = 0, 
        w = 2, color = "green", background = "white") {
        s <- w * sqrt(4/pi)
        p <- tkcreate(canvas, "oval", x - s, y - s, x + s, y + s, 
            fill = color(object), activefill = "IndianRed")
        s <- 0.6 * w * sqrt(4/pi)
        q <- tkcreate(canvas, "oval", x - s, y - s, x + s, y + s, 
            fill = background)
        return(list(dynamic = list(p), fixed = list(q)))
    })
    if (!isGeneric("propertyDialog")) {
        if (is.function("propertyDialog")) 
            fun <- propertyDialog
        else fun <- function(object, classes = NULL, title = class(object), 
                             sub.title = label(object), 
                             name.object = name(object), okReturn = TRUE, 
                             fixedSlots = NULL, difficultSlots = NULL, 
                             top = NULL, entryWidth = 20, do.grab = FALSE) 
            standardGeneric("propertyDialog")
        setGeneric("propertyDialog", fun)
    }
    setMethod("propertyDialog", "dg.Node", 
        function(object, classes = NULL, title = class(object), 
                 sub.title = label(object), name.object = name(object), 
                 okReturn = TRUE, fixedSlots = NULL, difficultSlots = NULL, 
                 top = NULL, entryWidth = 20, do.grab = FALSE) {
        .propertyDialog(object, classes = classes, title = title, 
                        sub.title = sub.title, name.object = name.object, 
                        okReturn = okReturn, fixedSlots = fixedSlots, 
                        difficultSlots = difficultSlots, 
                        top = top, entryWidth = entryWidth, do.grab = do.grab)
    })
    if (!isGeneric("addToPopups")) {
        if (is.function("addToPopups")) 
            fun <- draw
        else fun <- function(object, type, nodePopupMenu, i, 
            updateArguments, Args, ...) standardGeneric("addToPopups")
        setGeneric("addToPopups", fun)
    }
    setMethod("addToPopups", "dg.Node", 
        function(object, type, 
        nodePopupMenu, i, updateArguments, Args, ...) {
        if ((type == "Factor")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a factor!"), 
                command = function() {
                  str(i)
                })
        if ((type == "Vertex")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a vertex!"), 
                command = function() {
                  print(i)
                  updateArguments(NULL)
                  Arguments = Args()
                  print(Arguments$vertexList[[i]]@label)
                })
        else if ((type == "OpenBlock")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is an open block!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "ClosedBlock")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a closed block!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "Edge")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is an edge!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "VertexEdge")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a graph edge!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "FactorEdge")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a factor edge!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "ExtraEdge")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a extra edge!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "BlockEdge")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is a block edge!"), 
                command = function() {
                  str(i)
                })
        else if ((type == "Extra")) 
            tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is an extra!"), 
                command = function() {
                  str(i)
                })
        else tkadd(nodePopupMenu, "command", 
                label = paste(" --- This is ??? "), 
            command = function() {
                str(i)
            })
    })
    setMethod("draw", "dg.Generator", 
        function(object, canvas, position, x = position[1], 
                 y = position[2], stratum = 0, w = 2, color = "green", 
                 background = "white") {
        s <- w * sqrt(2)
        p <- tkcreate(canvas, "polygon", x - 0, y + s, x + s, 
            y + 0, x + 0, y - s, x - s, y - 0, fill = color(object), 
            activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.DiscreteGenerator", 
        function(object, canvas, position, x = position[1], 
                 y = position[2], stratum = 0, w = 2, color = "green", 
                 background = "white") {
        s <- w * sqrt(2)
        p <- tkcreate(canvas, "polygon", x - w, y + w, x + w, 
            y + w, x + w, y - w, x - w, y - w, fill = color(object), 
            activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.LinearGenerator", 
        function(object, canvas, position, x = position[1], 
                 y = position[2], stratum = 0, w = 2, color = "green", 
                 background = "white") {
        s <- w * sqrt(2)
        p <- tkcreate(canvas, "polygon", x - w, y + w, x + w, 
            y + w, x + w, y - w, x - w, y - w, fill = color(object), 
            activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.QuadraticGenerator", 
        function(object, canvas, position, x = position[1], 
                 y = position[2], stratum = 0, w = 2, color = "green", 
                 background = "white") {
        s <- w * sqrt(2)
        p <- tkcreate(canvas, "polygon", x - w, y + w, x + w, 
            y + w, x + w, y - w, x - w, y - w, fill = color(object), 
            activefill = "IndianRed")
        return(list(dynamic = list(p), fixed = NULL))
    })
    setMethod("draw", "dg.Edge", 
        function(object, canvas, position, 
                 x = lapply(position, function(e) e[1]), 
                 y = lapply(position, function(e) e[2]), 
                 stratum = as.vector(rep(0, length(position)), mode = "list"), 
                 w = 2, color = "green", background = "white", 
                 font.edge.label = "8x16") {
        "f" <- function(i, j) {
            arrowhead <- "none"
            if ((class(object) == "dg.VertexEdge") || (class(object) == 
                "dg.BlockEdge")) 
                if (!is.na(object@oriented) && object@oriented) 
                  arrowhead <- "last"
                else if (stratum[[i]] == stratum[[j]]) 
                  arrowhead <- "none"
                else if (stratum[[i]] < stratum[[j]]) 
                  arrowhead <- "last"
                else arrowhead <- "first"
            if (class(object) == "dg.DoubleArrowEdge") 
                arrowhead <- "both"
            dash <- dash(object)
            if (class(object) == "dg.DashedEdge") 
                dash <- "-"
            l <- function(xi, yi, xj, yj) 
                tkcreate(canvas, "line", xi, yi, xj, yj, width = w, 
                         arrow = arrowhead, dash = dash, fill = color(object), 
                         activefill = "DarkSlateGray")
            if ((class(object) == "dg.DoubleConnectedEdge") || 
                (class(object) == "dg.TripleConnectedEdge")) {
                dx <- x[[i]] - x[[j]]
                dy <- y[[i]] - y[[j]]
                ld <- sqrt(dx^2 + dy^2)
                dx <- w * dx/ld
                dy <- w * dy/ld
                line1 <- l(x[[i]] - dy - dx, y[[i]] + dx - dy, 
                  x[[j]] - dy + dx, y[[j]] + dx + dy)
                line2 <- l(x[[i]] + dy - dx, y[[i]] - dx - dy, 
                  x[[j]] + dy + dx, y[[j]] - dx + dy)
                if ((class(object) == "dg.TripleConnectedEdge")) {
                  line0 <- l(x[[i]], y[[i]], x[[j]], y[[j]])
                  lines <- list(line1, line0, line2)
                }
                else lines <- list(line1, line2)
            }
            else lines <- list(l(x[[i]], y[[i]], x[[j]], y[[j]]))
            label.position <- (position[[i]] + position[[j]])/2
            pos <- label.position + rep(0, length(label.position))
            label <- tkcreate(canvas, "text", pos[1], pos[2], 
                              text = object@label, anchor = "nw", 
                              font = font.edge.label, 
                              activefill = "DarkSlateGray")
            tags <- NULL
            if (class(object) == "dg.DottedEdge") {
                x. <- mean(unlist(x))
                y. <- mean(unlist(y))
                s <- w * sqrt(4/pi)
                p <- tkcreate(canvas, "oval", x. - s, y. - s, 
                              x. + s, y. + s, fill = color(object), 
                              activefill = "SeaGreen")
                tags <- list(p)
            }
            return(list(lines = lines, tags = tags, 
                        from = object@vertex.indices[i], 
                        to = object@vertex.indices[j], label = label, 
                        label.position = label.position))
        }
        result <- NULL
        edge <- object@vertex.indices
        m <- length(edge)
        for (j in seq(along = edge)) if (j < length(edge)) 
            for (k in (j + 1):length(edge)) result <- append(result, 
                list(f(j, k)))
        return(result)
    })
    if (!isGeneric("width")) {
        if (is.function("width")) 
            fun <- width
        else fun <- function(object) standardGeneric("width")
        setGeneric("width", fun)
    }
    setMethod("width", "dg.Edge", 
        function(object) object@width)
    setGeneric("width<-", 
        function(x, value) standardGeneric("width<-"))
    setReplaceMethod("width", "dg.Edge", 
        function(x, value) {
        x@width <- value
        x
    })
    if (!isGeneric("dash")) {
        if (is.function("dash")) 
            fun <- dash
        else fun <- function(object) standardGeneric("dash")
        setGeneric("dash", fun)
    }
    setMethod("dash", "dg.Edge", 
        function(object) object@dash)
    setGeneric("dash<-", 
        function(x, value) standardGeneric("dash<-"))
    setReplaceMethod("dash", "dg.Edge", 
        function(x, value) {
        .dashReplaceMethod(x, value)
    })
    if (!isGeneric("oriented")) {
        if (is.function("oriented")) 
            fun <- oriented
        else fun <- function(object) standardGeneric("oriented")
        setGeneric("oriented", fun)
    }
    setMethod("oriented", "dg.VertexEdge", 
        function(object) object@oriented)
    setGeneric("oriented<-", 
        function(x, value) standardGeneric("oriented<-"))
    setReplaceMethod("oriented", "dg.VertexEdge", 
        function(x, value) {
        x@oriented <- value
        x
    })
    setMethod("oriented", "dg.BlockEdge", 
        function(object) object@oriented)
    setReplaceMethod("oriented", "dg.BlockEdge", 
        function(x, value) {
        x@oriented <- value
        x
    })
    setMethod("oriented", "dg.graphedges", 
        function(object) object@oriented)
    if (!isGeneric("nodeTypesOfEdge")) {
        if (is.function("nodeTypesOfEdge")) 
            fun <- nodeTypesOfEdge
        else fun <- function(object) standardGeneric("nodeTypesOfEdge")
        setGeneric("nodeTypesOfEdge", fun)
    }
    setMethod("nodeTypesOfEdge", "dg.VertexEdge", 
        function(object) rep("Vertex", 
        length(object@vertex.indices)))
    setMethod("nodeTypesOfEdge", "dg.FactorEdge", 
        function(object) ifelse(object@vertex.indices > 0, "Vertex", "Factor"))
    setMethod("nodeTypesOfEdge", "dg.ExtraEdge", 
        function(object) ifelse(object@vertex.indices > 0, "Vertex", "Extra"))
    setMethod("nodeTypesOfEdge", "dg.BlockEdge", 
        function(object) ifelse(object@vertex.indices > 0, "Vertex", "Block"))
    if (!isGeneric("nodeIndicesOfEdge")) {
        if (is.function("nodeIndicesOfEdge")) 
            fun <- nodeIndicesOfEdge
        else fun <- function(object) standardGeneric("nodeIndicesOfEdge")
        setGeneric("nodeIndicesOfEdge", fun)
    }
    setMethod("nodeIndicesOfEdge", "dg.Edge", 
        function(object) object@vertex.indices)
    setGeneric("nodeIndicesOfEdge<-", 
        function(x, value) standardGeneric("nodeIndicesOfEdge<-"))
    setReplaceMethod("nodeIndicesOfEdge", "dg.Edge", 
        function(x, value) {
        x@vertex.indices <- value
        x
    })
    setMethod("draw", "dg.Block", 
        function(object, canvas, position, 
                 x = position[1], y = position[2], stratum = 0, w = 10, 
                 color = "green", background = "white") {
        s <- w
        p <- tkcreate(canvas, "rectangle", x - s, y - s, x + s, y + s, 
                      fill = color(object), activefill = "IndianRed")
        s <- w - 2
        q <- tkcreate(canvas, "rectangle", x - s, y - s, x + s, y + s, 
                      fill = background)
        return(list(dynamic = list(p), fixed = list(q)))
    })
    if (!isGeneric("Names")) {
        if (is.function("Names")) 
            fun <- Names
        else fun <- function(objects) standardGeneric("Names")
        setGeneric("Names", fun)
    }
    setMethod("Names", "dg.list", 
        function(objects) {
        NAMES <- lapply(objects, 
                        function(x) if (!is.null(x)) name(x))
        names(NAMES) <- NULL
        return(unlist(NAMES))
    })
    setGeneric("Names<-", 
        function(objectlist, value) standardGeneric("Names<-"))
    setReplaceMethod("Names", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) objectlist[[i]]@name <- value[i]
        }
        else warning("Invalid list of values for names")
        objectlist
    })
    if (!isGeneric("Colors")) {
        if (is.function("Colors")) 
            fun <- Colors
        else fun <- function(objectlist) standardGeneric("Colors")
        setGeneric("Colors", fun)
    }
    setMethod("Colors", "dg.list", 
        function(objectlist) {
        Colors <- lapply(objectlist, 
                         function(x) if (!is.null(x)) color(x))
        names(Colors) <- Names(objectlist)
        return(unlist(Colors))
    })
    setGeneric("Colors<-", 
        function(objectlist, value) standardGeneric("Colors<-"))
    setReplaceMethod("Colors", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@color <- value[i]
        }
        else warning("Invalid list of values for colors")
        objectlist
    })
    if (!isGeneric("Dashes")) {
        if (is.function("Dashes")) 
            fun <- Dashes
        else fun <- function(objectlist) standardGeneric("Dashes")
        setGeneric("Dashes", fun)
    }
    setMethod("Dashes", "dg.list", 
        function(objectlist) {
        Dashes <- lapply(objectlist, 
                         function(x) if (!is.null(x)) dash(x))
        names(Dashes) <- Names(objectlist)
        return(unlist(Dashes))
    })
    setGeneric("Dashes<-", 
        function(objectlist, value) standardGeneric("Dashes<-"))
    setReplaceMethod("Dashes", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@dash <- value[i]
        }
        else warning("Invalid list of values for dashs")
        objectlist
    })
    if (!isGeneric("Labels")) {
        if (is.function("Labels")) 
            fun <- Labels
        else fun <- function(objectlist) standardGeneric("Labels")
        setGeneric("Labels", fun)
    }
    setMethod("Labels", "dg.list", 
        function(objectlist) {
        Labels <- lapply(objectlist, 
                         function(x) if (!is.null(x)) label(x))
        names(Labels) <- Names(objectlist)
        return(unlist(Labels))
    })
    setGeneric("Labels<-", 
        function(objectlist, value) standardGeneric("Labels<-"))
    setReplaceMethod("Labels", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@label <- value[i]
        }
        else warning("Invalid list of values for labels")
        objectlist
    })
    if (!isGeneric("LabelPositions")) {
        if (is.function("LabelPositions")) 
            fun <- LabelPositions
        else fun <- function(objectlist) standardGeneric("LabelPositions")
        setGeneric("LabelPositions", fun)
    }
    setMethod("LabelPositions", "dg.list", 
        function(objectlist) {
        positions <- lapply(objectlist, 
                            function(x) if (!is.null(x)) labelPosition(x))
        N.list <- unlist(lapply(positions, function(x) length(x)))
        if ((N.list[1] > 0) && all(N.list == N.list[1])) {
            positions <- matrix(unlist(positions), 
                                ncol = length(positions[[1]]), byrow = TRUE)
            labels <- c("X", "Y")
            if (ncol(positions) > 2) 
                labels <- c(labels, paste("Z", 3:ncol(positions), 
                  sep = "-"))
            dimnames(positions) <- list(Names(objectlist), labels)
        }
        return(positions)
    })
    setGeneric("LabelPositions<-", 
        function(objectlist, value) standardGeneric("LabelPositions<-"))
    setReplaceMethod("LabelPositions", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == nrow(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@label.position <- value[i, ]
        }
        else warning("Invalid list of values for labelpositions")
        objectlist
    })
    if (!isGeneric("FixedPositions")) {
        if (is.function("FixedPositions")) 
            fun <- FixedPositions
        else fun <- function(objectlist) standardGeneric("FixedPositions")
        setGeneric("FixedPositions", fun)
    }
    setMethod("FixedPositions", "dg.list", 
        function(objectlist) {
        fixed.positions <- lapply(objectlist, 
                             function(x) if (!is.null(x)) fixed.positions(x))
        names(fixed.positions) <- Names(objectlist)
        return(unlist(fixed.positions))
    })
    setGeneric("FixedPositions<-", 
        function(objectlist, value) standardGeneric("FixedPositions<-"))
    setReplaceMethod("FixedPositions", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) if (("fixed.positions" %in% 
                slotNames(objectlist[[i]]))) 
                objectlist[[i]]@fixed.positions <- value[i]
        }
        else warning("Invalid list of values for fixed.positions")
        objectlist
    })
    if (!isGeneric("Constrained")) {
        if (is.function("Constrained")) 
            fun <- Constrained
        else fun <- function(objectlist) standardGeneric("Constrained")
        setGeneric("Constrained", fun)
    }
    setMethod("Constrained", "dg.list", 
        function(objectlist) {
        constrained <- lapply(objectlist, 
                              function(x) if (!is.null(x)) constrained(x))
        names(constrained) <- Names(objectlist)
        return(unlist(constrained))
    })
    setGeneric("Constrained<-", 
        function(objectlist, value) standardGeneric("Constrained<-"))
    setReplaceMethod("Constrained", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) if (("constrained" %in% 
                slotNames(objectlist[[i]]))) 
                objectlist[[i]]@constrained <- value[i]
        }
        else warning("Invalid list of values for constrained")
        objectlist
    })
    if (!isGeneric("Positions")) {
        if (is.function("Positions")) 
            fun <- Positions
        else fun <- function(objectlist) standardGeneric("Positions")
        setGeneric("Positions", fun)
    }
    setMethod("Positions", "dg.list", 
        function(objectlist) {
        if (length(objectlist) == 0) 
            return(NULL)
        positions <- lapply(objectlist, 
                            function(x) if (!is.null(x)) position(x))
        if (class(objectlist[[1]]) == "dg.Block") {
            N <- nrow(positions[[1]])
            N.list <- unlist(lapply(positions, function(x) nrow(x)))
            if ((N > 0) && all(N.list == N)) {
                positions <- matrix(unlist(positions), 
                                    ncol = 2 * N, byrow = TRUE)
                labels <- c("X", "Y", "x", "y")
                if (N > 2) 
                  labels <- c("X", "Y", paste("Z", 3:N, sep = "-"), 
                    "x", "y", paste("z", 3:N, sep = "-"))
                if (ncol(positions) > 2) 
                  dimnames(positions) <- list(Names(objectlist), labels)
            }
        }
        else {
            N.list <- unlist(lapply(positions, function(x) length(x)))
            if ((N.list[1] > 0) && all(N.list == N.list[1])) {
                positions <- matrix(unlist(positions), 
                                    ncol = length(positions[[1]]), 
                                    byrow = TRUE)
                labels <- c("X", "Y")
                if (ncol(positions) > 2) 
                  labels <- c(labels, paste("Z", 3:ncol(positions), 
                    sep = "-"))
                dimnames(positions) <- list(Names(objectlist), 
                  labels)
            }
        }
        return(positions)
    })
    setGeneric("Positions<-", 
        function(objectlist, value) standardGeneric("Positions<-"))
    setReplaceMethod("Positions", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == nrow(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@position <- value[i, ]
        }
        else warning("Invalid list of values for positions")
        objectlist
    })
    if (!isGeneric("Closed")) {
        if (is.function("Closed")) 
            fun <- Closed
        else fun <- function(objectlist) standardGeneric("Closed")
        setGeneric("Closed", fun)
    }
    setMethod("Closed", "dg.list", 
        function(objectlist) {
        closed <- lapply(objectlist, 
                         function(x) if (!is.null(x)) closed(x))
        names(closed) <- Names(objectlist)
        return(unlist(closed))
    })
    setGeneric("Closed<-", 
        function(objectlist, value) standardGeneric("Closed<-"))
    setReplaceMethod("Closed", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@closed <- value[i]
        }
        else warning("Invalid list of values for closed")
        objectlist
    })
    if (!isGeneric("Strata")) {
        if (is.function("Strata")) 
            fun <- Strata
        else fun <- function(objectlist) standardGeneric("Strata")
        setGeneric("Strata", fun)
    }
    setMethod("Strata", "dg.list", 
        function(objectlist) {
        strata <- lapply(objectlist, 
                         function(x) if (!is.null(x)) stratum(x))
        names(strata) <- Names(objectlist)
        return(unlist(strata))
    })
    setGeneric("Strata<-", 
        function(objectlist, value) standardGeneric("Strata<-"))
    setReplaceMethod("Strata", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@stratum <- value[i]
        }
        else warning("Invalid list of values for strata")
        objectlist
    })
    if (!isGeneric("Parents")) {
        if (is.function("Parents")) 
            fun <- Parents
        else fun <- function(objectlist) standardGeneric("Parents")
        setGeneric("Parents", fun)
    }
    setMethod("Parents", "dg.list", 
        function(objectlist) {
        parents <- lapply(objectlist, 
                          function(x) if (!is.null(x)) parent(x))
        names(parents) <- Names(objectlist)
        return(parents)
    })
    setGeneric("Parents<-", 
        function(objectlist, value) standardGeneric("Parents<-"))
    setReplaceMethod("Parents", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@parent <- value[i]
        }
        else warning("Invalid list of values for parents")
        objectlist
    })
    if (!isGeneric("Children")) {
        if (is.function("Children")) 
            fun <- Children
        else fun <- function(objectlist) standardGeneric("Children")
        setGeneric("Children", fun)
    }
    setMethod("Children", "dg.list", 
        function(objectlist) {
        parents <- lapply(objectlist, 
                          function(x) if (!is.null(x)) children(x))
        names(parents) <- Names(objectlist)
        return(parents)
    })
    setGeneric("Children<-", 
        function(objectlist, value) standardGeneric("Children<-"))
    setReplaceMethod("Children", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@children <- value[i]
        }
        else warning("Invalid list of values for children")
        objectlist
    })
    if (!isGeneric("NodeAncestors")) {
        if (is.function("NodeAncestors")) 
            fun <- NodeAncestors
        else fun <- function(objectlist) standardGeneric("NodeAncestors")
        setGeneric("NodeAncestors", fun)
    }
    setMethod("NodeAncestors", "dg.list", 
        function(objectlist) {
        ancestors <- lapply(objectlist, 
                            function(x) if (!is.null(x)) ancestors(x))
        names(ancestors) <- Names(objectlist)
        return(ancestors)
    })
    setGeneric("NodeAncestors<-", 
        function(objectlist, value) standardGeneric("NodeAncestors<-"))
    setReplaceMethod("NodeAncestors", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@ancestors <- value[i]
        }
        else warning("Invalid list of values for ancestors")
        objectlist
    })
    if (!isGeneric("NodeDescendants")) {
        if (is.function("NodeDescendants")) 
            fun <- NodeDescendants
        else fun <- function(objectlist) standardGeneric("NodeDescendants")
        setGeneric("NodeDescendants", fun)
    }
    setMethod("NodeDescendants", "dg.list", 
        function(objectlist) {
        descendants <- lapply(objectlist, 
                              function(x) if (!is.null(x)) descendants(x))
        names(descendants) <- Names(objectlist)
        return(descendants)
    })
    setGeneric("NodeDescendants<-", 
        function(objectlist, value) standardGeneric("NodeDescendants<-"))
    setReplaceMethod("NodeDescendants", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@descendants <- value[i]
        }
        else warning("Invalid list of values for descendants")
        objectlist
    })
    if (!isGeneric("Indices")) {
        if (is.function("Indices")) 
            fun <- Indices
        else fun <- function(objectlist) standardGeneric("Indices")
        setGeneric("Indices", fun)
    }
    setMethod("Indices", "dg.list", 
        function(objectlist) {
        indices <- lapply(objectlist, 
                          function(x) if (!is.null(x)) index(x))
        names(indices) <- Names(objectlist)
        return(unlist(indices))
    })
    setGeneric("Indices<-", 
        function(objectlist, value) standardGeneric("Indices<-"))
    setReplaceMethod("Indices", "dg.list", 
        function(objectlist, 
        value) {
        for (i in seq(along = objectlist)) 
            objectlist[[i]]@index <- value[i]
        objectlist
    })
    if (!isGeneric("Blockindices")) {
        if (is.function("Blockindices")) 
            fun <- Blockindices
        else fun <- function(objectlist) standardGeneric("Blockindices")
        setGeneric("Blockindices", fun)
    }
    setMethod("Blockindices", "dg.list", 
        function(objectlist) {
        blockindices <- lapply(objectlist, 
                               function(x) if (!is.null(x)) blockindex(x))
        names(blockindices) <- Names(objectlist)
        return(unlist(blockindices))
    })
    if (!isGeneric("NodeTypes")) {
        if (is.function("NodeTypes")) 
            fun <- NodeTypes
        else fun <- function(objectlist) standardGeneric("NodeTypes")
        setGeneric("NodeTypes", fun)
    }
    setMethod("NodeTypes", "dg.list", 
        function(objectlist) {
        indices <- lapply(objectlist, 
                          function(x) if (!is.null(x)) nodeTypesOfEdge(x))
        names(indices) <- Names(objectlist)
        return(indices)
    })
    if (!isGeneric("NodeIndices")) {
        if (is.function("NodeIndices")) 
            fun <- NodeIndices
        else fun <- function(objectlist) standardGeneric("NodeIndices")
        setGeneric("NodeIndices", fun)
    }
    setMethod("NodeIndices", "dg.list", 
        function(objectlist) {
        indices <- lapply(objectlist, 
                          function(x) if (!is.null(x)) nodeIndicesOfEdge(x))
        names(indices) <- Names(objectlist)
        return(indices)
    })
    if (!isGeneric("Widths")) {
        if (is.function("Widths")) 
            fun <- Widths
        else fun <- function(objectlist) standardGeneric("Widths")
        setGeneric("Widths", fun)
    }
    setMethod("Widths", "dg.list", 
        function(objectlist) {
        widths <- lapply(objectlist, 
                         function(x) if (!is.null(x)) width(x))
        names(widths) <- Names(objectlist)
        return(unlist(widths))
    })
    setGeneric("Widths<-", 
        function(objectlist, value) standardGeneric("Widths<-"))
    setReplaceMethod("Widths", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@width <- value[i]
        }
        else warning("Invalid list of values for widths")
        objectlist
    })
    if (!isGeneric("Oriented")) {
        if (is.function("Oriented")) 
            fun <- Oriented
        else fun <- function(objectlist) standardGeneric("Oriented")
        setGeneric("Oriented", fun)
    }
    setMethod("Oriented", "dg.list", 
        function(objectlist) {
        oriented <- lapply(objectlist, 
                           function(x) if (!is.null(x)) oriented(x))
        names(oriented) <- Names(objectlist)
        return(unlist(oriented))
    })
    setGeneric("Oriented<-", 
        function(objectlist, value) standardGeneric("Oriented<-"))
    setReplaceMethod("Oriented", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@oriented <- value[i]
        }
        else warning("Invalid list of values for oriented")
        objectlist
    })
    if (!isGeneric("Visible")) {
        if (is.function("Visible")) 
            fun <- Visible
        else fun <- function(objectlist) standardGeneric("Visible")
        setGeneric("Visible", fun)
    }
    setMethod("Visible", "dg.list", 
        function(objectlist) {
        visible <- lapply(objectlist, 
                          function(x) if (!is.null(x)) visible(x))
        names(visible) <- Names(objectlist)
        return(unlist(visible))
    })
    setGeneric("Visible<-", 
        function(objectlist, value) standardGeneric("Visible<-"))
    setReplaceMethod("Visible", "dg.list", 
        function(objectlist, 
        value) {
        if (length(objectlist) == length(value)) {
            for (i in seq(along = objectlist)) 
                objectlist[[i]]@visible <- value[i]
        }
        else warning("Invalid list of values for visible")
        objectlist
    })
    setClass("dg.Model", 
        representation(name = "character", dg = "dg.graphedges"))
    setMethod("setSlots", "dg.Model", 
        function(object, arguments) {
        for (i in seq(along = arguments)) {
            name <- names(arguments)[i]
            if (is.element(name, slotNames(object))) 
                slot(object, name) <- arguments[[i]]
            else message(paste("Argument '", name, "' not valid slot of '", 
                class(object), "', thus ignored.", sep = ""))
        }
        return(object)
    })
    setMethod("initialize", "dg.Model", 
        function(.Object, ...) {
        Args <- list(...)
        .Object <- setSlots(.Object, Args)
        return(.Object)
    })
    if (!isGeneric("graphComponents")) {
        if (is.function("graphComponents")) 
            fun <- graphComponents
        else fun <- function(object, viewType = NULL, ...) 
                standardGeneric("graphComponents")
        setGeneric("graphComponents", fun)
    }
    if (!isGeneric("graphEdges")) {
        if (is.function("graphEdges")) 
            fun <- graphEdges
        else fun <- function(object, viewType = NULL, ...) 
                standardGeneric("graphEdges")
        setGeneric("graphEdges", fun)
    }
    setMethod("graphEdges", "dg.Model", 
        function(object, viewType = NULL, ...) {
        dots <- list(...)
        localArguments <- dots$Arguments
        Vertices <- localArguments$vertexList
        Edges <- object@dg@edgeList
        VisibleVertices <- object@dg@visibleVertices
        if (viewType == "Factor") {
            factors <- .cliquesFromEdges(Edges, Vertices, VisibleVertices)
            if (is.null(factors) || (length(factors) == 0)) {
                FactorVertices <- new("dg.FactorVertexList")
                FactorEdges <- new("dg.FactorEdgeList")
            }
            else {
                result <- returnFactorVerticesAndEdges(Vertices, 
                  factors)
                FactorVertices <- result$FactorVertices
                FactorEdges <- result$FactorEdges
            }
            new("dg.graphedges", viewType = viewType,  
                oriented = object@dg@oriented, 
                edgeList = object@dg@edgeList, 
                blockEdgeList = object@dg@blockEdgeList, 
                factorVertexList = FactorVertices, 
                factorEdgeList = FactorEdges, 
                visibleVertices = object@dg@visibleVertices, 
                visibleBlocks = object@dg@visibleBlocks, 
                extraList = object@dg@extraList, 
                extraEdgeList = object@dg@extraEdgeList)
        }
        else if (viewType == "Moral") {
            message("Moral view not implemented; ")
            new("dg.graphedges", viewType = viewType, 
                oriented = object@dg@oriented, 
                edgeList = object@dg@edgeList, 
                visibleVertices = object@dg@visibleVertices, 
                visibleBlocks = numeric(), 
                extraList = object@dg@extraList, 
                extraEdgeList = object@dg@extraEdgeList)
        }
        else if (viewType == "Essential") {
            message("Essential view not implemented; ")
            new("dg.graphedges", viewType = viewType, 
                oriented = object@dg@oriented, 
                edgeList = object@dg@edgeList, 
                visibleVertices = object@dg@visibleVertices, 
                visibleBlocks = numeric(), 
                extraList = object@dg@extraList, 
                extraEdgeList = object@dg@extraEdgeList)
        }
        else if (viewType == "Simple") {
            new("dg.graphedges", viewType = viewType, 
                oriented = object@dg@oriented, 
                edgeList = object@dg@edgeList, 
                blockEdgeList = object@dg@blockEdgeList, 
                visibleVertices = object@dg@visibleVertices, 
                visibleBlocks = object@dg@visibleBlocks, 
                extraList = object@dg@extraList, 
                extraEdgeList = object@dg@extraEdgeList)
        }
        else message("View type not implemented; ")
    })
    if (!isGeneric("setGraphComponents")) {
        if (is.function("setGraphComponents")) 
            fun <- setGraphComponents
        else fun <- function(object, viewType = NULL, visibleVertices = NULL, 
            visibleBlocks = NULL, extraVertices = NULL, vertexEdges = NULL, 
            blockEdges = NULL, factorVertices = NULL, factorEdges = NULL, 
            extraEdges = NULL, ...) standardGeneric("setGraphComponents")
        setGeneric("setGraphComponents", fun)
    }
    if (!isGeneric("setGraphEdges")) {
        if (is.function("setGraphEdges")) 
            fun <- setGraphEdges
        else fun <- function(object, dg = NULL, ...) 
                standardGeneric("setGraphEdges")
        setGeneric("setGraphEdges", fun)
    }
    setMethod("setGraphEdges", signature(object = "dg.Model"), 
        function(object, dg = NULL, ...) {
            if (!is.null(dg)) 
                object@dg <- dg
            return(object)
        })
    setClass("dg.Test", representation(deviance = "numeric", 
                                       df = "numeric", p = "numeric"))
    setMethod("setSlots", "dg.Test", 
        function(object, arguments) {
        for (i in seq(along = arguments)) {
            name <- names(arguments)[i]
            if (is.element(name, slotNames(object))) 
                slot(object, name) <- arguments[[i]]
            else message(paste("Argument '", name, "' not valid slot of '", 
                class(object), "', thus ignored.", sep = ""))
        }
        return(object)
    })
    setMethod("initialize", "dg.Test", 
        function(.Object, ...) {
        Args <- list(...)
        if (!is.element("df", names(Args)) || !is.element("deviance", 
            names(Args))) {
            Args <- (Args[!names(Args) == "df"])
            Args <- (Args[!names(Args) == "deviance"])
            .Object@df <- round(runif(1, 1, 25))
            .Object@deviance <- rchisq(1, .Object@df)
            .Object@p <- 1 - pchisq(.Object@deviance, .Object@df)
            message("Just generating a random test!!!!")
        }
        .Object <- setSlots(.Object, Args)
        return(.Object)
    })
    if (!isGeneric("asDataFrame")) {
        if (is.function("asDataFrame")) 
            fun <- asDataFrame
        else fun <- function(objectlist, setRowLabels = FALSE, 
            ...) standardGeneric("asDataFrame")
        setGeneric("asDataFrame", fun)
    }
    setMethod("asDataFrame", "dg.list", 
        function(objectlist, 
        setRowLabels = FALSE, ...) {
        .asDataFrame(objectlist, setRowLabels, ...)
    })
    if (!isGeneric("Str")) {
        if (is.function("Str")) 
            fun <- Str
        else fun <- function(object, setRowLabels = FALSE, title = "", 
            ...) standardGeneric("Str")
        setGeneric("Str", fun)
    }
    setMethod("Str", "NULL", 
        function(object, setRowLabels = FALSE, 
        title = "", ...) message(paste(title, "NULL", sep = ": ")))
    setMethod("Str", "integer", 
        function(object, setRowLabels = FALSE, 
        title = "", ...) message(paste(title, paste(object, collapse = ", "), 
        sep = ": ")))
    setMethod("Str", "numeric", 
        function(object, setRowLabels = FALSE, 
        title = "", ...) message(paste(title, paste(object, collapse = ", "), 
        sep = ": ")))
    setMethod("show", "dg.list", 
        function(object) Str(object))
    setMethod("Str", "list", 
        function(object, setRowLabels = TRUE, title = "", ...) {
        if ((length(object) > 0) && (!is.null(object[[1]]))) 
            if ((extends(class(object[[1]]), "dg.Node")) || 
                (extends(class(object[[1]]), "dg.Vertex")) || 
                (extends(class(object[[1]]), "dg.VertexEdge"))) {
                if (class(object) == "list") 
                  message("<<<<< List-object not 'dg.list' !!! >>>>>")
                if ((extends(class(object[[1]]), "dg.Block")) && 
                  (length(object) == 2) && (is.list(object[[2]]))) {
                  .StrBlockTree(object, title)
                }
                else {
                  Y <- asDataFrame(object, setRowLabels = setRowLabels, ...)
                  if ("label" %in% dimnames(Y)[[2]]) {
                    labels <- t(Y["label"])
                    if (length(unique(match(labels, labels))) < 
                      length(labels)) 
                      labels <- paste(title, labels, 1:length(labels), 
                        sep = ".")
                    else labels <- paste(title, labels, sep = ".")
                    dimnames(Y)[[1]] <- labels
                  }
                  print(Y)
                }
            }
            else NextMethod("str", object, ...)
    })
    setMethod("show", "dg.graphedges", 
        function(object) Str(object))
    setMethod("Str", "dg.graphedges", 
        function(object, setRowLabels = FALSE, 
        title = "", m = 0, ...) {
        if (!is.null(object)) {
            message(object@viewType)
            message(object@oriented)
            graphEdges <- is.element("graphEdges", names(list(...)))
            if (!graphEdges) {
                if (is.element("vertexList", slotNames(object))) 
                  Str(object@vertexList, title = "vertices", ...)
                if (is.element("blockList", slotNames(object))) 
                  Str(object@blockList, title = "blocks", ...)
            }
            Str(object@visibleVertices, title = "visibleVertices", ...)
            Str(object@visibleBlocks, title = "visibleBlocks", ...)
            Str(object@extraList, title = "extraList", ...)
            Str(object@extraEdgeList, title = "extraEdgeList", ...)
            Str(object@edgeList, title = "edgeList", ...)
            Str(object@blockEdgeList, title = "blockEdgeList", ...)
            Str(object@factorVertexList, title = "factorVertexList", ...)
            Str(object@factorEdgeList, title = "factorEdgeList", ...)
        }
    })
    setMethod("show", "DynamicGraphView", 
        function(object) Str(object))
    setMethod("Str", "DynamicGraphView", 
        function(object, setRowLabels = FALSE, 
        title = "", m = 0, ...) {
        message(paste(rep("-", 80)))
        sub.title <- paste("<<", object@label, " | ", m, object@index, 
            ">>", collapse = " ")
        nc <- 78 - nchar(title) - nchar(sub.title)
        if (nc < 0) 
            nc <- 1
        cat(paste(title, paste(rep(" ", nc), collapse = ""), 
            sub.title, "\n", collapse = ""))
        if (is.element("vertexList", slotNames(object@dg))) 
            object@dg@vertexList <- new("dg.VertexList")
        if (is.element("blockList", slotNames(object@dg))) 
            object@dg@blockList <- new("dg.BlockList")
        Str(object@dg, title = "dg", ...)
    })
    setMethod("Str", "dg.Model", 
        function(object, setRowLabels = FALSE, 
        title = "", ...) {
        message(object@name)
    })
    setMethod("show", "DynamicGraphModel",
        function(object) Str(object))
    setMethod("Str", "DynamicGraphModel", 
        function(object, setRowLabels = FALSE, 
        title = "", ...) {
        if (!is.null(object)) {
            message(paste(rep("=", 80)))
            sub.title <- paste("<<", object@label, " | ", object@index, 
                ">>", collapse = " ")
            nc <- 78 - nchar(title) - nchar(sub.title)
            if (nc < 0) 
                nc <- 1
            cat(paste(title, paste(rep(" ", nc), collapse = ""), 
                sub.title, "\n", collapse = ""))
            if (hasMethod("Str", class(object@model[[1]]))) 
                Str(object@model[[1]])
            for (i in 1:length(object@graphs)) Str(object@graphs[[i]], 
                m = object@index, title = paste(title, "; Graph: ", 
                  i, sep = " "), ...)
        }
    })
    setMethod("show", "DynamicGraph", 
        function(object) Str(object))
    setMethod("Str", "DynamicGraph", 
        function(object, setRowLabels = FALSE, 
        title = "", ...) {
        message(paste(rep("#", 80)))
        message(paste("<<", object@label, ">>", sep = ""))
        Str(object@vertices, title = "vertices", ...)
        Str(object@blocks, title = "blocks", ...)
        for (i in 1:length(object@models)) Str(object@models[[i]], 
            title = paste("Model: ", i, sep = " "), ...)
        message(paste(rep("#", 80)))
    })
    if (!isGeneric("testEdge")) {
        if (is.function("testEdge")) 
            fun <- testEdge
        else fun <- function(object, action, name.1, name.2, 
            ...) standardGeneric("testEdge")
        setGeneric("testEdge", fun)
    }
    setMethod("testEdge", signature(object = "dg.Model"), 
        function(object, 
        action, name.1, name.2, ...) {
        args <- list(...)
        from.type <- args$from.type
        to.type <- args$to.type
        "f" <- function(type) if (is.null(type)) 
            ""
        else paste("(", type, ")")
        message(paste("Should return an object with the edge from", 
            name.1, f(from.type), "to", name.2, f(to.type), 
            "deleted from the argument object"))
        return(new("dg.Test"))
    })
    if (!isGeneric("modifyModel")) {
        if (is.function("modifyModel")) 
            fun <- modifyModel
        else fun <- function(object, action, name, name.1, name.2, ...) 
            standardGeneric("modifyModel")
        setGeneric("modifyModel", fun)
    }
    setMethod("modifyModel", signature(object = "dg.Model"), 
        function(object, action, name, name.1, name.2, ...) {
            args <- list(...)
            localArguments <- args$Arguments
            Edges <- args$newEdges$vertexEdges
            Vertices <- localArguments$vertexList
            viewType <- "Simple"
            DoFactors <- FALSE
            if (!is.null(args$Arguments) && 
                !is.null(args$Arguments$factorVertexList) && 
                (length(args$Arguments$factorVertexList) > 0) && 
                !is.null(args$Arguments$vertexList)) 
                DoFactors <- TRUE
            if (DoFactors) 
                viewType <- "Factor"
            FactorVertices <- new("dg.FactorVertexList")
            FactorEdges <- new("dg.FactorEdgeList")
            BlockEdges <- new("dg.BlockEdgeList")
            VisibleVertices <- localArguments$visibleVertices
            VisibleBlocks <- localArguments$visibleBlocks
            ExtraVertices <- new("dg.VertexList")
            ExtraEdges <- new("dg.ExtraEdgeList")
            "f" <- function(type) if (is.null(type)) 
                ""
            else paste("(", type, ")")
            "g" <- function(type) if (is.null(type)) 
                ""
            else type
            if (action == "dropEdge") {
                message(paste("Should return an object with the edge from", 
                  name.1, f(args$from.type), "to", name.2, f(args$to.type), 
                  "deleted from the argument object"))
                if ((g(args$from.type) == "Factor") || 
                    (g(args$from.type) == "Factor")) 
                  return(NULL)
            }
            else if (action == "addEdge") {
                message(paste("Should return an object with the edge from", 
                  name.1, f(args$from.type), "to", name.2, f(args$to.type), 
                  "added to the argument object"))
                if ((g(args$from.type) == "Factor") || 
                    (g(args$from.type) == "Factor")) 
                  return(NULL)
            }
            else if (action == "dropVertex") {
                message(paste("Should return an object with the vertex", 
                  name, f(args$type), "deleted from the argument object"))
                if ((g(args$type) == "Factor")) 
                  return(NULL)
                VisibleVertices <- VisibleVertices[VisibleVertices != 
                  args$index]
                if (DoFactors && (args$index > 0)) {
                  x <- (localArguments$factorVertexList)
                  factors <- lapply(x, function(i) i@vertex.indices)
                  types <- lapply(x, function(i) class(i))
                  factors <- lapply(factors, function(x) {
                    y <- x[x != args$index]
                    if (length(y) > 0) 
                      return(y)
                    else return(NULL)
                  })
                  if (!is.null(factors)) {
                    types <- types[unlist(lapply(factors, 
                                                 function(i) !is.null(i)))]
                    factors <- .removeNull(factors)
                  }
                  if (!is.null(factors)) {
                    subset <- function(x) lapply(x, 
                      function(a) 
                        any(unlist(lapply(x, 
                      function(A) 
                        all(!is.na(match(a, A))) && (length(a) < length(A))))))
                    s <- subset(factors)
                    types <- types[!unlist(s)]
                    factors <- factors[!unlist(s)]
                    if (!(is.null(factors))) {
                      result <- returnFactorVerticesAndEdges(
                        localArguments$vertexList, 
                        factors, types, factorClasses = validFactorClasses())
                      FactorVertices <- result$FactorVertices
                      FactorEdges <- result$FactorEdges
                    }
                  }
                  else {
                    DoFactors <- FALSE
                    FactorVertices <- new("dg.FactorVertexList")
                    FactorEdges <- new("dg.FactorEdgeList")
                  }
                }
            }
            else if (action == "addVertex") {
                VisibleVertices <- c(VisibleVertices, args$index)
                message(paste("Should return an object with the vertex", 
                  name, f(args$type), args$index, 
                  "added to the argument object"))
                if (DoFactors && (args$index > 0)) {
                  x <- (localArguments$factorVertexList)
                  factors <- lapply(x, function(i) i@vertex.indices)
                  types <- lapply(x, function(i) class(i))
                  if (!is.null(factors)) 
                    factors <- .removeNull(factors)
                  if (is.null(factors)) {
                    factors <- list(args$index)
                    types <- validFactorClasses()[1, 1]
                  }
                  else {
                    n <- length(types)
                    factors <- append(factors, list(args$index))
                    types <- append(types, types[n])
                  }
                  if (!(is.null(factors))) {
                    result <- returnFactorVerticesAndEdges(
                      localArguments$vertexList, 
                      factors, types, factorClasses = validFactorClasses())
                    FactorVertices <- result$FactorVertices
                    FactorEdges <- result$FactorEdges
                  }
                }
            }
            if (is.null(FactorVertices) && DoFactors && !is.null(Edges)) {
                factors <- .cliquesFromEdges(Edges, Vertices, 
                  VisibleVertices)
                if (is.null(factors) || (length(factors) == 0)) {
                  FactorVertices <- new("dg.FactorVertexList")
                  FactorEdges <- new("dg.FactorEdgeList")
                }
                else {
                  result <- returnFactorVerticesAndEdges(Vertices, 
                    factors)
                  FactorVertices <- result$FactorVertices
                  FactorEdges <- result$FactorEdges
                }
            }
            dg <- new("dg.graphedges", edgeList = Edges, 
                viewType = viewType, 
                blockEdgeList = BlockEdges, 
                factorVertexList = FactorVertices, 
                factorEdgeList = FactorEdges, 
                visibleVertices = VisibleVertices, 
                visibleBlocks = VisibleBlocks, 
                extraList = ExtraVertices, 
                extraEdgeList = ExtraEdges)
            if (.IsEmpty(FactorEdges) && (viewType == "Factor")) {
                object <- setGraphEdges(object, dg = dg)
                graphContent <- graphEdges(object, viewType = viewType, 
                  Arguments = localArguments)
                dg <- graphContent
            }
            return(list(object = object, dg = dg))
        })
    setMethod("label", "dg.Test", 
        function(object) format(object@p, digits = 4))
    setMethod("width", "dg.Test", 
        function(object) round(2 + 5 * (1 - object@p)))
