
# Edit the following to meet your needs:
#
# - Change the name "your.Model"
#
# - Work out how the get names, types and edges from the model object.
#
# - At "message", insert the relevant code for testing and modifying the model.
#
# - The slots visibleVertices, visibleBlocks, extraVertices, graphEdges, 
#    blockEdges, factorVertices, factorEdges should be eliminated, and you
#    should in "graphComponents" return relevant lists.
#

setClass("your.Model", 
         representation(name             = "character",
                        visibleVertices  = "numeric",
                        visibleBlocks    = "numeric",
                        extraVertices    = "dg.VertexList",
                        vertexEdges      = "dg.VertexEdgeList",
                        blockEdges       = "dg.BlockEdgeList",
                        factorVertices   = "dg.FactorVertexList",
                        factorEdges      = "dg.FactorEdgeList",
                        extraEdges       = "dg.ExtraEdgeList"))

"newYourModelObject"<-
  function(name)
  {
    result <- new("your.Model", name = name,
                  extraVertices  = .emptyDgList("dg.VertexList"),
                  vertexEdges    = .emptyDgList("dg.VertexEdgeList"),
                  blockEdges     = .emptyDgList("dg.BlockEdgeList"),
                  factorVertices = .emptyDgList("dg.FactorVertexList"),
                  factorEdges    = .emptyDgList("dg.FactorEdgeList"),
                  extraEdges     = .emptyDgList("dg.ExtraEdgeList"))
    return(result)
  }

if (!isGeneric("graphComponents")) {
  if (is.function("graphComponents"))
    fun <- graphComponents
  else
    fun <- function(object, viewType = NULL, ...)
    standardGeneric("graphComponents")
  setGeneric("graphComponents", fun)
}

setMethod("graphComponents", "your.Model",
          function(object, viewType = NULL, ...)
          { # print(viewType); print ("graphComponents")
            args <- list(...)
            Args <- args$Arguments
            Edges <- object@vertexEdges
            Vertices <- Args$vertexList
            VisibleVertices <- object@visibleVertices
            if (viewType == "Factor") {
              # require(ggm)
              # e <- NodeIndices(Edges)
              # if (length(e) > 0) {
              #   e <- lapply(e, function(egde) if (sum(abs(egde)) > 0) egde)
              #   e <- .removeNull(e)
              # } else
              #   e <- NULL
              # factors <- NULL
              # if (length(e) < 2) {
              #   if (length(e) == 1)
              #     factors <- append(e, as.list(VisibleVertices))
              #   else if (length(VisibleVertices) > 0)
              #     factors <- as.list(VisibleVertices)
              # } else {
              #   n <- Names(Vertices)
              #   X <- matrix(rep(0, length(n)^2), ncol = length(n))
              #   lapply(e, function(i) { X[i[1], i[2]] <<- 1 ;
              #             X[i[2], i[1]] <<- 1 } )
              #   dimnames(X) <- list(n, n)
              #   X <- X[VisibleVertices, VisibleVertices]
              #   # print ("graphComponents: ")
              #   # print(X)
              #   factors <- cliques(X)
              # }
              factors <- .cliquesFromEdges(Edges, Vertices, VisibleVertices)
              # print(factors)
              if (is.null(factors) || (length(factors) == 0)) {
                FactorVertices  <- .emptyDgList("dg.FactorVertexList")
                FactorEdges     <- .emptyDgList("dg.FactorEdgeList")
              } else {
                result <- returnFactorVerticesAndEdges(Vertices, factors)
                FactorVertices  <- result$FactorVertices
                FactorEdges     <- result$FactorEdges
              }
              list(vertexEdges     = object@vertexEdges, 
                   blockEdges      = object@blockEdges, 
                   factorVertices  = FactorVertices,
                   factorEdges     = FactorEdges,
                   visibleVertices = object@visibleVertices, 
                   visibleBlocks   = object@visibleBlocks, 
                   extraVertices   = object@extraVertices, 
                   extraEdges      = object@extraEdges)
            } else if (viewType == "Moral") {
              message("Moral view not implemented; ")
              list(vertexEdges     = object@vertexEdges, 
                   blockEdges      = .emptyDgList("dg.BlockEdgeList"),
                   factorVertices  = .emptyDgList("dg.FactorVertexList"),
                   factorEdges     = .emptyDgList("dg.FactorEdgeList"),
                   visibleVertices = object@visibleVertices, 
                   visibleBlocks   = numeric(), 
                   extraVertices   = object@extraVertices, 
                   extraEdges      = object@extraEdges)
            } else if (viewType == "Essential") {
              message("Essential view not implemented; ")
              list(vertexEdges      = object@vertexEdges, 
                   blockEdges      = .emptyDgList("dg.BlockEdgeList"),
                   factorVertices  = .emptyDgList("dg.FactorVertexList"),
                   factorEdges     = .emptyDgList("dg.FactorEdgeList"),
                   visibleVertices = object@visibleVertices, 
                   visibleBlocks   = numeric(), 
                   extraVertices   = object@extraVertices, 
                   extraEdges      = object@extraEdges)
            } else if (viewType == "Simple") {
              list(vertexEdges     = object@vertexEdges, 
                   blockEdges      = object@blockEdges, 
                   factorVertices  = .emptyDgList("dg.FactorVertexList"),
                   factorEdges     = .emptyDgList("dg.FactorEdgeList"),
                   visibleVertices = object@visibleVertices, 
                   visibleBlocks   = object@visibleBlocks, 
                   extraVertices   = object@extraVertices, 
                   extraEdges      = object@extraEdges)
            } else 
              message("View type not implemented; ")
          })

if (!isGeneric("setGraphComponents")) {
  if (is.function("setGraphComponents"))
    fun <- setGraphComponents
  else
    fun <- function(object, viewType = NULL,
                    visibleVertices = NULL,
                    visibleBlocks   = NULL,
                    extraVertices   = NULL,
                    vertexEdges     = NULL,
                    blockEdges      = NULL,
                    factorVertices  = NULL,
                    factorEdges     = NULL,
                    extraEdges      = NULL, ...)
      standardGeneric("setGraphComponents")
  setGeneric("setGraphComponents", fun)
}

setMethod("setGraphComponents", "your.Model",
          function(object, viewType = NULL,
                   visibleVertices = NULL,
                   visibleBlocks   = NULL,
                   extraVertices   = NULL,
                   vertexEdges     = NULL,
                   blockEdges      = NULL,
                   factorVertices  = NULL,
                   factorEdges     = NULL,
                   extraEdges      = NULL, ...)
 {
    if (!is.null(visibleVertices)) object@visibleVertices <- visibleVertices
    if (!(viewType == "Moral"))
      if (!is.null(visibleBlocks  )) object@visibleBlocks   <- visibleBlocks
    if (!is.null(extraVertices  )) object@extraVertices   <- extraVertices
    if (!is.null(extraEdges     )) object@extraEdges      <- extraEdges
    if (!is.null(vertexEdges    )) object@vertexEdges     <- vertexEdges
    if (!is.null(blockEdges     )) object@blockEdges      <- blockEdges
    if ((viewType == "Factor")) {
      if (!is.null(factorVertices )) object@factorVertices  <- factorVertices
      if (!is.null(factorEdges    )) object@factorEdges     <- factorEdges
    }
    return(object)
 })


if (!isGeneric("dynamic.Graph")) {
  if (is.function("dynamic.Graph"))
    fun <- dynamic.Graph
  else
    fun <- function(object, ...)
  standardGeneric("dynamic.Graph")
  setGeneric("dynamic.Graph", fun)
}

setMethod("dynamic.Graph", signature(object = "your.Model"),
          function(object, ...)
  {

    Names <- Your.function.for.extracting.variable.names.from.object(
             object = object)
    Types <- Your.function.for.extracting.variable.types.from.object(
             object = object)
    Edges <- Your.function.for.extracting.variable.edges.from.object(
             object = object)

    DynamicGraph(names = Names, types = Types, 
                 from = Edges[,1], to = Edges[,2], 
                 object = object, ...)
 })


if (!isGeneric("testEdge")) {
  if (is.function("testEdge"))
    fun <- testEdge
  else
    fun <- function(object, action, name.1, name.2, ...) 
           standardGeneric("testEdge")
  setGeneric("testEdge", fun)
}

setMethod("testEdge", signature(object = "your.Model"),
          function(object, action, name.1, name.2, ...)
 {
    args <- list(...)
    from.type <- args$from.type
    to.type <- args$to.type
    f <- function(type) if(is.null(type)) "" else paste("(", type, ")")
    message(paste("Should return an object with the edge from",
                  name.1, f(from.type), "to", name.2, f(to.type),
                  "deleted from the argument object"))
    return(newYourTestObject())
 })


if (!isGeneric("modifyModel")) {
  if (is.function("modifyModel"))
    fun <- modifyModel
  else
    fun <- function(object, action, name, name.1, name.2, ...)
                    standardGeneric("modifyModel")
  setGeneric("modifyModel", fun)
}

setMethod("modifyModel", signature(object = "your.Model"),
          function(object, action, name, name.1, name.2, ...)
 {
    args <- list(...)
    Args <- args$Arguments
    Edges <- args$newEdges$vertexEdges
    Vertices <- Args$vertexList

    DoFactors <- FALSE
    if (!is.null(args$Arguments)
        && !is.null(args$Arguments$factorVertexList)
        && (length(args$Arguments$factorVertexList) > 0)
        && !is.null(args$Arguments$vertexList))
      DoFactors <- TRUE

    # str(args)
    # print(names(args))
    # print(names(args$Arguments))
    # print(args$Arguments$visibleVertices)
    # str(args$Arguments$selectedNodes)

    # if (length(args$selectedNodes) > 0)
    #   str(args$selectedNodes)
    # if (length(args$selectedEdges) > 0)
    #   str(args$selectedEdges)

    FactorVertices  <- NULL
    FactorEdges     <- NULL
    BlockEdges      <- NULL
    VisibleVertices <- Args$visibleVertices
    VisibleBlocks   <- Args$visibleBlocks
    ExtraVertices   <- NULL
    ExtraEdges      <- NULL

    f <- function(type) if (is.null(type)) "" else paste("(", type, ")")
    g <- function(type) if (is.null(type)) "" else type
    if (action == "dropEdge") {
      message(paste("Should return an object with the edge from",
                    name.1, f(args$from.type), "to", name.2, f(args$to.type),
                    "deleted from the argument object"))
      if ((g(args$from.type) == "Factor") || (g(args$from.type) == "Factor"))
        return(NULL)
    } else if (action == "addEdge") {
       message(paste("Should return an object with the edge from",
                     name.1, f(args$from.type), "to", name.2, f(args$to.type),
                     "added to the argument object"))
      if ((g(args$from.type) == "Factor") || (g(args$from.type) == "Factor"))
        return(NULL)
    } else if (action == "dropVertex")  {
       message(paste("Should return an object with the vertex", 
                     name, f(args$type),
                     "deleted from the argument object"))
      if ((g(args$type) == "Factor"))
        return(NULL)
      VisibleVertices <- VisibleVertices[VisibleVertices != args$index]
      if (DoFactors && (args$index > 0)) {
        x <- (args$Arguments$factorVertexList)
        factors <- lapply(x, function(i) i@vertex.indices)
        types   <- lapply(x, function(i) class(i))
        factors <- lapply(factors, 
                          function(x) { 
                            y <- x[x != args$index]
                            if (length(y) > 0) return(y) else return(NULL) } )


        if (!is.null(factors)) {
          types   <- types[unlist(lapply(factors, function(i) !is.null(i)))]
          factors <- .removeNull(factors)
        }
        if (!is.null(factors)) {
          subset <- function(x)
            lapply(x, function(a) 
                        any(unlist(lapply(x, 
                                          function(A) 
                                            all(!is.na(match(a, A))) &&
                                            (length(a) < length(A))))))
          s <- subset(factors)
          types   <- types[!unlist(s)]
          factors <- factors[!unlist(s)]
          if (!(is.null(factors))) {
            result <- returnFactorVerticesAndEdges(
                            args$Arguments$vertexList, factors, types, 
                            factorClasses = validFactorClasses())
            FactorVertices <- result$FactorVertices
            FactorEdges <- result$FactorEdges
          }
        } else { 
          DoFactors <- FALSE
          FactorVertices <- .emptyDgList("dg.FactorVertexList")
          FactorEdges    <- .emptyDgList("dg.FactorEdgeList")
        }
      }
    } else if (action == "addVertex") {
      VisibleVertices <- c(VisibleVertices, args$index)
      message(paste("Should return an object with the vertex", 
                    name, f(args$type), args$index, 
                    "added to the argument object"))
      if (DoFactors && (args$index > 0)) {
        x <- (args$Arguments$factorVertexList)
        factors <- lapply(x, function(i) i@vertex.indices)
        types   <- lapply(x, function(i) class(i))
        if (!is.null(factors))
          factors <- .removeNull(factors)
        if (is.null(factors)) {
          factors <- list(args$index)
          types   <- validFactorClasses()[1, 1]
        } else { 
          n <- length(types)
          factors <- append(factors, list(args$index))
          types   <- append(types, types[n])
        }
        if (!(is.null(factors))) {
          result <- returnFactorVerticesAndEdges(
                          args$Arguments$vertexList, factors, types, 
                          factorClasses = validFactorClasses())
          FactorVertices <- result$FactorVertices
          FactorEdges <- result$FactorEdges
        }
      }
    }
    if (is.null(FactorVertices) && DoFactors && !is.null(Edges)) {

      factors <- .cliquesFromEdges(Edges, Vertices, VisibleVertices)

      if (is.null(factors) || (length(factors) == 0)) {
        FactorVertices <- .emptyDgList("dg.FactorVertexList")
        FactorEdges    <- .emptyDgList("dg.FactorEdgeList")
      } else {
        result <- returnFactorVerticesAndEdges(Vertices, factors)
        FactorVertices  <- result$FactorVertices
        FactorEdges     <- result$FactorEdges
      }
    }
    return(list(object          = object,
                BlockEdges      = BlockEdges, 
                FactorVertices  = FactorVertices,
                FactorEdges     = FactorEdges,
                VisibleVertices = VisibleVertices, 
                VisibleBlocks   = VisibleBlocks, 
                ExtraVertices   = ExtraVertices,
                ExtraEdges      = ExtraEdges))
 })

setMethod("Str", "your.Model",
          function(object, setRowLabels = FALSE, title = "", ...) {
              message(object@name) })


setClass("your.Test", 
         representation(deviance = "numeric", df = "numeric", p = "numeric"))

"newYourTestObject"<-
  function(name)
  {
    df <- round(runif(1, 1, 25))
    message("Just generating a random test!!!!!")
    # warning("Just generating a random test!!!!!")
    deviance <- rchisq(1, df)
    p <- 1 - pchisq(deviance, df)
    result <- new("your.Test",
                   df = df, deviance = deviance, p = p)
    return(result)
  }

if (!isGeneric("label") && !isGeneric("label", where = 2)) {
  if (is.function("label"))
    fun <- label
  else
    fun <- function(object) standardGeneric("label")
  setGeneric("label", fun)
}

setMethod("label", "your.Test",
          function(object) format(object@p, digits = 4))

if (!isGeneric("width") && !isGeneric("width", where = 2)) {
  if (is.function("width"))
    fun <- width
  else
    fun <- function(object) standardGeneric("width")
  setGeneric("width", fun)
}

setMethod("width", "your.Test",
          function(object) round(2 + 5 * (1 - object@p)))

