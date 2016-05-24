
# Edit the following to meet your needs:
#
# - Change the name "your.Model"
#
# - Work out how the get names, types and edges from the model object.
#
# - At "message", insert the relevant code for testing and modifying the model.
#
# - The slots visibleVertices, visibleBlocks, extraVertices, edges, 
#    blockEdges, factorVertices, factorEdges should be eliminated,
#    and you should in "graphEdges" return relevant lists.
#

setClass("your.Model", 
         representation(name = "character",
                        dg   = "dg.graphedges"))

if (!isGeneric("setSlots") &&
      (length(attr(isGeneric("setSlots", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  warning("Method 'setSlots' should be defined and exported from dynamicGraph")
  if (is.function("setSlots"))
    fun <- setSlots
  else
    fun <- function(object, arguments) standardGeneric("setSlots")
  setGeneric("setSlots", fun)
}

setMethod("setSlots", "your.Model",
          function(object, arguments) { 
            for (i in seq(along = arguments)) {
              name <- names(arguments)[i]
              if (is.element(name, slotNames(object)))
                slot(object, name) <- arguments[[i]]
              else
                message(paste("Argument '", name, "' not valid slot of '", 
                              class(object), "', thus ignored.",
                              sep = "")) }
            return(object)
          })

setMethod("initialize", "your.Model",
          function(.Object, ...) {
              # print(c("initialize", "your.Model", class(.Object)))
              Args <- list(...)
              .Object <- setSlots(.Object, Args)
              return(.Object)
            }
          )

if (!isGeneric("graphEdges") &&
      (length(attr(isGeneric("graphEdges", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  message("Method 'graphEdges' set as generic")
  if (is.function("graphEdges"))
    fun <- graphEdges
  else
    fun <- function(object, viewType = NULL, ...)
    standardGeneric("graphEdges")
  setGeneric("graphEdges", fun)
}

setMethod("graphEdges", "your.Model",
          function(object, viewType = NULL, ...)
          { # print(viewType); print ("graphEdges")

            dots           <- list(...)
            localArguments <-      dots$Arguments
            Vertices       <-      localArguments$vertexList

            Edges           <- object@dg@edgeList
            VisibleVertices <- object@dg@visibleVertices

            if (viewType == "Factor") {
              factors <- .cliquesFromEdges(Edges, Vertices, VisibleVertices)
              # print(factors)
              if (is.null(factors) || (length(factors) == 0)) {
                FactorVertices  <- new("dg.FactorVertexList")
                FactorEdges     <- new("dg.FactorEdgeList")
              } else {
                result <- returnFactorVerticesAndEdges(Vertices, factors)
                FactorVertices  <- result$FactorVertices
                FactorEdges     <- result$FactorEdges
              }
              new("dg.graphedges", 
                   viewType         = viewType, 
                   oriented         = object@dg@oriented, 
                   edgeList         = object@dg@edgeList, 
                   blockEdgeList    = object@dg@blockEdgeList, 
                   factorVertexList = FactorVertices,
                   factorEdgeList   = FactorEdges,
                   visibleVertices  = object@dg@visibleVertices, 
                   visibleBlocks    = object@dg@visibleBlocks, 
                   extraList        = object@dg@extraList, 
                   extraEdgeList    = object@dg@extraEdgeList)
            } else if (viewType == "Moral") {
              message("Moral view not implemented; ")
              new("dg.graphedges", 
                   viewType         = viewType, 
                   oriented         = object@dg@oriented, 
                   edgeList         = object@dg@edgeList, 
                 # blockEdgeList    = new("dg.BlockEdgeList"),
                 # factorVertexList = new("dg.FactorVertexList"),
                 # factorEdgeList   = new("dg.FactorEdgeList"),
                   visibleVertices  = object@dg@visibleVertices, 
                   visibleBlocks    = numeric(), 
                   extraList        = object@dg@extraList, 
                   extraEdgeList    = object@dg@extraEdgeList)
            } else if (viewType == "Essential") {
              message("Essential view not implemented; ")
              new("dg.graphedges", 
                   viewType         = viewType, 
                   oriented         = object@dg@oriented, 
                   edgeList         = object@dg@edgeList, 
                 # blockEdgeList    = new("dg.BlockEdgeList"),
                 # factorVertexList = new("dg.FactorVertexList"),
                 # factorEdgeList   = new("dg.FactorEdgeList"),
                   visibleVertices  = object@dg@visibleVertices, 
                   visibleBlocks    = numeric(), 
                   extraList        = object@dg@extraList, 
                   extraEdgeList    = object@dg@extraEdgeList)
            } else if (viewType == "Simple") {
              new("dg.graphedges", 
                   viewType         = viewType, 
                   oriented         = object@dg@oriented, 
                   edgeList         = object@dg@edgeList, 
                   blockEdgeList    = object@dg@blockEdgeList, 
                 # factorVertexList = new("dg.FactorVertexList"),
                 # factorEdgeList   = new("dg.FactorEdgeList"),
                   visibleVertices  = object@dg@visibleVertices, 
                   visibleBlocks    = object@dg@visibleBlocks, 
                   extraList        = object@dg@extraList, 
                   extraEdgeList    = object@dg@extraEdgeList)
            } else 
              message("View type not implemented; ")
          })

if (!isGeneric("setGraphEdges") &&
      (length(attr(isGeneric("setGraphEdges", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  message("Method 'setGraphEdges' set as generic")
  if (is.function("setGraphEdges"))
    fun <- setGraphEdges
  else
    fun <- function(object, dg = NULL, ...)
      standardGeneric("setGraphEdges")
  setGeneric("setGraphEdges", fun)
}

setMethod("setGraphEdges", signature(object = "your.Model"),
          function(object, dg = NULL, ...)
 {
    if (!is.null(dg)) object@dg <- dg
    return(object)
 })

if (!isGeneric("dg") &&
      (length(attr(isGeneric("dg", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  warning("Method 'dg' should be defined and exported from dynamicGraph")
  if (is.function("dg"))
    fun <- position
  else
    fun <- function(object, 
                    # modelObject = NULL,
                    # modelObjectName = NULL,
                    # control = dg.control(...), 
                    ...) standardGeneric("dg")
  setGeneric("dg", fun)
}

setMethod("dg", signature(object = "your.Model"),
          function(object, 
                   # modelObject = NULL,
                   # modelObjectName = NULL,
                   # control = dg.control(...), 
                   ...) 
  {

    Names <- Your.function.for.extracting.variable.names.from.object(
             object = object)
    Types <- Your.function.for.extracting.variable.types.from.object(
             object = object)
    Edges <- Your.function.for.extracting.variable.edges.from.object(
             object = object)

    simpleGraph <- new("dg.simple.graph", vertex.names = Names, 
                      types = Types, # edge.list = Edges,
                      from = Edges[,1], to = Edges[,2])

    graph <- simpleGraphToGraph(simpleGraph)

    dg(graph, object = object, ...)

 })

if (!isGeneric("testEdge") &&
      (length(attr(isGeneric("testEdge", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  message("Method 'testEdge' set as generic")
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
    dots <- list(...)
    from.type <- dots$from.type
    to.type <- dots$to.type
    f <- function(type) if(is.null(type)) "" else paste("(", type, ")")
    message(paste("Should return an object with the edge from",
                  name.1, f(from.type), "to", name.2, f(to.type),
                  "deleted from the argument object"))
    return(new("your.Test"))
 })


if (!isGeneric("modifyModel") &&
      (length(attr(isGeneric("modifyModel", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  message("Method 'modifyModel' set as generic")
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
    dots           <- list(...)
    localArguments <-      dots$Arguments
    Edges          <-      dots$newEdges$vertexEdges
    Vertices       <-      localArguments$vertexList

    viewType <- "Simple"

    DoFactors <- FALSE
    if (!is.null(dots$Arguments)
        && !is.null(dots$Arguments$factorVertexList)
        && (length(dots$Arguments$factorVertexList) > 0)
        && !is.null(dots$Arguments$vertexList))
      DoFactors <- TRUE

    if (DoFactors)
      viewType <- "Factor"

    # print(names(dots))
    # str(dots)

    # print(names(localArguments))

    # print(localArguments$visibleVertices)

    # str(localArguments$selectedNodes)
    # if (length(dots$selectedNodes) > 0)
    #   str(dots$selectedNodes)

    # str(localArguments$selectedEdges)
    # if (length(dots$selectedEdges) > 0)
    #   str(dots$selectedEdges)

    FactorVertices  <- new("dg.FactorVertexList")
    FactorEdges     <- new("dg.FactorEdgeList")
    BlockEdges      <- new("dg.BlockEdgeList")
    VisibleVertices <- localArguments$visibleVertices
    VisibleBlocks   <- localArguments$visibleBlocks
    ExtraVertices   <- new("dg.VertexList")
    ExtraEdges      <- new("dg.ExtraEdgeList")

    f <- function(type) if (is.null(type)) "" else paste("(", type, ")")
    g <- function(type) if (is.null(type)) "" else type
    if (action == "dropEdge") {
      message(paste("Should return an object with the edge from",
                    name.1, f(dots$from.type), "to", name.2, f(dots$to.type),
                    "deleted from the argument object"))
      if ((g(dots$from.type) == "Factor") || (g(dots$from.type) == "Factor"))
        return(NULL)
    } else if (action == "addEdge") {
       message(paste("Should return an object with the edge from",
                     name.1, f(dots$from.type), "to", name.2, f(dots$to.type),
                     "added to the argument object"))
      if ((g(dots$from.type) == "Factor") || (g(dots$from.type) == "Factor"))
        return(NULL)
    } else if (action == "dropVertex")  {
       message(paste("Should return an object with the vertex", 
                     name, f(dots$type),
                     "deleted from the argument object"))
      if ((g(dots$type) == "Factor"))
        return(NULL)
      VisibleVertices <- VisibleVertices[VisibleVertices != dots$index]
      if (DoFactors && (dots$index > 0)) {
        x <- (localArguments$factorVertexList)
        factors <- lapply(x, function(i) i@vertex.indices)
        types   <- lapply(x, function(i) class(i))
        factors <- lapply(factors, 
                          function(x) { 
                            y <- x[x != dots$index]
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
                            localArguments$vertexList, factors, types, 
                            factorClasses = validFactorClasses())
            FactorVertices <- result$FactorVertices
            FactorEdges <- result$FactorEdges
          }
        } else { 
          DoFactors <- FALSE
          FactorVertices <- new("dg.FactorVertexList")
          FactorEdges    <- new("dg.FactorEdgeList")
        }
      }
    } else if (action == "addVertex") {
      VisibleVertices <- c(VisibleVertices, dots$index)
      message(paste("Should return an object with the vertex", 
                    name, f(dots$type), dots$index, 
                    "added to the argument object"))
      if (DoFactors && (dots$index > 0)) {
        x <- (localArguments$factorVertexList)
        factors <- lapply(x, function(i) i@vertex.indices)
        types   <- lapply(x, function(i) class(i))
        if (!is.null(factors))
          factors <- .removeNull(factors)
        if (is.null(factors)) {
          factors <- list(dots$index)
          types   <- validFactorClasses()[1, 1]
        } else { 
          n <- length(types)
          factors <- append(factors, list(dots$index))
          types   <- append(types, types[n])
        }
        if (!(is.null(factors))) {
          result <- returnFactorVerticesAndEdges(
                          localArguments$vertexList, factors, types, 
                          factorClasses = validFactorClasses())
          FactorVertices <- result$FactorVertices
          FactorEdges <- result$FactorEdges
        }
      }
    }
    if (is.null(FactorVertices) && DoFactors && !is.null(Edges)) {

      factors <- .cliquesFromEdges(Edges, Vertices, VisibleVertices)

      if (is.null(factors) || (length(factors) == 0)) {
        FactorVertices <- new("dg.FactorVertexList")
        FactorEdges    <- new("dg.FactorEdgeList")
      } else {
        result <- returnFactorVerticesAndEdges(Vertices, factors)
        FactorVertices  <- result$FactorVertices
        FactorEdges     <- result$FactorEdges
      }
    }
    dg <- new("dg.graphedges", 
              edgeList         = Edges,
              viewType         = viewType, 
            # oriented         = oriented, 
              blockEdgeList    = BlockEdges, 
              factorVertexList = FactorVertices,
              factorEdgeList   = FactorEdges,
              visibleVertices  = VisibleVertices, 
              visibleBlocks    = VisibleBlocks, 
              extraList        = ExtraVertices,
              extraEdgeList    = ExtraEdges)
    ".IsEmpty" <- function(x) {
      if (is.null(x) || (length(x) == 0) ||
          (length(x) == 1) && is.null(x[[1]]))
        return(TRUE)
      else
        return(FALSE)
      }
    if (.IsEmpty(FactorEdges) && (viewType == "Factor")) {
      object <- setGraphEdges(object, dg = dg)
      graphContent <- graphEdges(object, viewType = viewType, 
                                 Arguments = localArguments)
      dg <- graphContent
    }
    return(list(object  = object, dg = dg))
 })

if (!isGeneric("Str") &&
      (length(attr(isGeneric("Str", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  warning("Method 'Str' should be defined and exported from dynamicGraph")
  if (is.function("Str"))
    fun <- Str
  else
    fun <- function(object, setRowLabels = FALSE, title = "", ...)
             standardGeneric("Str")
  setGeneric("Str", fun)
}

setMethod("Str", "your.Model",
          function(object, setRowLabels = FALSE, title = "", ...) {
              message(object@name) })

setClass("your.Test", 
         representation(deviance = "numeric", df = "numeric", p = "numeric"))

setMethod("setSlots", "your.Test",
          function(object, arguments) { 
            for (i in seq(along = arguments)) {
              name <- names(arguments)[i]
              if (is.element(name, slotNames(object)))
                slot(object, name) <- arguments[[i]]
              else
                message(paste("Argument '", name, "' not valid slot of '", 
                              class(object), "', thus ignored.",
                              sep = "")) }
            return(object)
          })

setMethod("initialize", "your.Test",
          function(.Object, ...) {
              # print(c("initialize", "your.Test", class(.Object)))
              Args <- list(...)
              if (!is.element("df", names(Args)) || 
                  !is.element("deviance", names(Args))) {
                Args <- (Args[!names(Args) == "df"])
                Args <- (Args[!names(Args) == "deviance"])
                .Object@df <- round(runif(1, 1, 25))
                .Object@deviance <- rchisq(1, .Object@df)
                .Object@p <- 1 - pchisq(.Object@deviance, .Object@df)
                message("Just generating a random test!!!!")
              }
              .Object <- setSlots(.Object, Args)
              return(.Object)
            }
          )


if (!isGeneric("label") &&
      (length(attr(isGeneric("label", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  warning("Method 'label' should be defined and exported from dynamicGraph")
  if (is.function("label"))
    fun <- label
  else
    fun <- function(object) standardGeneric("label")
  setGeneric("label", fun)
}

setMethod("label", "your.Test",
          function(object) format(object@p, digits = 4))

if (!isGeneric("width") &&
      (length(attr(isGeneric("width", getName = TRUE),
                             "package") == "dynamicGraph") > 0)) {
  warning("Method 'width' should be defined and exported from dynamicGraph")
  if (is.function("width"))
    fun <- width
  else
    fun <- function(object) standardGeneric("width")
  setGeneric("width", fun)
}

setMethod("width", "your.Test",
          function(object) round(2 + 5 * (1 - object@p)))

