your.DrawModel <- function(object, slave = FALSE, viewType = "Simple", ...) {
    args <- list(...)
    Args <- args$Arguments

    # Here you should make your new model (this is just a copy):

    Object <- object
    title <- Object@name

    # and compute edges (here 'NULL' if the model not has been updated):

    Edges <- graphComponents(Object, viewType, Arguments = Args)

    EdgeList       <- Edges$vertexEdges
    ExtraVertices  <- Edges$extraVertices
    ExtraEdges     <- Edges$extraEdges
    FactorVertices <- Edges$factorVertices
    FactorEdges    <- Edges$factorEdges
    BlockEdges     <- Edges$blockEdges
    visualVertices <- Edges$visualVertices
    visualBlocks   <- Edges$visualBlocks

    if (slave) {
      # Drawing ''an other model'' in a new window:
      DynamicGraph(addModel = TRUE,                   # <-
                   frameModels = Args$frameModels, 
                   frameViews = NULL,                 # <- not used here
                   graphWindow = NULL,                # <- not used here
                   edgeList = EdgeList,
                   object = Object, 
                   extraList = ExtraVertices, 
                   extraEdgeList = ExtraEdges, 
                   factorVertexList = FactorVertices, 
                   factorEdgeList = FactorEdges, 
                   blockEdgeList = BlockEdges, 
                   visualVertices = visualVertices,
                   visualBlocks = visualBlocks,
                   title = title, 
                   Arguments = Args)
    }
    else {
      # Overwriting with ''an other model'' in same view:
      DynamicGraph(overwrite = TRUE,                 # <-
                   addModel = TRUE,                  # <-
                   frameModels = Args$frameModels, 
                   frameViews = Args$frameViews, 
                   graphWindow = Args$graphWindow,   # <-
                   edgeList = EdgeList, 
                   object = Object, 
                   extraList = ExtraVertices, 
                   extraEdgeList = ExtraEdges, 
                   factorVertexList = FactorVertices, 
                   factorEdgeList = FactorEdges, 
                   blockEdgeList = BlockEdges, 
                   visualVertices = visualVertices,
                   visualBlocks = visualBlocks,
                   title = "Now used!", 
                   Arguments = Args)
    }
}

your.LabelAllEdges <- function(object, slave = FALSE, ...) 
 {
  args <- list(...)
  Args <- args$Arguments

  getNodeName <- function(index, type)
    if (type == "Vertex")
      name(Args$vertexList[[index]])
    else if (type == "Factor")
      name(Args$factorVertexList[[abs(index)]])
    else if (type == "Extra")
      name(Args$extraList[[abs(index)]])
    else if (type == "Block")
      label(Args$blockList[[abs(index)]])
    else
      NULL

  visitEdges <- function(edges) {
    for (i in seq(along = edges)) {
      vertices <- nodeIndicesOfEdge(edges[[i]])
      types    <- nodeTypesOfEdge(edges[[i]])

      name.f <- getNodeName(vertices[1], types[1])
      name.t <- getNodeName(vertices[2], types[2])

      R <- testEdge(object, action = "remove",
                    name.1 = name.f, name.2 = name.t,
                    from = vertices[1], to = vertices[2],
                    from.type = types[1], to.type = types[2],
                    edge.index = i, force = force, Arguments = Args)

      if (!is.null(R)) {
        if (TRUE || (hasMethod("label", class(R))))
          label(edges[[i]]) <- label(R)
        if (TRUE || (hasMethod("width", class(R))))
          width(edges[[i]]) <- width(R)
      }
    }
    return(edges)
  }

  edgeList       <- visitEdges(Args$edgeList)
  factorEdgeList <- visitEdges(Args$factorEdgeList)
  blockEdgeList  <- visitEdges(Args$blockEdgeList)

  # if (length(formals((Args$frameView@redrawView))) > 0) {
    if (slave) {
      # Adding an other view of the same model:
      DynamicGraph(addModel = TRUE,                 # <-
                   # addView = TRUE,                # <-
                   frameModels = Args$frameModels, 
                   frameViews = Args$frameViews, 
                   graphWindow = NULL,              # <- not used here
                   edgeList = edgeList, 
                   factorEdgeList = factorEdgeList, 
                   blockEdgeList = blockEdgeList, 
                   title = "A slave window", 
                   Arguments = Args)
    }
    else {
      # Overwriting with an other view of the same model:
      DynamicGraph(overwrite = TRUE,                # <-
                   addModel = TRUE,                 # <-
                   # addView = TRUE,                # <-
                   frameModels = Args$frameModels, 
                   frameViews = Args$frameViews, 
                   graphWindow = Args$graphWindow,  # <-
                   edgeList = edgeList, 
                   factorEdgeList = factorEdgeList, 
                   blockEdgeList = blockEdgeList, 
                   title = "NoW used!", 
                   Arguments = Args)
    } 
 }

 test.function <- function(...) print("Test.Function")
 test.function <- function(...) print(list(...)$Arguments$object@name)

 Menus <- 
 list(MainUser = 
      list(label = "Transformation by 'prcomp' on position of \"vertices\", and redraw",
           command = function(object, ...) {
             Args <- list(...)$Arguments
             transformation <- t(prcomp(Positions(Args$vertexList))$rotation)
             Args$redrawView(graphWindow = Args$graphWindow,
                             transformation = transformation, Arguments = Args)
             }),
      MainUser = 
      list(label = "Position of \"vertices\" by 'cmdscale', and redraw",
           command = function(object, ...) {
             Args <- list(...)$Arguments
             Vertices <- Args$vertexList
             Edges <- Args$edgeList
             positions <- Positions(Args$vertexList)
             N <- dim(positions)[2]
             e <- NodeIndices(Edges)
             n <- Names(Vertices)
             X <- matrix(rep(-1, length(n)^2), ncol = length(n))
             for (i in 1:length(e)) {
               suppressWarnings(w <- as.numeric(names(e)[i]))
               if (is.na(w)) w <- .5
               X[e[[i]][1], e[[i]][2]] <- w
               X[e[[i]][2], e[[i]][1]] <- w
             }
             dimnames(X) <- list(n, n)
             d <- 1.25
             X[X==-1] <- d
             X <- X - d * diag(length(n))
             mdsX <- cmdscale(X, k = N, add = TRUE, eig = TRUE, x.ret = TRUE)
             # mdsX <- isoMDS(X, k = N)
             M <- max(abs(mdsX$points))
             Positions(Args$vertexList) <<- mdsX$points / M * 45
             Args$redrawView(graphWindow = Args$graphWindow, 
                             # Positions = Positions(Args$vertexList), 
                             vertexList = Args$vertexList, Arguments = Args)
             }),
      MainUser = 
      list(label = "Position of \"vertices\"",
           command = function(object, ...) 
             print(Positions(list(...)$Arguments$vertexList))),
      MainUser = 
      list(label = "Label all edges, in this window",
           command = function(object, ...) 
                       your.LabelAllEdges(object, slave = FALSE, ...)),
      MainUser = 
      list(label = "Label all edges, in slave window",
           command = function(object, ...) 
                       your.LabelAllEdges(object, slave = TRUE, ...)),
      MainUser = 
      list(label = "Draw model, in this window",
           command = function(object, ...) 
                       your.DrawModel(object, slave = FALSE, ...)),
      MainUser = 
      list(label = "Draw model, in slave window",
           command = function(object, ...) 
                       your.DrawModel(object, slave = TRUE, ...)),
      MainUser = 
      list(label = "Call of function 'modalDialog', result on 'title' at top",
           command = function(object, ...) 
           {
             Args <- list(...)$Arguments
             ReturnVal <- modalDialog("Test modalDialog Entry",
                                      "Enter name", Args$control$title, 
                                      top = Args$top)
             print(ReturnVal)
             if (ReturnVal == "ID_CANCEL")
               return()
             tktitle(Args$top) <- ReturnVal
           }
          ),
      MainUser = 
      list(label = "Call of function 'test.function', result on 'viewLabel' at bottom",
           command = function(object, ...)
           {
             Args <- list(...)$Arguments
             tkconfigure(Args$viewLabel, 
                         text = paste(Args$viewType, " | ", 
                                      test.function(...))) } ),
      Vertex = 
      list(label = "Test of user popup menu for vertices: Label",
           command = function(object, name, ...) 
           {
             # print(name)
             args <- list(...)
             # print(names(args))
             # print(c(args$type))
             # print(c(args$index))
             Args <- args$Arguments
             # print(names(Args))
             # str(Args$selectedNodes)
             # str(Args$selectedEdges)
             print(Args$vertexList[[args$index]]@label)
           }
          ),
      Edge = 
      list(label = "Test of user popup menu for edges: Class",
           command = function(object, name1, name2, ...) 
           {
             args <- list(...)
             # print(c(name1, name2))
             # print(c(args$edge.index, args$which.edge, args$from, args$to))
             # print(c(args$from.type, args$to.type, args$edge.type))
             Args <- list(...)$Arguments
             # print(names(Args))
             # str(Args$selectedNodes)
             # str(Args$selectedEdges)
             ReturnVal <- selectDialog("Test selectDialog Entry", "Select name", 
                                       Args$control$edgeClasses[,1], top = Args$top)
             print(ReturnVal)
             if (ReturnVal == "ID_CANCEL")
               return()
             if ((args$from > 0) && (args$to > 0)) {
               edgeList <- Args$edgeList
               class(edgeList[[args$edge.index]]) <- 
                                              Args$control$edgeClasses[ReturnVal, 2]
               # vertexEdges (Args$object) <<- edgeList # Not working !!!
               Args$redrawView(graphWindow = Args$graphWindow,
                               edgeList = edgeList, title = "Not used!", 
                               width = NULL, height = NULL, Arguments = Args)
             }
           }
          ),
      ClosedBlock = 
      list(label = "Test of user popup menu for blocks",
           command = function(object, name, ...) 
           {
             print(name)
             print(c(list(...)$index))
           }
          ),
     )
