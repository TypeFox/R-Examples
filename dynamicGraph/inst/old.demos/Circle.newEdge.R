
# Test with new edge class:

# demo("startup.0", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/startup.0.R", sep = "/"))

setClass("NewEdge", contains = c("dg.Node", "dg.Edge", "dg.VertexEdge"))

myEdgeClasses <- rbind(validEdgeClasses(), 
                       NewEdge = c("NewEdge", "NewEdge"))

setMethod("draw", "NewEdge",
          function(object, canvas, position,
                   x = lapply(position, function(e) e[1]),
                   y = lapply(position, function(e) e[2]),
                   stratum = as.vector(rep(0, length(position)),
                                       mode = "list"), 
                   w = 2, color = "green", background = "white")
          {
            f <- function(i, j) {
              dash <- "."
              arrowhead <- "both"
              l <- function(xi, yi, xj, yj)
                tkcreate(canvas, "line", xi, yi, xj, yj, width = w, 
                         arrow = arrowhead, dash = dash,
                         # arrowshape = as.list(c(2, 5, 3) * w),
                         fill = color(object), activefill = "DarkSlateGray")
              lines <- list(l(x[[i]], y[[i]], x[[j]], y[[j]]))
              label.position <- (position[[i]] + position[[j]]) / 2
              pos <- label.position + rep(0, length(label.position))
              label <- tkcreate(canvas, "text", pos[1], pos[2],
                                text = object@label, anchor = "nw", 
				font = "8x16", activefill = "DarkSlateGray")
              tags <- NULL
              x. <- mean(unlist(x))
              y. <- mean(unlist(y))
              s <- 4 * w * sqrt(4 / pi)
              p <- tkcreate(canvas, "rectangle",
			    x. - s, y. - s, x. + s, y. + s, 
                            fill = color(object), activefill = "SeaGreen")
              tags <- list(p)
              return(list(lines = lines, tags = tags,
                          from = object@vertex.indices[i],
                          to = object@vertex.indices[j],
                          label = label, label.position = label.position))
            }
            result <- NULL
            edge <- object@vertex.indices
            m <- length(edge)
            for (j in seq(along = edge))
              if (j < length(edge))
                for (k in (j+1):length(edge))
                    result <- append(result, list(f(j, k)))
            return(result) 
  })

setMethod("addToPopups", "NewEdge",
          function(object, type, nodePopupMenu, i,
			   updateArguments, Args, ...)
          {
               tkadd(nodePopupMenu, "command",
                     label = paste(" --- This is a my new vertex!"),
                     command = function() { print(name(object))})
          })

# Why are these 2 * 7 methods not avaliable from "dg.VertexEdge" ?

# setMethod("color", "NewEdge", function(object) object@color)
# setReplaceMethod("color", "NewEdge",
#                  function(x, value) {x@color <- value; x} )

# setMethod("label", "NewEdge", function(object) object@label)
# setReplaceMethod("label", "NewEdge",
#                  function(x, value) {x@label <- value; x} )

# setMethod("name", "NewEdge", function(object) object@label)
# setReplaceMethod("name", "NewEdge",
#                  function(x, value) {x@label <- value; x} )

# setMethod("labelPosition", "NewEdge",
#           function(object) object@label.position)
# setReplaceMethod("labelPosition", "NewEdge",
#                  function(x, value) {x@label.position <- value; x} )

# setMethod("nodeIndices", "NewEdge", function(object) object@vertex.indices)
# setReplaceMethod("nodeIndices", "NewEdge",
#                  function(x, value) {x@vertex.indices <- value; x} )

# setMethod("width", "NewEdge", function(object) object@width)
# setReplaceMethod("width", "NewEdge",
#                  function(x, value) {x@width <- value; x} )
  
# setMethod("dash", "NewEdge", function(object) object@dash)
# setReplaceMethod("dash", "NewEdge",
#   function(x, value) { .dashReplaceMethod(x, value) } )
# 
 
# setMethod("propertyDialog", "NewEdge",
#           function(object, classes = NULL, title = class(object),
#                    sub.title = label(object), name.object = name(object),
#                    okReturn = TRUE,
#                    fixedSlots = NULL, difficultSlots = NULL,
#                    top = NULL, entryWidth = 20, do.grab = FALSE) {
#   .propertyDialog(object, classes = classes, title = title,
#                   sub.title = sub.title, name.object = name.object,
#                   okReturn = okReturn, 
#                   fixedSlots = fixedSlots, difficultSlots = difficultSlots,
#                   top = top, entryWidth = entryWidth, do.grab = do.grab)
#   })


V.Types <- c("Discrete", "Ordinal", "Discrete",
             "Continuous", "Discrete", "Continuous")

V.Names <- c("Sex", "Age", "Eye", "FEV", "Hair", "Shosize")
V.Labels <- paste(V.Names, 1:6, sep ="/")

From <- c(1, 2, 3, 4, 5, 6, 3)
To   <- c(2, 3, 4, 5, 6, 1, 6)

Z <- DynamicGraph(V.Names, V.Types, From, To, texts = c("Gryf", "gaf"),
                  edge.types = c("NewEdge",
                                 "VertexEdge",
                                 "Dashed",
                                 "Dotted",
                                 "DoubleArrow",
                                 "DoubleConnected",
                                 "TripleConnected"),
                  labels = V.Labels, object = Object, UserMenus = Menus, 
                  updateEdgeLabels = FALSE,
                  edgeColor = "green", vertexColor = "blue", 
                  debug.strata = debug.strata, debug.edges = debug.edges, 
                  debug.position = debug.position, debug.update = debug.update,
                  edgeClasses = myEdgeClasses)
