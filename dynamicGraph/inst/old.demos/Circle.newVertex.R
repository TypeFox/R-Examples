
# Test with new vertex class:

# demo("startup.0", package = "dynamicGraph", verbose = FALSE)

source(paste(system.file(package = "dynamicGraph"), 
             "Demo/startup.0.R", sep = "/"))

setClass("NewVertex", contains = c("dg.Node", "dg.Vertex"),
         representation(my.text   = "character",
                        my.number = "numeric"), 
         prototype(my.text    = "",
                   my.number  = 2))

myVertexClasses <- rbind(validVertexClasses(), 
                         NewVertex = c("NewVertex", "NewVertex"))

setMethod("draw", "NewVertex",
          function(object, canvas, position,
                   x = position[1], y = position[2], stratum = 0,
                   w = 2, color = "green", background = "white")
          {
            s <- w * sqrt(4 / pi) / 2
            p1 <- tkcreate(canvas, "oval",
                           x - s - s, y - s,
                           x + s - s, y + s,
                           fill = color(object), activefill = "IndianRed")
            p2 <- tkcreate(canvas, "oval",
                           x - s + s, y - s,
                           x + s + s, y + s, 
                           fill = color(object), activefill = "IndianRed")
            p3 <- tkcreate(canvas, "oval",
                           x - s, y - s - s,
                           x + s, y + s - s, 
                           fill = color(object), activefill = "IndianRed")
            p4 <- tkcreate(canvas, "poly", 
                           x - 1.5 * s, y + 3 * s,
                           x + 1.5 * s, y + 3 * s, 
                           x, y, 
                           fill = color(object), activefill = "SteelBlue")
            return(list(dynamic = list(p1, p2, p3, p4), fixed = NULL)) })

setMethod("addToPopups", "NewVertex",
          function(object, type, nodePopupMenu, i,
			   updateArguments, Args, ...)
          {
               tkadd(nodePopupMenu, "command",
                     label = paste(" --- This is a my new vertex!"),
                     command = function() { print(name(object))})
          })

if (!isGeneric("my.text")) {
  if (is.function("my.text"))
    fun <- my.text
  else
    fun <- function(object) standardGeneric("my.text")
  setGeneric("my.text", fun)
}
setGeneric("my.text<-",
           function(x, value) standardGeneric("my.text<-"))

setMethod("my.text", "NewVertex",
          function(object) object@my.text)
setReplaceMethod("my.text", "NewVertex",
                 function(x, value) {x@my.text <- value; x} )

if (!isGeneric("my.number")) {
  if (is.function("my.number"))
    fun <- my.number
  else
    fun <- function(object) standardGeneric("my.number")
  setGeneric("my.number", fun)
}
setGeneric("my.number<-",
           function(x, value) standardGeneric("my.number<-"))

setMethod("my.number", "NewVertex",
          function(object) object@my.number)
setReplaceMethod("my.number", "NewVertex",
                 function(x, value) {x@my.number <- value; x} )

# Why are these 2 * 7 methods not avaliable from "dg.Vertex" ?

# setMethod("color", "NewVertex",
#           function(object) object@color)
# setReplaceMethod("color", "NewVertex",
#                  function(x, value) {x@color <- value; x} )

# setMethod("label", "NewVertex",
#           function(object) object@label)
# setReplaceMethod("label", "NewVertex",
#                  function(x, value) {x@label <- value; x} )

# setMethod("labelPosition", "NewVertex",
#           function(object) object@label.position)
# setReplaceMethod("labelPosition", "NewVertex",
#                  function(x, value) {x@label.position <- value; x} )

# setMethod("name", "NewVertex",
#           function(object) object@name)
# setReplaceMethod("name", "NewVertex",
#                  function(x, value) {x@name <- value; x} )

# setMethod("index", "NewVertex",
#           function(object) object@index)
# setReplaceMethod("index", "NewVertex",
#                  function(x, value) {x@index <- value; x} )

# setMethod("position", "NewVertex", 
#           function(object) object@position)
# setReplaceMethod("position", "NewVertex",
#                  function(x, value) {x@position <- value; x} )

# setMethod("stratum", "NewVertex",
#           function(object) object@stratum)
# setReplaceMethod("stratum", "NewVertex",
#                  function(x, value) {x@stratum <- value; x} )

# setMethod("propertyDialog", "NewVertex",
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

# V.Types <- rep("NewVertex", 6)


V.Types <- c(rep("NewVertex", 3), "Discrete", "Ordinal", "Continuous")

V.Names <- c("Sex", "Age", "Eye", "FEV", "Hair", "Shosize")
V.Labels <- paste(V.Names, 1:6, sep ="/")

From <- c(1, 2, 3, 4, 5, 6)
To   <- c(2, 3, 4, 5, 6, 1)

Z <- DynamicGraph(V.Names, V.Types, From, To, texts = c("Gryf", "gaf"),
                  labels = V.Labels, object = Object, UserMenus = Menus, 
                  updateEdgeLabels = FALSE,
                  edgeColor = "green", vertexColor = "blue", 
                  debug.strata = debug.strata, debug.edges = debug.edges, 
                  debug.position = debug.position, debug.update = debug.update,
                  vertexClasses = myVertexClasses)
