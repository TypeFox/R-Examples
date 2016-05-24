#' Class \code{"SpatialNetwork"}
#' 
#' Allow to store spatial networks, especially for rendering them
#'
#' @rdname SpatialNetwork
#' @export
#' @keywords classes spatial network sp
#' @family spnet-class
#' @slot .Data object of class \code{"list"}
#' @slot map object of class \code{"SpatialPolygons"}
#' @slot networks object of class \code{"list"}
#' @slot plot.title object of class \code{"list"}
#' @slot plot.label object of class \code{"list"}
#' @slot plot.color object of class \code{"list"}
#' @slot plot.blackwhite object of class \code{"list"}
#' @slot plot.symbol object of class \code{"list"}
#' @slot plot.arrow object of class \code{"list"}
#' @slot plot.barplot object of class \code{"list"}
#' @slot plot.legend object of class \code{"list"}
#' @slot plot.layout object of class \code{"list"}
#' @slot plot.par object of class \code{"list"}
#' @slot infos object of class \code{"list"}
#' @slot meta object of class \code{"list"}
#' @slot warnings object of class \code{"list"}
#' @slot names object of class \code{"character"}
#' @slot row.names object of class \code{"data.frameRowLabels"}
#' @slot .S3Class object of class \code{"character"}
#'    
#' @section Objects from the Class:
#' Objects can be created with the \code{\link{spnet}} function (official class builder).
#' 
#'   
#' @examples
#' people <- c("John", "Elsa", "Brian", "Kate")
#' position <- c(2,4,6,8)
#' 
#' net1.df <- data.frame(
#'   'NODE' = people,
#'   'POSITION' = position
#' )
#' 
#' net1 <- spnet.create(
#'   x = net1.df
#' )
#' net1
#' 
#' net2 <- spnet.create(
#'   x = people
#' )
#' net2
setClass(
  Class = "SpatialNetwork",
  contains = 'data.frame',
  slots = c(
    #     'edges' = 'data.frame',
    'map' = 'SpatialPolygons',
    'networks' = 'list',
    'plot.title' = 'list',
    'plot.label' = 'list',
    'plot.color' = 'list',
    'plot.blackwhite' = 'list',
    'plot.symbol' = 'list',
    'plot.arrow' = 'list',
    'plot.barplot' = 'list',
    'plot.legend' = 'list',
    'plot.layout' = 'list',
    'plot.par' = 'list',
    'infos' = 'list', # for the user
    'meta' = 'list', # for the dev
    'warnings' = 'list' # for the dev
  ),
  prototype = prototype(
    meta = list(
      date.created = Sys.time(),
      plot.color.default = list(
        background = "transparent",
        region = "#ffffff",
        node = "#f0f0f0",
        border = "#000000"
      ),
      plot.symbol.default = list(
        color = 'grey10',
        cex = 4,
        space = 0.5,
        shift.x = 0,
        shift.y = 0
      ),
      plot.arrow.default = list(
        max.networks = 5,
        color = c('dodgerblue2', 'brown3', 'darkorange3', 'olivedrab', 'hotpink3'),
        shift.x = c(0,0.2,-0.2,0,0),
        shift.y = c(0,0,0,-0.2,0.2),
        opacity = 0.9,
        thickness = 2.00,
        length.rate = 1,
        shorten = 0.3,
        line.type = 1:5,
        head.length = 0.20,
        head.type = 'curved'
      )
    )
  ),
  #   contains=character(),
  #   sealed = FALSE,
  validity = function(object) {
    flag = TRUE
    
    # .Data: NODE
    if(flag && (!'NODE' %in% names(object))){
      stop("You have to provide a 'NODE' column.")
    }
    # .Data: If POSITION exists
    if(flag && ('POSITION' %in% names(object))){
      if (length(object@map) > 0) { # there is a map
        exist.in.map <- object[,'POSITION'] %in% row.names(coordinates(object@map))
        if(!all(exist.in.map)) {
          stop(paste(
            "Some elements of the 'POSITION' column are not referenced on the map:",
            paste(object[which(!exist.in.map),'POSITION'], collapse = ', ')
          )
          )
        }
      }
    }
    
    # layout
    
    # colors
    #     color <- object@plot.color
    #     if(flag && (length(color) > 0)){
    #       if(!all(names(color) %in% c('variable', 'legend'))) stop("Elements in 'plot.color' have to be named by one of the following names: 'variable', 'legend'.")
    #       if(!'variable' %in% names(color)) stop("The 'plot.color' list should contain a 'variable' element")
    #       if(!'legend' %in% names(color)) stop("The 'plot.color' list should contain a 'legend' element")
    #       if(!color$variable %in% names(object)) stop("The 'variable' element of 'plot.color' doesn't exist in data.")
    #       exist.in.data <- names(color$legend) %in% unique(object[,color$variable])
    #       if(!all(exist.in.data)) {
    #         stop(paste(
    #           "Some values in the 'legend' referenced in 'plot.color' doesn't exist in the variable ",
    #           color$variable,
    #           ': ',
    #           paste(names(color$legend)[which(!exist.in.data)], collapse = ', '),
    #           sep = ''
    #           )
    #         )
    #       }
    #     }
    
    # symbol
    #     symbol <- object@plot.symbol
    #     if(flag && (length(symbol) > 0)){
    #       if(!all(names(symbol) %in% c('variable', 'legend', 'color', 'cex', 'space', 'shift.x', 'shift.y'))) stop("Elements in 'plot.symbol' have to be named by one of the following names: 'variable', 'legend', 'color', 'cex', 'space', 'shift.x', 'shift.y'.")
    #       if(!'variable' %in% names(symbol)) stop("The 'plot.symbol' list should contain a 'variable' element")
    #       if(!'legend' %in% names(symbol)) stop("The 'plot.symbol' list should contain a 'legend' element")
    #       if(!symbol$variable %in% names(object)) stop("The 'variable' element of 'plot.symbol' doesn't exist in data.")
    #       values.of.symbol.column <- unique(.extract.multiple.strings(object[,symbol$variable]))
    #       exist.in.data <- names(symbol$legend) %in% values.of.symbol.column
    #       if(!all(exist.in.data)) {
    #         stop(paste(
    #           "Some values in the 'legend' referenced in 'plot.symbol' doesn't exist in the variable ",
    #           symbol$variable,
    #           ': ',
    #           paste(names(symbol$legend)[which(!exist.in.data)], collapse = ', '),
    #           sep = ''
    #           )
    #         )
    #       }
    #       exist.in.symbol <- symbol$legend %in% names(spnet::.graph.symbol.list)
    #       if(!all(exist.in.symbol)) stop(
    #         paste(
    #           "Some symbol names you provided doesn't exist:",
    #           paste(symbol$legend[which(!exist.in.symbol)])
    #         ))
    #     }
    
    # barplots
    # v??rif variable num??rique, 
    
    # networks
    nets <- object@networks
    if(flag && (length(nets) > 0)){
      
      if(!(length(nets) <= object@meta$plot.arrow.default$max.networks)) {
        stop(paste("The number of networks is limited to", object@meta$plot.arrow.default$max.networks))
      }
      
      for (net in nets) {
        #         if(!'data' %in% names(net)) {
        #           stop("One of the network doesn't contain a 'data' element")
        #         }
        net.data <- net$data
        if(is.matrix(net.data)) {
          if(flag && (nrow(net.data)>0 || ncol(net.data)>0)) {
            if(nrow(net.data) != ncol(net.data)) {
              message("When network data are provided by a 'matrix', it is expected to be squared.")
              message("Here ncol=", ncol(net.data), "and nrow=", nrow(net.data))
              stop("Invalid 'networks' matrix dimensions")
            }
          }
        }
        if(flag && 'opacity' %in% net) {
          if(net$opacity < 0 || net$opacity > 1)
            stop("In each network the opacity has to be in [0;1]")
        }
      }
    }
    
    return(flag)
    
    # COLOR
    ## legend in list
    ## variable in list
    ##
  }
)


#' Extract or replace parts of a SpatialNetwork object
#'
#' @name [
#' @aliases [,SpatialNetwork-method
#' @docType methods
#' @rdname extract-methods
NULL
# setMethod(
#   "[", signature(x = "SpatialNetwork", i = "ANY", j="ANY"),
#   function (x, i, j, ..., drop) {
#     x[i,j]
#   }
# )

#' set parts of SpatialNetwork
#'
#' @name [<-
#' @aliases [<-,SpatialNetwork-method
#' @docType methods
#' @rdname extract-methods
NULL
# setReplaceMethod(
#   "[", signature(x = "SpatialNetwork", i = "ANY", j="ANY"),
#   function (x, i, j, ..., drop, value) {
#     x[i,j] <- value
#   }
# )














#' Get the map to a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the map object. Currently only \code{SpatialPolygons} from the \code{sp} package are supported.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get the map.
#' @param value the map.
#' @export
setGeneric("graph.map", function(object){ standardGeneric("graph.map") })

#' @describeIn graph.map method for \code{SpatialPolygons} objects.
setMethod(
  f = "graph.map",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "map"))
  }
)

#' Set the map to a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the map object. Currently only \code{SpatialPolygons} from the \code{sp} package are supported.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set the map.
#' @param value the map.
#' @export
setGeneric("graph.map<-", function(object, value){ standardGeneric("graph.map<-") })

#' @describeIn graph.map method for \code{SpatialPolygons} objects.
setMethod(
  f = "graph.map<-" ,
  signature = c("SpatialNetwork", 'SpatialPolygons'),
  definition = function(object, value){
    object@map <- value
    if (is.null(graph.color.background(object)))
      graph.color.background(object) <- object@meta$plot.color$background
    if (is.null(graph.color.region(object)))
      graph.color.region(object) <- object@meta$plot.color$region
    if (is.null(graph.color.node(object)))
      graph.color.node(object) <- object@meta$plot.color$node
    if (is.null(graph.color.border(object)))
      graph.color.border(object) <- object@meta$plot.color$border
    validObject(object)
    return(object)
  }
)











#' Get the list of all networks parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract networks parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.networks.list", function(object){ standardGeneric("graph.networks.list") })

#' @describeIn graph.networks.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.networks.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "networks"))
  }
)

#' Set the list of all networks parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace networks parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.networks.list<-", function(object, value){ standardGeneric("graph.networks.list<-") })

#' @describeIn graph.networks.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.networks.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@networks <- value
    validObject(object)
    return(object)
  }
)




#' Get the list of all parameters of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract all parameters of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param network.name character; the name of the network.
#' @param value a list of parameters.
#' @export
setGeneric("graph.network.list", function(object, network.name){ standardGeneric("graph.network.list") })

#' @describeIn graph.network.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.list",
  signature = c("SpatialNetwork", "character"),
  definition = function (object, network.name) { 
    return(slot(object, "networks")[[network.name]])
  }
)

#' Set the list of all parameters of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace all parameters of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param network.name character; the name of the network.
#' @param value a list of parameters.
#' @export
setGeneric("graph.network.list<-", function(object, network.name, value){ standardGeneric("graph.network.list<-") })

#' @describeIn graph.network.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.list<-" ,
  signature = c("SpatialNetwork", "character", 'list'),
  definition = function(object, network.name, value){
    object@networks[[network.name]] <- value
    validObject(object)
    return(object)
  }
)












#' Test if a network exist
#' 
#' This function tests if the network name given in parameter match the name of a network defined within a \code{SpatialNetwork} object.
#' @param object a \code{SpatialNetwork} object.
#' @param network.name a character; the name of the network.
#' @export
graph.network.exists <- function(object, network.name) {
  stopifnot(inherits(object, 'SpatialNetwork'))
  return(is.element(network.name, names(graph.networks.list(object))))
}






#' Add a network
#' 
#' This function defines a new network item in a \code{SpatialNetwork} object.
#' @param object a \code{SpatialNetwork} object.
#' @param value a character; the name of the network.
#' @export
setGeneric("graph.networks.add<-", function(object, value){ standardGeneric("graph.networks.add<-") })

#' @rdname graph.networks.add-set
setMethod(
  f = "graph.networks.add<-",
  signature = c("SpatialNetwork", "character"),
  definition = function(object, value){
    network.name <- value
    if(make.names(network.name) != network.name) {
      stop("The name is not valid. Please check it with 'make.names()'.")
    }
    if(graph.network.exists(object, network.name)) {
      stop("This network name is already defined.")
    }
    graph.networks.list(object) <- eval(parse(text = paste0("c(graph.networks.list(object), list(", network.name, " = list()))")))
    return(object)
  }
)


#' Remove a network
#' 
#' This function remove a network item in a \code{SpatialNetwork} object.
#' @param object a \code{SpatialNetwork} object.
#' @param value a character; the name of the network.
#' @export
setGeneric("graph.networks.remove<-", function(object, value){ standardGeneric("graph.networks.remove<-") })

#' @rdname graph.networks.remove-set
setMethod(
  f = "graph.networks.remove<-",
  signature = c("SpatialNetwork", "character"),
  definition = function(object, value){
    network.name <- value
    if(make.names(network.name) != network.name) {
      stop("The name is not valid. Please check it with 'make.names()'.")
    }
    if(!graph.network.exists(object, network.name)) {
      stop("This network name doesn't exist.")
    }
    if(length(graph.networks.list(object)) == 1) {
      graph.networks.list(object) <- list()
    } else {
      graph.networks.list(object) <- graph.networks.list(object)[-which(names(graph.networks.list(object)) == network.name)]
    }
    return(object)
  }
)







#' Get the data of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the data of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the network data. Currently only support a \code{matrix} object.
#' @export
setGeneric("graph.network.data", function(object, network.name){ standardGeneric("graph.network.data") })

#' @describeIn graph.network.data method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.data",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'. Please use the 'graph.network.add' function to define a network before trying to add data.")
    }
    return(object@networks[[network.name]]$data)
  }
)


#' Set the data of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the data of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the network data. Currently only support a \code{matrix} object.
#' @export
setGeneric("graph.network.data<-", function(object, network.name, value){ standardGeneric("graph.network.data<-") })

#' @describeIn graph.network.data method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.data<-" ,
  signature = c("SpatialNetwork", "character", "matrix"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'. Please use the 'graph.networks.add' function to define a network before trying to add data.")
    }
    object@networks[[network.name]]$data <- value
    validObject(object)
    return(object)
  }
)













#' Get the label of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the label of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the network label.
#' @export
setGeneric("graph.network.label", function(object, network.name){ standardGeneric("graph.network.label") })

#' @describeIn graph.network.label method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.label",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$label)
  }
)


#' Set the label of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the label of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the network label.
#' @export
setGeneric("graph.network.label<-", function(object, network.name, value){ standardGeneric("graph.network.label<-") })

#' @describeIn graph.network.label method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.label<-" ,
  signature = c("SpatialNetwork", "character", "character"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$label <- value
    validObject(object)
    return(object)
  }
)












#' Get the arrow color of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow color of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow color.
#' @export
setGeneric("graph.network.arrow.color", function(object, network.name){ standardGeneric("graph.network.arrow.color") })

#' @describeIn graph.network.arrow.color method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.color",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$color)
  }
)


#' Set the arrow color of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow color of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow color.
#' @export
setGeneric("graph.network.arrow.color<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.color<-") })

#' @describeIn graph.network.arrow.color method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.color<-" ,
  signature = c("SpatialNetwork", "character", "character"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$color <- value
    validObject(object)
    return(object)
  }
)















#' Get the arrow line type of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow line type of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value a numeric; the arrow line type.
#' @export
setGeneric("graph.network.arrow.line.type", function(object, network.name){ standardGeneric("graph.network.arrow.line.type") })

#' @describeIn graph.network.arrow.line.type method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.line.type",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$line.type)
  }
)


#' Set the arrow line type of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow line type of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value a numeric; the arrow line type.
#' @export
setGeneric("graph.network.arrow.line.type<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.line.type<-") })

#' @describeIn graph.network.arrow.line.type method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.line.type<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$line.type <- value
    validObject(object)
    return(object)
  }
)













#' Get the arrow opacity of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow opacity of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow opacity.
#' @export
setGeneric("graph.network.arrow.opacity", function(object, network.name){ standardGeneric("graph.network.arrow.opacity") })

#' @describeIn graph.network.arrow.opacity method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.opacity",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$opacity)
  }
)


#' Set the arrow opacity of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow opacity of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow opacity.
#' @export
setGeneric("graph.network.arrow.opacity<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.opacity<-") })

#' @describeIn graph.network.arrow.opacity method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.opacity<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no arrow called '", network.name, "'.")
    }
    object@networks[[network.name]]$opacity <- value
    validObject(object)
    return(object)
  }
)














#' Get the arrow thickness of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow thickness of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow thickness.
#' @export
setGeneric("graph.network.arrow.thickness", function(object, network.name){ standardGeneric("graph.network.arrow.thickness") })

#' @describeIn graph.network.arrow.thickness method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.thickness",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$thickness)
  }
)


#' Set the arrow thickness of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow thickness of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow thickness.
#' @export
setGeneric("graph.network.arrow.thickness<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.thickness<-") })

#' @describeIn graph.network.arrow.thickness method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.thickness<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$thickness <- value
    validObject(object)
    return(object)
  }
)













#' Get the arrow shortening of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow shortening of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow shortening.
#' @export
setGeneric("graph.network.arrow.shorten", function(object, network.name){ standardGeneric("graph.network.arrow.shorten") })

#' @describeIn graph.network.arrow.shorten method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.shorten",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$shorten)
  }
)


#' Set the arrow shortening of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow shortening of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow shortening.
#' @export
setGeneric("graph.network.arrow.shorten<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.shorten<-") })

#' @describeIn graph.network.arrow.shorten method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.shorten<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'. Please use the 'graph.networks.add' function to define a network before trying to add shorten.")
    }
    object@networks[[network.name]]$shorten <- value
    validObject(object)
    return(object)
  }
)

















#' Get the arrow head type of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow head type of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value type of arrowhead to draw, one of "simple", "curved", "triangle", "circle", "ellipse" or "T". See \code{\link[shape]{Arrows}} for details.
#' @export
setGeneric("graph.network.arrow.head.type", function(object, network.name){ standardGeneric("graph.network.arrow.head.type") })

#' @describeIn graph.network.arrow.head.type method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.head.type",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$head.type)
  }
)


#' Set the arrow head type of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow head type of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value type of arrowhead to draw, one of "simple", "curved", "triangle", "circle", "ellipse" or "T". See \code{\link[shape]{Arrows}} for details.
#' @export
setGeneric("graph.network.arrow.head.type<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.head.type<-") })

#' @describeIn graph.network.arrow.head.type method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.head.type<-" ,
  signature = c("SpatialNetwork", "character", "character"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$head.type <- value
    validObject(object)
    return(object)
  }
)












#' Get the arrow head length of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow head length of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow head length.
#' @export
setGeneric("graph.network.arrow.head.lth", function(object, network.name){ standardGeneric("graph.network.arrow.head.lth") })

#' @describeIn graph.network.arrow.head.lth method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.head.lth",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$head.length)
  }
)


#' Set the arrow head length of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow head length of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow head length.
#' @export
setGeneric("graph.network.arrow.head.lth<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.head.lth<-") })

#' @describeIn graph.network.arrow.head.lth method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.head.lth<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$head.length <- value
    validObject(object)
    return(object)
  }
)












#' Get the arrow shift on the x axis of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow shift on the x axis of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow shift on the x axis.
#' @export
setGeneric("graph.network.arrow.shift.x", function(object, network.name){ standardGeneric("graph.network.arrow.shift.x") })

#' @describeIn graph.network.arrow.shift.x method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.shift.x",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$shift.x)
  }
)


#' Set the arrow shift on the x axis of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow shift on the x axis of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow shift on the x axis.
#' @export
setGeneric("graph.network.arrow.shift.x<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.shift.x<-") })

#' @describeIn graph.network.arrow.shift.x method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.shift.x<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$shift.x <- value
    validObject(object)
    return(object)
  }
)










#' Get the arrow shift on the y axis of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the arrow shift on the y axis of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow shift on the y axis.
#' @export
setGeneric("graph.network.arrow.shift.y", function(object, network.name){ standardGeneric("graph.network.arrow.shift.y") })

#' @describeIn graph.network.arrow.shift.y method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.shift.y",
  signature = c("SpatialNetwork", "character"), 
  definition = function (object, network.name) {
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    return(object@networks[[network.name]]$shift.y)
  }
)


#' Set the arrow shift on the y axis of a given network of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the arrow shift on the y axis of a given network of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param network.name character; the name of the network.
#' @param value the arrow shift on the y axis.
#' @export
setGeneric("graph.network.arrow.shift.y<-", function(object, network.name, value){ standardGeneric("graph.network.arrow.shift.y<-") })

#' @describeIn graph.network.arrow.shift.y method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.network.arrow.shift.y<-" ,
  signature = c("SpatialNetwork", "character", "numeric"),
  definition = function(object, network.name, value){
    if(!graph.network.exists(object, network.name)) {
      stop("There is no network called '", network.name, "'.")
    }
    object@networks[[network.name]]$shift.y <- value
    validObject(object)
    return(object)
  }
)
















#' Get the list of all title parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract title parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @export
setGeneric("graph.title.list", function(object){ standardGeneric("graph.title.list") })

#' @describeIn graph.title.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.title.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.title"))
  }
)

#' Set the list of all title parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace title parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.title.list<-", function(object, value){ standardGeneric("graph.title.list<-") })

#' @describeIn graph.networks.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.title.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.title <- value
    validObject(object)
    return(object)
  }
)








#' Get the main title of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the main title of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new title.
#' @export
setGeneric("graph.title.main", function(object){ standardGeneric("graph.title.main") })

#' @describeIn graph.title.main method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.title.main",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.title")$main)
  }
)


#' Set the main title  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the main title of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new title.
#' @export
setGeneric("graph.title.main<-", function(object, value){ standardGeneric("graph.title.main<-") })

#' @describeIn graph.title.main method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.title.main<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.title$main <- value
    validObject(object)
    return(object)
  }
)













#' Get the sub title of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the sub title of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new title.
#' @export
setGeneric("graph.title.sub", function(object){ standardGeneric("graph.title.sub") })

#' @describeIn graph.title.sub method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.title.sub",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.title")$sub)
  }
)


#' Set the sub title  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the sub title of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new title.
#' @export
setGeneric("graph.title.sub<-", function(object, value){ standardGeneric("graph.title.sub<-") })

#' @describeIn graph.title.sub method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.title.sub<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.title$sub <- value
    validObject(object)
    return(object)
  }
)












#' Get the list of all label parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract label parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.label.list", function(object){ standardGeneric("graph.label.list") })

#' @describeIn graph.label.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.label"))
  }
)

#' Set the list of all label parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace label parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.label.list<-", function(object, value){ standardGeneric("graph.label.list<-") })

#' @describeIn graph.label.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.label <- value
    validObject(object)
    return(object)
  }
)















#' Get the label variable of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the label variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new label, for example "#000000".
#' @export
setGeneric("graph.label.variable", function(object){ standardGeneric("graph.label.variable") })

#' @describeIn graph.label.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.variable",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.label")$variable)
  }
)

#' Set the label variable  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the label variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new label, for example "#000000".
#' @export
setGeneric("graph.label.variable<-", function(object, value){ standardGeneric("graph.label.variable<-") })

#' @describeIn graph.label.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.variable<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.label$variable <- value
    validObject(object)
    return(object)
  }
)













#' Get the label cex of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the label cex of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value numeric; the cex parameter.
#' @export
setGeneric("graph.label.cex", function(object){ standardGeneric("graph.label.cex") })

#' @describeIn graph.label.cex method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.cex",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.label")$cex)
  }
)

#' Set the label cex  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the label cex of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value numeric; the cex parameter.
#' @export
setGeneric("graph.label.cex<-", function(object, value){ standardGeneric("graph.label.cex<-") })

#' @describeIn graph.label.cex method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.cex<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.label$cex <- value
    validObject(object)
    return(object)
  }
)















#' Get the label color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the label color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new label, for example "#000000".
#' @export
setGeneric("graph.label.color", function(object){ standardGeneric("graph.label.color") })

#' @describeIn graph.label.color method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.color",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.label")$col)
  }
)

#' Set the label color  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the label color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new label, for example "#000000".
#' @export
setGeneric("graph.label.color<-", function(object, value){ standardGeneric("graph.label.color<-") })

#' @describeIn graph.label.color method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.label.color<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.label$col <- value
    validObject(object)
    return(object)
  }
)


















#' Get the list of all color parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract color parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.color.list", function(object){ standardGeneric("graph.color.list") })

#' @describeIn graph.color.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color"))
  }
)

#' Set the list of all color parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace color parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.color.list<-", function(object, value){ standardGeneric("graph.color.list<-") })

#' @describeIn graph.color.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.color <- value
    validObject(object)
    return(object)
  }
)














#' Get the color variable of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the color variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new color, for example "#000000".
#' @export
setGeneric("graph.color.variable", function(object){ standardGeneric("graph.color.variable") })

#' @describeIn graph.color.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.variable",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color")$variable)
  }
)

#' Set the color variable  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the color variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new color, for example "#000000".
#' @export
setGeneric("graph.color.variable<-", function(object, value){ standardGeneric("graph.color.variable<-") })

#' @describeIn graph.color.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.variable<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.color$variable <- value
    validObject(object)
    return(object)
  }
)













#' Get the color legend of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the color legend of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the color legend.
#' @export
setGeneric("graph.color.legend", function(object){ standardGeneric("graph.color.legend") })

#' @describeIn graph.color.legend method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.legend",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color")$legend)
  }
)

#' Set the color legend  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the color legend of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the color legend.
#' @export
setGeneric("graph.color.legend<-", function(object, value){ standardGeneric("graph.color.legend<-") })

#' @describeIn graph.color.legend method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.legend<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.color$legend <- value
    validObject(object)
    return(object)
  }
)












#' Get the background color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the background color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.background", function(object){ standardGeneric("graph.color.background") })

#' @describeIn graph.color.background method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.background",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color")$background)
  }
)

#' Set the background color  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the background color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.background<-", function(object, value){ standardGeneric("graph.color.background<-") })

#' @describeIn graph.color.background method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.background<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.color$background <- value
    validObject(object)
    return(object)
  }
)











#' Get the default color of a region of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the default color of a region of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.region", function(object){ standardGeneric("graph.color.region") })

#' @describeIn graph.color.region method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.region",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color")$region)
  }
)

#' Set the default color of a region  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the default color of a region of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.region<-", function(object, value){ standardGeneric("graph.color.region<-") })

#' @describeIn graph.color.region method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.region<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.color$region <- value
    validObject(object)
    return(object)
  }
)









#' Get the default color of a node of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the default color of a node of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.node", function(object){ standardGeneric("graph.color.node") })

#' @describeIn graph.color.node method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.node",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color")$node)
  }
)

#' Set the default color of a node  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the default color of a node of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.node<-", function(object, value){ standardGeneric("graph.color.node<-") })

#' @describeIn graph.color.node method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.node<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.color$node <- value
    validObject(object)
    return(object)
  }
)










#' Get the border color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the border color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.border", function(object){ standardGeneric("graph.color.border") })

#' @describeIn graph.color.border method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.border",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.color")$border)
  }
)

#' Set the border color  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the border color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{character}, the color.
#' @export
setGeneric("graph.color.border<-", function(object, value){ standardGeneric("graph.color.border<-") })

#' @describeIn graph.color.border method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.color.border<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.color$border <- value
    validObject(object)
    return(object)
  }
)















#' Get the list of all black and white mode parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract black and white mode parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.blackwhite.list", function(object){ standardGeneric("graph.blackwhite.list") })

#' @describeIn graph.blackwhite.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.blackwhite"))
  }
)

#' Set the list of all black and white mode parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace black and white mode parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.blackwhite.list<-", function(object, value){ standardGeneric("graph.blackwhite.list<-") })

#' @describeIn graph.blackwhite.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.blackwhite <- value
    validObject(object)
    return(object)
  }
)










#' Get the black and white mode status of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the black and white mode status of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{logical}, the black and white mode status.
#' @export
setGeneric("graph.blackwhite.enable", function(object){ standardGeneric("graph.blackwhite.enable") })

#' @describeIn graph.blackwhite.enable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.enable",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.blackwhite")$enable)
  }
)

#' Set the black and white mode status  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the black and white mode status of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{logical}, the black and white mode status.
#' @export
setGeneric("graph.blackwhite.enable<-", function(object, value){ standardGeneric("graph.blackwhite.enable<-") })

#' @describeIn graph.blackwhite.enable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.enable<-" ,
  signature = c("SpatialNetwork", 'logical'),
  definition = function(object, value){
    object@plot.blackwhite$enable <- value
    validObject(object)
    return(object)
  }
)











#' Get the black and white mode minimal gray value of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the black and white mode minimal gray value (from 0 to 1) of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{logical}, the black and white mode minimal gray value.
#' @export
setGeneric("graph.blackwhite.min", function(object){ standardGeneric("graph.blackwhite.min") })

#' @describeIn graph.blackwhite.min method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.min",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.blackwhite")$min)
  }
)

#' Set the black and white mode minimal gray value of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the black and white mode minimal gray value (from 0 to 1) of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{numeric}, the black and white mode minimal gray value.
#' @export
setGeneric("graph.blackwhite.min<-", function(object, value){ standardGeneric("graph.blackwhite.min<-") })

#' @describeIn graph.blackwhite.min method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.min<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.blackwhite$min <- value
    validObject(object)
    return(object)
  }
)












#' Get the black and white mode maximal gray value of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the black and white mode maximal gray value (from 0 to 1) of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{logical}, the black and white mode maximal gray value.
#' @export
setGeneric("graph.blackwhite.max", function(object){ standardGeneric("graph.blackwhite.max") })

#' @describeIn graph.blackwhite.max method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.max",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.blackwhite")$max)
  }
)

#' Set the black and white mode maximal gray value of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the black and white mode maximal gray value (from 0 to 1) of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a \code{numeric}, the black and white mode maximal gray value.
#' @export
setGeneric("graph.blackwhite.max<-", function(object, value){ standardGeneric("graph.blackwhite.max<-") })

#' @describeIn graph.blackwhite.max method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.blackwhite.max<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.blackwhite$max <- value
    validObject(object)
    return(object)
  }
)
















#' Get the list of all symbol parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract symbol parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.symbol.list", function(object){ standardGeneric("graph.symbol.list") })

#' @describeIn graph.symbol.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol"))
  }
)

#' Set the list of all symbol parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace symbol parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.symbol.list<-", function(object, value){ standardGeneric("graph.symbol.list<-") })

#' @describeIn graph.symbol.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.symbol <- value
    validObject(object)
    return(object)
  }
)















#' Get the symbol variable of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the symbol variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the symbol variable.
#' @export
setGeneric("graph.symbol.variable", function(object){ standardGeneric("graph.symbol.variable") })

#' @describeIn graph.symbol.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.variable",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol")$variable)
  }
)

#' Set the symbol variable  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the symbol variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the symbol variable.
#' @export
setGeneric("graph.symbol.variable<-", function(object, value){ standardGeneric("graph.symbol.variable<-") })

#' @describeIn graph.symbol.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.variable<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.symbol$variable <- value
    validObject(object)
    return(object)
  }
)

















#' Get the symbol legend of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the symbol legend of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new legend.
#' @export
setGeneric("graph.symbol.legend", function(object){ standardGeneric("graph.symbol.legend") })

#' @describeIn graph.symbol.legend method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.legend",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol")$legend)
  }
)

#' Set the symbol legend  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the symbol legend of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new legend.
#' @export
setGeneric("graph.symbol.legend<-", function(object, value){ standardGeneric("graph.symbol.legend<-") })

#' @describeIn graph.symbol.legend method for \code{SpatialNetwork} objects.

setMethod(
  f = "graph.symbol.legend<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.symbol$legend <- value
    validObject(object)
    return(object)
  }
)












#' Get the symbol cex parameter of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the symbol cex parameter of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new cex parameter.
#' @export
setGeneric("graph.symbol.cex", function(object){ standardGeneric("graph.symbol.cex") })

#' @describeIn graph.symbol.cex method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.cex",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol")$cex)
  }
)

#' Set the symbol cex parameter  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the symbol cex parameter of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new cex parameter.
#' @export
setGeneric("graph.symbol.cex<-", function(object, value){ standardGeneric("graph.symbol.cex<-") })

#' @describeIn graph.symbol.cex method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.cex<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.symbol$cex <- value
    validObject(object)
    return(object)
  }
)

















#' Get the symbol color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the symbol color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the color.
#' @export
setGeneric("graph.symbol.color", function(object){ standardGeneric("graph.symbol.color") })

#' @describeIn graph.symbol.color method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.color",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol")$color)
  }
)

#' Set the symbol color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the symbol color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the color.
#' @export
setGeneric("graph.symbol.color<-", function(object, value){ standardGeneric("graph.symbol.color<-") })

#' @describeIn graph.symbol.color method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.color<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.symbol$color <- value
    validObject(object)
    return(object)
  }
)














#' Get the symbol shift on the x axis of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the value of symbol shift on the x axis of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric; the value of the shift.s
#' @export
setGeneric("graph.symbol.shift.x", function(object){ standardGeneric("graph.symbol.shift.x") })

#' @describeIn graph.symbol.shift.x method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.shift.x",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol")$shift.x)
  }
)

#' Set the symbol shift on the x axis of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the value of symbol shift on the x axis of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric; the value of the shift.
#' @export
setGeneric("graph.symbol.shift.x<-", function(object, value){ standardGeneric("graph.symbol.shift.x<-") })

#' @describeIn graph.symbol.shift.x method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.shift.x<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.symbol$shift.x <- value
    validObject(object)
    return(object)
  }
)












#' Get the symbol shift on the y axis of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the value of the symbol shift on the y of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric; the value of the shift.
#' @export
setGeneric("graph.symbol.shift.y", function(object){ standardGeneric("graph.symbol.shift.y") })

#' @describeIn graph.symbol.shift.y method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.shift.y",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.symbol")$shift.y)
  }
)

#' Set the symbol shift on the y axis of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the value of the symbol shift on the y axis of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric; the value of the shift.
#' @export
setGeneric("graph.symbol.shift.y<-", function(object, value){ standardGeneric("graph.symbol.shift.y<-") })

#' @describeIn graph.symbol.shift.y method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.symbol.shift.y<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.symbol$shift.y <- value
    validObject(object)
    return(object)
  }
)












#' Get the list of all barplot parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract barplot parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.barplot.list", function(object){ standardGeneric("graph.barplot.list") })

#' @describeIn graph.barplot.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot"))
  }
)

#' Set the list of all barplot parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace barplot parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.barplot.list<-", function(object, value){ standardGeneric("graph.barplot.list<-") })

#' @describeIn graph.barplot.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.barplot <- value
    validObject(object)
    return(object)
  }
)

















#' Get the barplot variable of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the barplot variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the name of the variable to use for plotting barplots.
#' @export
setGeneric("graph.barplot.variable", function(object){ standardGeneric("graph.barplot.variable") })

#' @describeIn graph.barplot.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.variable",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot")$variable)
  }
)

#' Set the barplot variable  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the barplot variable of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the name of the variable to use for plotting barplots.
#' @export
setGeneric("graph.barplot.variable<-", function(object, value){ standardGeneric("graph.barplot.variable<-") })

#' @describeIn graph.barplot.variable method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.variable<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.barplot$variable <- value
    validObject(object)
    return(object)
  }
)
















#' Get the barplot foreground color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the barplot foreground color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the color.
#' @export
setGeneric("graph.barplot.fgcolor", function(object){ standardGeneric("graph.barplot.fgcolor") })

#' @describeIn graph.barplot.fgcolor method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.fgcolor",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot")$fgcolor)
  }
)


#' Set the barplot foreground color  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the barplot foreground color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the color.
#' @export
setGeneric("graph.barplot.fgcolor<-", function(object, value){ standardGeneric("graph.barplot.fgcolor<-") })

#' @describeIn graph.barplot.fgcolor method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.fgcolor<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.barplot$fgcolor <- value
    validObject(object)
    return(object)
  }
)


















#' Get the barplot background color of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the barplot background color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new color.
#' @export
setGeneric("graph.barplot.bgcolor", function(object){ standardGeneric("graph.barplot.bgcolor") })

#' @describeIn graph.barplot.bgcolor method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.bgcolor",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot")$bgcolor)
  }
)

#' Set the barplot background color  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the barplot background color of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value the new color.
#' @export
setGeneric("graph.barplot.bgcolor<-", function(object, value){ standardGeneric("graph.barplot.bgcolor<-") })

#' @describeIn graph.barplot.bgcolor method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.bgcolor<-" ,
  signature = c("SpatialNetwork", 'character'),
  definition = function(object, value){
    object@plot.barplot$bgcolor <- value
    validObject(object)
    return(object)
  }
)



















#' Get the barplot lower bound position of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the barplot lower bound position of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric vector of coordinates, c(x,y), specifying a shift from the center of each country.
#' @export
setGeneric("graph.barplot.bound.lower", function(object){ standardGeneric("graph.barplot.bound.lower") })

#' @describeIn graph.barplot.bound.lower method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.bound.lower",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot")$bound.lower)
  }
)

#' Set the barplot lower bound position  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the barplot lower bound position of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric vector of coordinates, c(x,y), specifying a shift from the center of each country.
#' @export
setGeneric("graph.barplot.bound.lower<-", function(object, value){ standardGeneric("graph.barplot.bound.lower<-") })

#' @describeIn graph.barplot.bound.lower method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.bound.lower<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.barplot$bound.lower <- value
    validObject(object)
    return(object)
  }
)
















#' Get the barplot upper bound position of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the barplot upper bound position of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric vector of coordinates, c(x,y), specifying a shift from the center of each country.
#' @export
setGeneric("graph.barplot.bound.upper", function(object){ standardGeneric("graph.barplot.bound.upper") })

#' @describeIn graph.barplot.bound.upper method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.bound.upper",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot")$bound.upper)
  }
)

#' Set the barplot upper bound position  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the barplot upper bound position of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric vector of coordinates, c(x,y), specifying a shift from the center of each country.
#' @export
setGeneric("graph.barplot.bound.upper<-", function(object, value){ standardGeneric("graph.barplot.bound.upper<-") })

#' @describeIn graph.barplot.bound.upper method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.bound.upper<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.barplot$bound.upper <- value
    validObject(object)
    return(object)
  }
)















#' Get the barplot width of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the barplot width of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric.
#' @export
setGeneric("graph.barplot.width", function(object){ standardGeneric("graph.barplot.width") })

#' @describeIn graph.barplot.width method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.width",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.barplot")$width)
  }
)

#' Set the barplot width  of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the barplot width of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric.
#' @export
setGeneric("graph.barplot.width<-", function(object, value){ standardGeneric("graph.barplot.width<-") })

#' @describeIn graph.barplot.width method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.barplot.width<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.barplot$width <- value
    validObject(object)
    return(object)
  }
)
















#' Get the list of all legend parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract legend parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export

setGeneric("graph.legend.list", function(object){ standardGeneric("graph.legend.list") })

#' @describeIn graph.legend.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.legend"))
  }
)

#' Set the list of all legend parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace legend parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.legend.list<-", function(object, value){ standardGeneric("graph.legend.list<-") })

#' @describeIn graph.legend.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.legend <- value
    validObject(object)
    return(object)
  }
)














#' Get the legend print (yes/no) status of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the legend print (yes/no) status of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a logical.
#' @export
setGeneric("graph.legend.print", function(object){ standardGeneric("graph.legend.print") })

#' @describeIn graph.legend.print method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.print",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.legend")$print)
  }
)

#' Set the legend print (yes/no) status of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the legend print (yes/no) status of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a logical.
#' @export
setGeneric("graph.legend.print<-", function(object, value){ standardGeneric("graph.legend.print<-") })

#' @describeIn graph.legend.print method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.print<-" ,
  signature = c("SpatialNetwork", 'logical'),
  definition = function(object, value){
    object@plot.legend$print <- value
    validObject(object)
    return(object)
  }
)








#' Get the legend cex parameter of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the legend cex parameter of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric.
#' @export
setGeneric("graph.legend.cex", function(object){ standardGeneric("graph.legend.cex") })

#' @describeIn graph.legend.cex method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.cex",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.legend")$cex)
  }
)

#' Set the legend cex parameter of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the legend cex parameter of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric.
#' @export
setGeneric("graph.legend.cex<-", function(object, value){ standardGeneric("graph.legend.cex<-") })

#' @describeIn graph.legend.cex method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.cex<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.legend$cex <- value
    validObject(object)
    return(object)
  }
)





#' Get the legend number of columns of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the legend number of columns of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric.
#' @export
setGeneric("graph.legend.ncol", function(object){ standardGeneric("graph.legend.ncol") })

#' @describeIn graph.legend.ncol method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.ncol",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.legend")$ncol)
  }
)

#' Set the legend number of columns of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the legend number of columns of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a numeric.
#' @export
setGeneric("graph.legend.ncol<-", function(object, value){ standardGeneric("graph.legend.ncol<-") })

#' @describeIn graph.legend.ncol method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.ncol<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.legend$ncol <- value
    validObject(object)
    return(object)
  }
)






#' Get the legend horizontal or vertical setting of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the legend horizontal or vertical setting of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a logical.
#' @export
setGeneric("graph.legend.horiz", function(object){ standardGeneric("graph.legend.horiz") })

#' @describeIn graph.legend.horiz method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.horiz",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.legend")$horiz)
  }
)

#' Set the legend horizontal or vertical setting of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the legend horizontal or vertical setting of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a logical.
#' @export
setGeneric("graph.legend.horiz<-", function(object, value){ standardGeneric("graph.legend.horiz<-") })

#' @describeIn graph.legend.horiz method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.horiz<-" ,
  signature = c("SpatialNetwork", 'logical'),
  definition = function(object, value){
    object@plot.legend$horiz <- value
    validObject(object)
    return(object)
  }
)





#' Get the legend line width parameter of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract the legend line width parameter of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a logical.
#' @export
setGeneric("graph.legend.line.width", function(object){ standardGeneric("graph.legend.line.width") })

#' @describeIn graph.legend.line.width method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.line.width",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.legend")$lwd)
  }
)

#' Set the legend line width parameter of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace the legend line width parameter of a \code{SpatialNetwork} object.
#' 
#' @param object a \code{SpatialNetwork} object.
#' @param value a logical.
#' @export
setGeneric("graph.legend.line.width<-", function(object, value){ standardGeneric("graph.legend.line.width<-") })

#' @describeIn graph.legend.line.width method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.legend.line.width<-" ,
  signature = c("SpatialNetwork", 'numeric'),
  definition = function(object, value){
    object@plot.legend$lwd <- value
    validObject(object)
    return(object)
  }
)







#' Get the list of all layout parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract layout parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.layout.list", function(object){ standardGeneric("graph.layout.list") })

#' @describeIn graph.layout.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.layout.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.layout"))
  }
)

#' Set the list of all layout parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace layout parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.layout.list<-", function(object, value){ standardGeneric("graph.layout.list<-") })

#' @describeIn graph.layout.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.layout.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.layout <- value
    validObject(object)
    return(object)
  }
)

















#' Get the list of all par parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to extract par parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to get parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.par.list", function(object){ standardGeneric("graph.par.list") })

#' @describeIn graph.par.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.par.list",
  signature = "SpatialNetwork", 
  definition = function (object) { 
    return(slot(object, "plot.par"))
  }
)

#' Set the list of all par parameters of a \code{SpatialNetwork} object
#' 
#' This generic method intends to set or replace par parameters of a \code{SpatialNetwork} object.
#' 
#' @param object the \code{SpatialNetwork} object for which we want to set parameters.
#' @param value a list of parameters.
#' @export
setGeneric("graph.par.list<-", function(object, value){ standardGeneric("graph.par.list<-") })

#' @describeIn graph.par.list method for \code{SpatialNetwork} objects.
setMethod(
  f = "graph.par.list<-" ,
  signature = c("SpatialNetwork", 'list'),
  definition = function(object, value){
    object@plot.par <- value
    validObject(object)
    return(object)
  }
)









#' Create a \code{SpatialNetwork} object
#' 
#' The \code{spnet.create} function is the official builder for creating \code{SpatialNetwork} objects.
#' 
#' @author Emmanuel Rousseaux
#' 
#' @param x a \code{data.frame} containing at least two columns: \code{NODE} and \code{POSITION}.
#' @param map a \code{\link[sp]{SpatialPolygons}} object.
#' @param networks a list of the networks to plot.
#' @param plot.title a list of parameters for setting the title.
#' @param plot.label a list of parameters to be passed to the \code{\link{text}} function for setting labels.
#' @param plot.color a list of parameters for setting colors.
#' @param plot.blackwhite a list of parameters for setting the black and white mode.
#' @param plot.symbol a list of parameters for setting symbols.
#' @param plot.barplot a list of parameters for setting barplots.
#' @param plot.arrow a list of parameters for setting arrows.
#' @param plot.legend a list of parameters for setting the legend.
#' @param plot.layout a list of parameters for setting the layout.
#' @param plot.par a list of graphical parameters.
#' @param infos a list of meta information about the instance of the object.
#' @param quiet = FALSE a logical, suppress all messages.
#' @export
#' @examples
#' people <- c("John", "Elsa", "Brian", "Kate")
#' position <- c(2,4,6,8)
#' 
#' net1.df <- data.frame(
#'   'NODE' = people,
#'   'POSITION' = position
#' )
#' 
#' net1 <- spnet.create(
#'   x = net1.df
#' )
#' net1
#' 
#' net2 <- spnet.create(
#'   x = people
#' )
#' net2
spnet.create <- function(
  x,
  map,
  networks,
  plot.title = list(main = "Untitled SPNET object", sub = "", cex = 2, col = "#333333"),
  plot.label = list(cex = 1, col = '#333333'),
  plot.color,
  plot.blackwhite = list(enable = FALSE, min = 0.02, max = 0.98),
  plot.symbol,
  plot.barplot = list(variable = "", bound.lower = c(-0.5,-0.5), bound.upper = c(0.5,-0.5), fgcolor = "#666666", bgcolor = "#eeeeee", width = 8),
  plot.arrow,
  plot.legend = list(print = TRUE, cex = 1, ncol = 1, horiz = FALSE, lwd = 1),
  plot.layout = list(ratios = c('title' = 1/10, 'graphic' = 7/10, 'legend' = 2/10), mat = NULL, reset = TRUE),
  plot.par = list(mar = c(1,1,1,1)), # par(mar = c(5,4,2,2))
  infos,
  quiet = FALSE
) {
  
  if(inherits(x, 'data.frame')) {
    df <- x
  } else if(length(x) > 0) {
    df <- data.frame('NODE' = x)
  } else stop("Invalid 'x' argument")
  
  out <- new(
    Class = 'SpatialNetwork',
    .Data = df,
    row.names = 1:nrow(df),
    plot.title =plot.title,
    plot.label = plot.label,
    plot.blackwhite = plot.blackwhite,
    plot.barplot = plot.barplot,
    plot.legend = plot.legend,
    plot.layout = plot.layout,
    plot.par = plot.par
  )
  
  if(!missing(networks)) {
    graph.networks.list(out) <- networks
  }
  
  return(out)
}


#' @importFrom utils head
setMethod(
  f = 'show',
  signature = 'SpatialNetwork',
  definition = function(object) {
    cat("This is a valid 'SpatialNetwork' object.\n\n")
    
    cat("- Data: (first rows) \n\n")
    print.data.frame(head(object))
    cat("\n")
    
    if(length(graph.map(object)) > 0) {
      cat("- Map:\n")
      cat("    Length:", length(graph.map(object)), "\n\n")
    }
    
    if(length(object@networks) > 0) {
      cat("- Network data:\n")
      #       for (net in nets) {
      cat("    ", "Number of network(s): ", length(object@networks), sep = "")
      #       }
      cat("\n\n")
    }
    
    color <- object@plot.color
    symbol <- object@plot.symbol
    barplot <- object@plot.barplot
    
    if('variable' %in% c(names(color), names(symbol), names(barplot))) {
      cat("- Plotting options:\n")
      if('variable' %in% names(color))
        cat("    ", "Variable used to colorize: '", color$variable, "'\n", sep = "")
      if('variable' %in% names(symbol))
        cat("    ", "Variable used to draw symbols: '", symbol$variable, "'\n", sep = "")
      if(nzchar(barplot$variable))
        cat("    ", "Variable used to draw barplots: '", barplot$variable, "'\n", sep = "")
      cat("\n")
    }
  }
)

setMethod(
  f = 'print',
  signature = 'SpatialNetwork',
  definition = function(x, ...) {
    show(x)
  }
)


#' @importFrom graphics layout
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom grDevices col2rgb
#' @importFrom grDevices rgb
setMethod(
  f = 'plot',
  signature = 'SpatialNetwork',
  definition = function(x, ...) {
    if(length(graph.map(x)) == 0)
      stop("The map is empty. Please define a valid map.")
    
    tit <- x@plot.title
    
    color <- x@plot.color
    flag.color <- ifelse(!is.null(color$variable) && !is.null(color$legend), T, F)
    
    symbol <- x@plot.symbol
    flag.symbol <- ifelse(!is.null(symbol$variable) && !is.null(symbol$legend), T, F)
    if(flag.symbol) {
      if('color' %in% names(symbol)) {
        symb.color <- symbol$color
      } else {
        symb.color <- x@meta$plot.symbol.default$color
      }
    }
    
    barplot <- x@plot.barplot
    flag.barplot <- ifelse(nzchar(barplot$variable), T, F)
    
    nets <- x@networks
    if(length(nets) > 0) {flag.arrow <- T} else {flag.arrow <- F}
    
    lay <- x@plot.layout
    
    arg.col <- numeric()
    
    # black and white mode
    plot.blackwhite <- x@plot.blackwhite
    flag.blackwhite <- plot.blackwhite$enable
    
    if(flag.blackwhite) {
      allcolors = character(0)
      
      # we build the vector of all colors that will be printed
      allcolors = c(allcolors, graph.color.background(x)) # background
      allcolors = c(allcolors, graph.color.border(x)) # border
      allcolors = c(allcolors, graph.color.node(x)) # node (default)
      allcolors = c(allcolors, graph.color.region(x)) # region
      if(flag.color)
        allcolors = c(allcolors, color$legend) # node (default)
      if(flag.symbol)
        allcolors = c(allcolors, symb.color) # symbol
      
      allcolors = sort(unique(allcolors))
      
      allcolorsbw = color2blackwhite(
        allcolors,
        contrast.min = plot.blackwhite$min,
        contrast.max = plot.blackwhite$max
      )
      
      # we replace each color to its corresponding black and white color
      graph.color.background(x) <- allcolorsbw[match(graph.color.background(x), allcolors)] # background
      graph.color.border(x) <- allcolorsbw[match(graph.color.border(x), allcolors)] # border
      graph.color.node(x) <- allcolorsbw[match(graph.color.node(x), allcolors)] # node
      graph.color.region(x) <- allcolorsbw[match(graph.color.region(x), allcolors)] # region
      if(flag.color){
        oldcolors <- color$legend
        newcolors <- allcolorsbw[match(color$legend, allcolors)]
        color$legend <- newcolors
        names(color$legend) <- names(oldcolors)
      }
      if(flag.symbol)
        symb.color <- allcolorsbw[match(symb.color, allcolors)] # symbol
    }
    
    ## PLOT
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    #     nf <- layout(
    #       mat = matrix(c(1,2),nrow=2,byrow = TRUE),
    #       widths = c(6.8),
    #       heights = c(4,1),
    #       respect = TRUE
    #     )
    
    if(!is.null(lay$mat)) {
      #       nf <- layout(
      #         mat = matrix(c(1,3,2,2),nrow=2,byrow = TRUE),
      #         widths = c(3,3),
      #         heights = c(4,1), 
      #         respect = TRUE
      #       )
    } else {
      nf <- layout(
        mat = matrix(c(1,2,3),nrow=3,byrow = TRUE),
        widths = c(1),
        heights = lay$ratios, 
        respect = TRUE
      )
    }
    #     layout.show(nf)
    par(graph.par.list(x))
    
    plot.new()
    # plot the title
    if(any(c('main', 'sub') %in% names(tit))) {
      tit.main.clean = tit[-which(names(tit) %in% c('main', 'sub'))]
    } else {
      tit.main.clean = tit
    }
    main.call = as.call(c(
      list(
        fun=text,
        x=0.5,
        y=0.65,
        labels = tit$main
      ),
      tit.main.clean
    ))
    eval(main.call)
    if(any(c('main', 'sub', 'cex') %in% names(tit))) {
      tit.sub.clean = tit[-which(names(tit) %in% c('main', 'sub', 'cex'))]
    } else {
      tit.sub.clean = tit
    }
    sub.call = as.call(c(
      list(
        fun=text,
        x=0.5,
        y=0.20,
        labels = tit$sub,
        cex = tit$cex / 1.5
      ),
      tit.sub.clean
    ))
    eval(sub.call)
    
    if(!'POSITION' %in% names(x)) { # we only plot the map
      plot(graph.map(x), ... = ...)
      plot.new()
    } else { # we plot position referenced
      
      coord <- coordinates(x@map)
      ids <- row.names(coord)
      seats <- x[, 'POSITION']
      seats.which <- match(seats, ids)
      
      col <- rep(graph.color.region(x), nrow(coord))
      col[seats.which] <- graph.color.node(x)
      
      if(flag.color) {
        names(seats.which) <- x[, color$variable] # we store color labels in the names
        names(seats.which)[!names(seats.which) %in% names(color$legend)] <- graph.color.node(x)
        for(k in names(color$legend)){
          names(seats.which)[names(seats.which) == k] <- color$legend[[k]]
        }
        col[seats.which] <- names(seats.which)
      }
      
      
      
      plot(
        x@map,
#         col = arg.col,
        col = col,
        border = graph.color.border(x),
        bg = graph.color.background(x),
        ... = ...
      )
      
      # Then we also plot the NODE names to corresponding positions
      coord <- coordinates(x@map)
      ids <- row.names(coord)
      seats <- x[, 'POSITION']
      seats.which <- match(seats, ids)
      lab <- rep("", nrow(coord))
      
      lab.opt <- graph.label.list(x)
      if(any(c('variable') %in% names(lab.opt))) {
        lab.opt.clean = lab.opt[-which(names(lab.opt) %in% c('variable'))]
        lab[seats.which] <- as.character(x[, lab.opt$variable])
      } else {
        lab.opt.clean = lab.opt
        lab[seats.which] <- as.character(x[, 'NODE'])
      }
      
      label.call = as.call(c(
        list(
          fun=text,
          x=coord,
          labels = lab
        ),
        lab.opt.clean
      ))
      eval(label.call)
      
      
      
      if(flag.symbol) {
        coord <- coordinates(x@map)
        ids <- row.names(coord)
        seats <- x[, 'POSITION']
        symb.coord <- coord[match(seats, as.numeric(ids)),] # symb.coord give the coordinates of each existing seat
        names(seats) <- x[,symbol$variable] # we store the data.frame column to match to symbols in seat names
        
        seats <- .expand.multiple.names(seats)
        
        role.match <- names(seats) %in% names(symbol$legend)
        seats <- seats[role.match] # then if no role match we remove
        
        #symb.coord <- symb.coord[seats,]
        
        if('cex' %in% names(symbol)) {
          symb.cex <- symbol$cex
        } else {
          symb.cex <- x@meta$plot.symbol.default$cex
        }
        if('space' %in% names(symbol)) {
          symb.space <- symbol$space
        } else {
          symb.space <- x@meta$plot.symbol.default$space
        }
        
        symb.coord.multiple.flag <- TRUE
        for(l in unique(seats)){
          nb <- sum(seats == l)
          val <- .position.multiple.symbols(symb.coord[match(l, rownames(symb.coord)),], n = nb, cex = symb.cex, space = symb.space)
          dimnames(val) <- list(rep(l, nb), NULL) # rownames
          if(symb.coord.multiple.flag) {
            symb.coord.multiple <- val
            symb.coord.multiple.flag <- FALSE
          }
          else {
            symb.coord.multiple <- rbind(symb.coord.multiple, val)
          }
        }      
        
        for(k in names(symbol$legend)){ # we replace by the symbol code
          role.match2 <- names(seats) == k # we select elements which match
          if(any(role.match2)){
            names(seats)[role.match2] <- symbol$legend[[k]] # and we replace
          }
        }
        #         print(seats)
        arg.pch <- names(seats)
        allsymb <- .graph.symbol.list
        arg.pch <- allsymb[match(arg.pch, names(allsymb))]
        
#         MOVED UP FOR SETTING THE B&W MODE
#         if('color' %in% names(symbol)) {
#           symb.color <- symbol$color
#         } else {
#           symb.color <- x@meta$plot.symbol.default$color
#         }
        
        if('shift.x' %in% names(symbol)) {
          symb.shift.x <- symbol$shift.x
        } else {
          symb.shift.x <- x@meta$plot.symbol.default$shift.x
        }
        if('shift.y' %in% names(symbol)) {
          symb.shift.y <- symbol$shift.y
        } else {
          symb.shift.y <- x@meta$plot.symbol.default$shift.y
        }
        
        symb.coord.multiple[,1] <- symb.coord.multiple[,1] + symb.shift.x
        symb.coord.multiple[,2] <- symb.coord.multiple[,2] + symb.shift.y
        
        points(
          symb.coord.multiple,
          pch = arg.pch,
          cex = symb.cex,
          col = symb.color
        )
      }
      
      if(flag.barplot) {
        coord <- coordinates(x@map)
        ids <- row.names(coord)
        seats <- x[, 'POSITION']
        seats.which <- match(seats, ids)
        
        values <- x[, graph.barplot.variable(x)]
        
        for(i in 1:length(values)) {
          value <- values[i]
          if(!is.na(value)) {
            lines.barplot(
              value = value,
              bound.lower = coord[seats.which[i],] + graph.barplot.bound.lower(x),
              bound.upper = coord[seats.which[i],] + graph.barplot.bound.upper(x),
              bgcolor = graph.barplot.bgcolor(x),
              fgcolor = graph.barplot.fgcolor(x),
              lwd = graph.barplot.width(x)
            )
          }
        }
        
      }
      
      if(flag.arrow) {
        arrow.col.list <- c()
        arrow.shift.x.list <- c()
        arrow.shift.y.list <- c()
        arrow.label.list <- names(nets)
        
        coord <- coordinates(x@map)
        ids <- row.names(coord)
        seats <- x[, 'POSITION']
        seats.which <- match(seats, ids)
        names(seats.which) <- x[, 'NODE']
        
        default.color = x@meta$plot.arrow.default$color
        default.shift.x = x@meta$plot.arrow.default$shift.x
        default.shift.y = x@meta$plot.arrow.default$shift.y
        default.opacity = x@meta$plot.arrow.default$opacity
        default.thickness = x@meta$plot.arrow.default$thickness
        default.length.rate = x@meta$plot.arrow.default$length.rate
        default.shorten = x@meta$plot.arrow.default$shorten
        default.head.length = x@meta$plot.arrow.default$head.length
        default.head.type = x@meta$plot.arrow.default$head.type
        
        for (k in 1:length(nets)) {
          net.list <- nets[[k]]
          net <- net.list$data
          
          if('opacity' %in% names(net.list)) {
            arrow.opacity <- net.list$opacity
          } else {
            arrow.opacity <- default.opacity
          }
          
          if('color' %in% names(net.list)) {
            arrow.col <- net.list$color
          } else {
            arrow.col <- x@meta$plot.arrow.default$color[k]
          }
          arrow.col <- rgb(t(col2rgb(arrow.col)), alpha = round(arrow.opacity*255), maxColorValue = 255)
          arrow.col.list <- c(arrow.col.list, arrow.col)
          
          if('shift.x' %in% names(net.list)) {
            arrow.shift.x <- net.list$shift.x
          } else {
            arrow.shift.x <- x@meta$plot.arrow.default$shift.x[k]
          }
          arrow.shift.x.list <- c(arrow.shift.x.list, arrow.shift.x)
          
          if('shift.y' %in% names(net.list)) {
            arrow.shift.y <- net.list$shift.y
          } else {
            arrow.shift.y <- x@meta$plot.arrow.default$shift.y[k]
          }
          arrow.shift.y.list <- c(arrow.shift.x.list, arrow.shift.y)
          
          if('length.rate' %in% names(net.list)) {
            arrow.length.rate <- net.list$length.rate
          } else {
            arrow.length.rate <- default.length.rate
          }
          if('shorten' %in% names(net.list)) {
            arrow.shorten <- net.list$shorten
          } else {
            arrow.shorten <- default.shorten
          }
          
          if('label' %in% names(net.list)) {
            arrow.label.list[k] <- net.list$label
          }
          
          if('head.length' %in% names(net.list)) {
            arrow.head.length <- net.list$head.length
          } else {
            arrow.head.length <- default.head.length
          }
          
          if('head.type' %in% names(net.list)) {
            arrow.head.type <- net.list$head.type
          } else {
            arrow.head.type <- default.head.type
          }
          
          if('thickness' %in% names(net.list)) {
            arrow.thickness <- net.list$thickness
          } else {
            arrow.thickness <- default.thickness
          }
          
          for (i in dimnames(net)[[1]]){
            for(j in dimnames(net)[[2]]) {
              if (i != j){
                if (net[i,j] > 0) {
                  arrow.start <- coord[seats.which[i],]
                  arrow.stop <- coord[seats.which[j],]
                  arrow.coords <- .arrow.resize(
                    arrow.start[1],
                    arrow.start[2],
                    arrow.stop[1],
                    arrow.stop[2],
                    size = arrow.length.rate
                  )
                  arrow.coords <- .arrow.cut(
                    x0 = arrow.coords['x0'],
                    y0 = arrow.coords['y0'],
                    x1 = arrow.coords['x1'],
                    y1 = arrow.coords['y1'],
                    cut = arrow.shorten
                  )
                  #                   print(arrow.coords)
                  Arrows(
                    x0 = arrow.coords['x0'] + arrow.shift.x,
                    y0 = arrow.coords['y0'] + arrow.shift.y,
                    x1 = arrow.coords['x1'] + arrow.shift.x,
                    y1 = arrow.coords['y1'] + arrow.shift.y,
                    col=arrow.col,
                    arr.col=arrow.col,
                    arr.length=arrow.head.length,
                    lwd=net[i,j] * arrow.thickness,
                    arr.type = arrow.head.type
                  )
                }
              }
            }
          }
        }
      }
      
      
      
      
      ## LEGEND
      leg.pring = graph.legend.list(x)$print
      leg.cex = graph.legend.list(x)$cex
      leg.ncol = graph.legend.list(x)$ncol
      leg.horiz = graph.legend.list(x)$horiz
      leg.lwd = graph.legend.list(x)$lwd
      
      par(graph.par.list(x))
      plot.new()
      if(leg.pring) {
        if(flag.color) {
          legend(
            x = "topleft",
            legend = names(color$legend),
            fill = color$legend,
            bty = 'n',
            cex = leg.cex,
            ncol = leg.ncol,
            horiz = leg.horiz
          )
        }
        if(flag.symbol) {
          legend(
            x = "top",
            legend = names(x@plot.symbol$legend),
            pch = .graph.symbol.list[match(x@plot.symbol$legend, names(.graph.symbol.list))],
            col = symb.color,
            bty = 'n',
            cex = leg.cex
          )
        }
        if(flag.arrow) {
          legend(
            x = "topright",
            legend = arrow.label.list,
            col = arrow.col.list,
            lty = 1,
            bty = 'n',
            cex = leg.cex,
            lwd = leg.lwd
          )
        }
        #     if(flag.arrow){
        #       dn <- x@plot.arrow$legend
        #       if(nchar(dn) > 0) {
        #         xx <- 0.8
        #         yy <- 0.66
        #         tt <- 0.04
        #         arrows(xx, yy, xx + tt, yy, col=arrow.col, length=arrow.length.rate)
        #         
        #         text(x= xx+tt+0.11, y = yy, labels = dn, ...=...)
        #       }
        #     }
      }
    }
    
    reset = TRUE
    if(is.element('reset',names(lay)))
      if(lay$reset == FALSE)
        reset = FALSE
    if(reset)
      par(def.par)  #- reset to default
  }
)


# ---------------------------------------------------------------
# Symbol list
# ---------------------------------------------------------------
.graph.symbol.list <- c(
  "circle" = 1,
  "triangle.up" = 2,
  "triangle.down" = 6,
  "square" = 22,
  "square.rotated" = 5,
  "square.triangle.down" = 14,
  "cross" = 3,
  "times" = 4  
)

#' @importFrom graphics text
plot.symbol.list <- function(){
  l <- .graph.symbol.list
  coord <- 2:(length(l)+1)
  plot(
    coord,
    pch = l,
    ylim = c(1, length(l)+1),
    xlim = c(0, length(l)+1),
    xaxt = 'n',
    xlab = '',
    yaxt = 'n',
    ylab = '',
    bty = "n"
  )
  text(
    coord - 0.5,
    names(l)
  )
}
# plot.symbol.list()
.plot.symbol.list.all <- function(){
  n.symbols <- 25
  plot(
    1:n.symbols,
    pch = 1:n.symbols
  )
}
# .plot.symbol.list.all()
.arrow.resize <- function(x0, y0, x1, y1, size = 0.9) {
  stopifnot(size > 0 && size <= 1)
  
  segment <- c(0,1)
  breaks.size <- (1-size)/2
  breaks <- c(segment[1] + breaks.size, segment[2] - breaks.size)
  
  new.x0 <- (1-breaks[1]) * x0 + breaks[1] * x1
  new.y0 <- (1-breaks[1]) * y0 + breaks[1] * y1
  new.x1 <- (1-breaks[2]) * x0 + breaks[2] * x1
  new.y1 <- (1-breaks[2]) * y0 + breaks[2] * y1
  
  out <- c(new.x0, new.y0, new.x1, new.y1)
  names(out) <- c('x0', 'y0', 'x1', 'y1')
  return(out)
}
# .arrow.resize(0,0,1,1)
# .arrow.resize(0,0,1,1)['x0']
.arrow.cut <- function(x0, y0, x1, y1, cut = 0.2) {
  
  arrow.norm <- sqrt((x1-x0)^2 + (y1-y0)^2)
  unitaire.x <- (x1-x0)/arrow.norm
  unitaire.y <- (y1-y0)/arrow.norm
  
  new.x0 <- x0 + unitaire.x * cut
  new.y0 <- y0 + unitaire.y * cut
  new.x1 <- x1 - unitaire.x * cut
  new.y1 <- y1 - unitaire.y * cut
  
  out <- c(new.x0, new.y0, new.x1, new.y1)
  names(out) <- c('x0', 'y0', 'x1', 'y1')
  return(out)
}
# .arrow.cut('x0' = 0,'y0' = 0,'x1' = 1,'y1' = 1)
# .arrow.cut(0,0,2,2)
# .arrow.cut(1,1,2,2)
# .arrow.resize(0,0,1,1)['x0']


# .base1to256 <- function(x) {
#   return(round(x*255))
# }
# 
# .base256tohex <- function(x) {
#   sprintf("%X", x) 
# }
# .base256tohex(0)
# .base256tohex(10)
# .base256tohex(255)
# .base1tohex <- function(x) {
#   return(.base256tohex(.base1to256(x)))
# }
# .addzero <- function(x) {
#   if(nchar(x) == 1)
#     x <- paste('0', x, sep = '')
#   
#   return(x)
#   }
# }
# .addzero(.base1tohex(0))
# .addzero(.base1tohex(0.5))
# .addzero(.base1tohex(1))