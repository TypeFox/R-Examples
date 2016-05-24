#' @title Object of class ClientRecommender
#' @description An implementation of MyrrixRecommender which accesses a remote Serving Layer instance over HTTP or HTTPS. 
#' This is like a local "handle" on the remote recommender.
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{recommender}:}{A java object of class net.myrrix.client.ClientRecommender}
#'    \item{\code{clientConfiguration}:}{An object of class MyrrixClientConfiguration, which holds the java object with the 
#'    connection settings to Myrrix}
#'  }
#' @param config An object of class MyrrixClientConfiguration
#' @name ClientRecommender-class 
#' @rdname ClientRecommender-class
#' @aliases ClientRecommender-class
#' @exportClass ClientRecommender
#' @importFrom rJava .jnew
#' @usage new("ClientRecommender", config)
#' @examples
#' myconfig <- new("MyrrixClientConfiguration")
#' myconfig
#' recommendationengine <- new("ClientRecommender", config=myconfig)
#' recommendationengine
setClass(Class="ClientRecommender",
         representation=representation(recommender="jobjRef", clientConfiguration="MyrrixClientConfiguration"))
setMethod(f="initialize", signature="ClientRecommender",
          definition = function(.Object, config) {            
            .Object@recommender <- .jnew("net.myrrix.client.ClientRecommender", config@config)
            .Object@clientConfiguration <- config
            .Object
          })
