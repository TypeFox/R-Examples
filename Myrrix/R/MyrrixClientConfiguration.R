#' @title Object of class MyrrixClientConfiguration
#' @description An object of class MyrrixClientConfiguration describes how access options to the ClientRecommender.
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{config}:}{A java object representing a MyrrixClientConfiguration}
#'  }
#' @name MyrrixClientConfiguration-class 
#' @rdname MyrrixClientConfiguration-class
#' @aliases MyrrixClientConfiguration-class
#' @exportClass MyrrixClientConfiguration
#' @importFrom rJava .jnew
#' @examples
#' myconfig <- new("MyrrixClientConfiguration")
#' myconfig
setClass(Class="MyrrixClientConfiguration",
         representation=representation(config="jobjRef"))
setMethod(f="initialize", signature="MyrrixClientConfiguration",
          definition = function(.Object) {
            .Object@config <- .jnew("net/myrrix/client/MyrrixClientConfiguration")
            .Object
          })



#' @title Methods to apply on objects of class MyrrixClientConfiguration
#' @description Methods to apply on objects of class MyrrixClientConfiguration. These allow to define which Myrrix service
#' to listen to. See the methods section of this doc for an enumerated list of function to apply.\cr
#' More information on http://myrrix.com/docs/client/java/javadoc/net/myrrix/client/MyrrixClientConfiguration.html
#'
#' @section Methods: 
#'  \describe{
#'    \item{\code{show(object)}:}{Show method for a MyrrixClientConfiguration object:
#'    Prints configuration settings, indicating which server information Myrrix talks to}
#'    \item{\code{getMyrrixOptions(object)}:}{Returns a list of configuration settings, indicating which information Myrrix talks to}
#'    \item{\code{setHost(object, host)}}{Sets the host of the Myrrix Serving layer}
#'    \item{\code{setPort(object, port)}}{Sets the port on which to access the Serving Layer}
#'    \item{\code{setUserName(object, userName)}}{Sets the user name needed to access the Serving Layer}
#'    \item{\code{setPassword(object, password)}}{Sets the password needed to access the Serving Layer}
#'    \item{\code{setContextPath(object, contextPath)}}{Sets the context path under which the target Serving Layer app is deployed}
#'    \item{\code{setKeystoreFile(object, keystoreFile)}}{Sets the the keystore file containing the server's SSL keys.}
#'    \item{\code{setKeystorePassword(object, keystorePassword)}}{Sets the the password for keystorefile}
#'    \item{\code{setSecure(object, secure)}}{Set if this client is accessing the Serving Layer over HTTPS, not HTTP}
#'    \item{\code{setAllPartitionsSpecification(object, allPartitionsSpecification)}}{Sets the specification for all servers that have partitions}
#'  }
#' @param object An object of class MyrrixClientConfiguration
#' @param host host containing the Serving Layer, if not in distributed mode
#' @param port port on which to access the Serving Layer, if not in distributed mode.
#' @param username user name needed to access the Serving Layer, if any
#' @param password password needed to access the Serving Layer, if any
#' @param allpartitionsspecification specification for all servers that have partitions. 
#' Only applicable in distributed mode and returns null otherwise. 
#' May be specified as "auto", in which case getHost() and getPort() must be valid, 
#' since this host will be queried for partition details. 
#' Otherwise, Serving Layers are specified explicitly as "host:port" pairs. 
#' Replicas are specified as many Serving Layers, separated by commas, 
#' like "rep1:port1,rep2:port2,...". Finally, partitions are specified as multiple 
#' replicas separated by semicolon, 
#' like "part1rep1:port11,part1rep2:port12;part2rep1:port21,part2rep2:port22;...". 
#' Example: "foo:80,foo2:8080;bar:8080;baz2:80,baz3:80"
#' @param contextpath the context path under which the target Serving Layer app is deployed (e.g. http://example.org/contextPath/...), 
#' or null if the default root context should be used.
#' @param keystorefile the keystore file containing the server's SSL keys. 
#' Only necessary when accessing a server with a temporary self-signed key, 
#' which is not by default trusted by the Java SSL implementation
#' @param keystorepassword password for keystorefile
#' @param secure if true, this client is accessing the Serving Layer over HTTPS, not HTTP
#' @param ... other arguments passed on the the methods
#' @name MyrrixClientConfiguration-methods
#' @rdname MyrrixClientConfiguration-methods
#' @aliases getMyrrixOptions getMyrrixOptions,MyrrixClientConfiguration-method 
#' @exportMethod getMyrrixOptions
#' @examples
#' myconfig <- new("MyrrixClientConfiguration")
#' getMyrrixOptions(myconfig)
#' setHost(myconfig, "myhostname")
#' setPort(myconfig, 20L)
setGeneric("getMyrrixOptions", function(object, ...) standardGeneric("getMyrrixOptions"))
setMethod(f="getMyrrixOptions", signature=signature(object="MyrrixClientConfiguration"),
          definition = function(object){
            config <- list()
            config$allPartitionsSpecification <- object@config$getAllPartitionsSpecification()
            config$contextPath <- object@config$getContextPath()
            config$host <- object@config$getHost()
            config$keystoreFile <- object@config$getKeystoreFile()
            config$keystorePassword <- object@config$getKeystorePassword()
            config$allPartitions <- object@config$getPartitions()
            config$password <- object@config$getPassword()
            config$port <- object@config$getPort()
            config$userName <- object@config$getUserName()
            config$secure <- object@config$isSecure()
            return(config)
          }
)
setMethod(f="show", signature=signature(object = "MyrrixClientConfiguration"),
          definition = function(object){
            print(getMyrrixOptions(object))
          }
)

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setHost setHost,MyrrixClientConfiguration,character-method 
#' @exportMethod setHost
setGeneric("setHost", function(object, host, ...) standardGeneric("setHost"))
setMethod("setHost", signature=signature(object = "MyrrixClientConfiguration", host="character"),
          definition = function(object, host){
             object@config$setHost(host)
             object
           })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setContextPath setContextPath,MyrrixClientConfiguration,character-method 
#' @exportMethod setContextPath
setGeneric("setContextPath", function(object, contextPath, ...) standardGeneric("setContextPath"))
setMethod("setContextPath", signature=signature(object = "MyrrixClientConfiguration", contextPath="character"),
          definition = function(object, contextPath){
            object@config$setContextPath(contextPath)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setKeystorePassword setKeystorePassword,MyrrixClientConfiguration,character-method 
#' @exportMethod setKeystorePassword
setGeneric("setKeystorePassword", function(object, keystorePassword, ...) standardGeneric("setKeystorePassword"))
setMethod("setKeystorePassword", signature=signature(object = "MyrrixClientConfiguration", keystorePassword="character"),
          definition = function(object, keystorePassword){
            object@config$setKeystorePassword(keystorePassword)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setPassword setPassword,MyrrixClientConfiguration,character-method 
#' @exportMethod setPassword
setGeneric("setPassword", function(object, password, ...) standardGeneric("setPassword"))
setMethod("setPassword", signature=signature(object = "MyrrixClientConfiguration", password="character"),
          definition = function(object, password){
            object@config$setPassword(password)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setPort setPort,MyrrixClientConfiguration,integer-method 
#' @exportMethod setPort
setGeneric("setPort", function(object, port, ...) standardGeneric("setPort"))
setMethod("setPort", signature=signature(object = "MyrrixClientConfiguration", port="integer"),
          definition = function(object, port){
            object@config$setPort(port)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setSecure setSecure,MyrrixClientConfiguration,logical-method 
#' @exportMethod setSecure
setGeneric("setSecure", function(object, secure, ...) standardGeneric("setSecure"))
setMethod("setSecure", signature=signature(object = "MyrrixClientConfiguration", secure="logical"),
          definition = function(object, secure){
            object@config$setSecure(secure)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setUserName setUserName,MyrrixClientConfiguration,character-method 
#' @exportMethod setUserName
setGeneric("setUserName", function(object, userName, ...) standardGeneric("setUserName"))
setMethod("setUserName", signature=signature(object = "MyrrixClientConfiguration", userName="character"),
          definition = function(object, userName){
            object@config$setUserName(userName)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setKeystoreFile setKeystoreFile,MyrrixClientConfiguration,character-method 
#' @exportMethod setKeystoreFile
setGeneric("setKeystoreFile", function(object, keystoreFile, ...) standardGeneric("setKeystoreFile"))
setMethod("setKeystoreFile", signature=signature(object = "MyrrixClientConfiguration", keystoreFile="character"),
          definition = function(object, keystoreFile){
            keystoreFile <- .jnew("java.io.File", keystoreFile)
            object@config$setKeystoreFile(keystoreFile)
            object
          })

#' @rdname MyrrixClientConfiguration-methods
#' @aliases setAllPartitionsSpecification setAllPartitionsSpecification,MyrrixClientConfiguration,character-method 
#' @exportMethod setAllPartitionsSpecification
setGeneric("setAllPartitionsSpecification", function(object, allPartitionsSpecification, ...) standardGeneric("setAllPartitionsSpecification"))
setMethod("setAllPartitionsSpecification", signature=signature(object = "MyrrixClientConfiguration", allPartitionsSpecification="character"),
          definition = function(object, allPartitionsSpecification){
            object@config$setAllPartitionsSpecification(allPartitionsSpecification)
            object
          })


