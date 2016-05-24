#
#   This work was created by the National Center for Ecological Analysis and Synthesis.
#
#     Copyright 2015 Regents of the University of California
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#' @title An RDF Serializer object.
#' @description The Serializer class provides methods to convert a Model object
#' to other forms, for example, write out a Model to a file.
#' @slot librdf_serializer A redland statement object
#' @rdname Serializer-class
#' @aliases Serializer
#' @include redland.R
#' @include World.R
#' @include Model.R
#' @keywords classes
#' @export
#' @section Methods:
#' \itemize{
#'   \item{\code{\link{Serializer-initialize}}}{: Initialize a Serializer object.}
#'   \item{\code{\link{setNameSpace}}}{: Set a namespace for the serializer.}
#'   \item{\code{\link{serializeToCharacter}}}{: Serialize a model to a character vector.}
#'   \item{\code{\link{serializeToFile}}}{: Serialize a model to a file.}
#'   \item{\code{\link{freeSerializer}}}{: Free memory used by a librdf serializer.}
#' }
#' @seealso \code{\link{redland}}{: redland package}
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' filePath <- system.file("extdata/example.rdf", package="redland")
#' parser <- new("Parser", world)
#' parseFileIntoModel(parser, world, filePath, model)
#' # Creat the default "rdfxml" serizlizer
#' serializer <- new("Serializer", world)
#' # Add a namespace definition to the serializer
#' status <- setNameSpace(serializer, world, namespace="http://purl.org/dc/elements/1.1/", prefix="dc")
#' rdf <- serializeToCharacter(serializer, world, model, baseUri="")
setClass("Serializer", slots = c(librdf_serializer = "_p_librdf_serializer_s"))

#' Construct a Serializer object.
#' @rdname Serializer-initialize
#' @aliases Serializer-initialize
#' @param .Object the Serializer object
#' @param world a World object
#' @param name name of a previously created serializer factory to use
#' @param mimeType a mime type of the syntax of the model
#' @param typeUri a URI for the syntax of the model
#' @return the Serializer object
#' @export
setMethod("initialize", signature = "Serializer", definition = function(.Object, world, name="rdfxml", mimeType="application/rdf+xml", typeUri=as.character(NA)) {
  # Ensure that all provided params are not null
  stopifnot(!is.null(world))
  
  if(is.na(typeUri)) {
    librdf_uri <- NULL
  } else {
    librdf_uri <- librdf_new_uri(world@librdf_world, typeUri)
  }
  
  # Create the underlying redland statement object
  .Object@librdf_serializer <- librdf_new_serializer(world@librdf_world, 
                                                     name,
                                                     mimeType,
                                                     librdf_uri);
  return(.Object)
})

#' Set a namespace for the serializer.
#' @rdname setNameSpace
#' @param .Object a Serializer object
#' @param world a World object
#' @param namespace the namespace to add to the serializer
#' @param prefix the namespace prefix to associate with the namespace
#' @export
setGeneric("setNameSpace", function(.Object, world, namespace, prefix) {
  standardGeneric("setNameSpace")
})

#' @rdname setNameSpace
setMethod("setNameSpace", signature("Serializer", "World", "character", "character"), function(.Object, world, namespace, prefix) {
  
  stopifnot(!is.null(world))
  
  librdf_uri <- librdf_new_uri(world@librdf_world, namespace)
  librdf_serializer_set_namespace(.Object@librdf_serializer, librdf_uri, prefix)
})

#' Serialize a model to a character vector.
#' @rdname serializeToCharacter
#' @param .Object a Serializer object
#' @param world a World object
#' @param model a Model object
#' @param ... Additional parameters
#' @return a character vector containing the serialized model
#' @export
setGeneric("serializeToCharacter", function(.Object, world, model, ...) {
  standardGeneric("serializeToCharacter")
})

#' @rdname serializeToCharacter
#' @param baseUri a URI to prepend to relative URIs in the document
setMethod("serializeToCharacter", signature("Serializer", "World", "Model"), function(.Object, world, model, baseUri=as.character(NA)) {
  
  stopifnot(!is.null(world))
  stopifnot(!is.null(model))
  
  # Convert baseUri to a librdf_uri
  if(is.na(baseUri)) {
    librdf_uri <- NULL
  } else {
    librdf_uri <- librdf_new_uri(world@librdf_world, baseUri)
  }
  
  RDFstring <- librdf_serializer_serialize_model_to_string(.Object@librdf_serializer, librdf_uri, model@librdf_model)
  return(RDFstring)
})

#' Serialize a model to a file.
#' @rdname serializeToFile
#' @param .Object a Serializer object
#' @param world a World object
#' @param model a Model object
#' @param filePath a file path that the serialized model will be written to
#' @param ... Additional parameters
#' @return an integer containing the return status where non zero indicates an error occured during serialization
#' @export
setGeneric("serializeToFile", function(.Object, world, model, filePath, ...) {
  standardGeneric("serializeToFile")
})

#' @rdname serializeToFile
#' @param baseUri a base URI to use for the serialization
setMethod("serializeToFile", signature("Serializer", "World", "Model", "character"), function(.Object, world, model, filePath, baseUri=as.character(NA)) {
  
  stopifnot(!is.null(world))
  stopifnot(!is.null(model))
  stopifnot(!is.null(filePath))
  
  # Convert baseUri to a librdf_uri
  if(is.na(baseUri)) {
    librdf_uri <- NULL
  } else {
    librdf_uri <- librdf_new_uri(world@librdf_world, baseUri)
  }
  
  status <-librdf_serializer_serialize_model_to_file (.Object@librdf_serializer, filePath, librdf_uri, model@librdf_model);
  return(status)
})

#' Free memory used by a librdf serializer.
#' @rdname freeSerializer
#' @param .Object a Serializer object
#' @examples 
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' filePath <- system.file("extdata/example.rdf", package="redland")
#' parser <- new("Parser", world)
#' parseFileIntoModel(parser, world, filePath, model)
#' # Creat the default "rdfxml" serizlizer
#' serializer <- new("Serializer", world)
#' # At this point, some operations would be performed with the Serializer object. 
#' # See '?Serializer' for a complete example.
#' # When the serializer object is no longer needed, the resources it had allocated can be freed.
#' freeSerializer(serializer)
#' rm(serializer)
#' @export
setGeneric("freeSerializer", function(.Object) {
  standardGeneric("freeSerializer")
})

#' @rdname freeSerializer
setMethod("freeSerializer", signature("Serializer"), function(.Object) {
  librdf_free_serializer(.Object@librdf_serializer)
})
