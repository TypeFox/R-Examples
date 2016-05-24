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

#' @title An RDF Parser object
#' @description The Parser class provides methods to parse RDF content into a Redland
#' RDF model.
#' @slot librdf_parser A redland parser object
#' @rdname Parser-class
#' @aliases Parser
#' @include redland.R
#' @include World.R
#' @include Model.R
#' @keywords classes
#' @export
#' @section Methods:
#' \itemize{
#'   \item{\code{\link{Parser-initialize}}}{: Initialize a Parser object.}
#'   \item{\code{\link{parseFileIntoModel}}}{: Parse the contents of a file into a model.}
#'   \item{\code{\link{freeParser}}}{: Free memory used by a librdf parser.}
#' }
#' @seealso \code{\link{redland}}{: redland package}
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' # Create the default "rdfxml" parser
#' parser <- new("Parser", world)
#' filePath <- system.file("extdata/example.rdf", package="redland")
#' parseFileIntoModel(parser, world, filePath, model)
setClass("Parser", slots = c(librdf_parser = "_p_librdf_parser_s"))

#' Initialize a Parser object.
#' @description A Parser object is initialized for a specific RDF serialization.
#' @details The serialization format that are supported by 
#' @rdname Parser-initialize
#' @aliases Parser-initialize
#' @param .Object the Parser object
#' @param world a World object
#' @param name name of the parser factory to use
#' @param mimeType a mime type of the syntax of the model
#' @param typeUri a URI for the syntax of the model
#' @return the Parser object
#' @export
setMethod("initialize", signature = "Parser", definition = function(.Object, world, name="rdfxml", mimeType="application/rdf+xml", typeUri=as.character(NA)) {
  # Ensure that all provided params are not null
  stopifnot(!is.null(world))
  
  if(is.na(typeUri)) {
    librdf_uri <- librdf_new_uri(world@librdf_world, "")
  } else {
    librdf_uri <- librdf_new_uri(world@librdf_world, typeUri)
  }
  
  # Create the underlying redland statement object
  .Object@librdf_parser <- librdf_new_parser(world@librdf_world,
                                             name,
                                             mimeType,
                                             librdf_uri);
  return(.Object)
})

#' Parse the contents of a file into a model
#' @description The contents of a the specified file are read and parsed into the initialized
#' Parser object
#' @details The parser factory name specified during initialization determines how the content is
#' parsed, for example, if 'rdfxml' was specified during parser initialization, then the parser
#' expects RDF/XML content as specified in the W3C recommendation (http://www.we3.org/TR/REC-rdf-syntax)
#' @rdname parseFileIntoModel
#' @param .Object a Parser object 
#' @param world a World object
#' @param filePath a file that contains the RDF content
#' @param model a Model object to parse the RDF content into
#' @param ... (Additional parameters)
#' @examples 
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' # Create the default "rdfxml" parser
#' parser <- new("Parser", world)
#' filePath <- system.file("extdata/example.rdf", package="redland")
#' parseFileIntoModel(parser, world, filePath, model)
#' @export
setGeneric("parseFileIntoModel", function(.Object, world, filePath, model, ...) {
  standardGeneric("parseFileIntoModel")
})

#' @rdname parseFileIntoModel
#' @param baseUri a base URI (i.e. XML base) to apply to the model
setMethod("parseFileIntoModel", signature("Parser", "World", "character", "Model"), function(.Object, world, filePath, model, baseUri=as.character(NA)) {
  stopifnot(!is.null(model))
  
  if(is.na(baseUri)) {
    librdf_uri <- NULL
  } else {
    librdf_uri <- librdf_new_uri(world@librdf_world, baseUri)
  }
  
  # Construct a file URI of the input file name
  err <- try(absFilePath <- normalizePath(filePath, mustWork = TRUE))
  stopifnot (!class(err) == "try-error")
  
  # Absolute paths on Windows require leading slash: /C:/foo/bar
  if(grepl("^[a-zA-Z]:", absFilePath))
    absFilePath <- paste0("/", absFilePath)
  
  fileUri = sprintf("file://%s", absFilePath)
  contentUri <- librdf_new_uri(world@librdf_world, fileUri)
  status <- librdf_parser_parse_into_model(.Object@librdf_parser, contentUri, NULL, model@librdf_model)
})

#' Free memory used by a librdf parser
#' @details After freeNode is called, the Node object is no longer usable and should
#' be deleted  \code{"rm(nodeName)"} and a new object created.
#' @rdname freeParser
#' @param .Object a Node object
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' parser <- new("Parser", world)
#' filePath <- system.file("extdata/example.rdf", package="redland")
#' parseFileIntoModel(parser, world, filePath, model)
#' # At this point, some operations would be performed with the Model that has been populated
#' # with the parser.
#' # See '?redland' for a complete example.
#' # When the parser object is no longer needed, the resources it had allocated can be freed.
#' freeParser(parser)
#' rm(parser)
#' @export
setGeneric("freeParser", function(.Object) {
  standardGeneric("freeParser")
})

#' @rdname freeParser
setMethod("freeParser", signature("Parser"), function(.Object) {
  librdf_free_parser(.Object@librdf_parser)
})
