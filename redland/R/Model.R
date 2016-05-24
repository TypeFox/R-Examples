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
#' @title A Redland Model object
#' @description A Model object is used to store the statements (triples) of an RDF model.
#' @details A Model may be created manually by creating \code{\link{Statement}} and adding
#' them to the Model using \code{\link{addStatement}}, or a Model may be read in from a
#' previously saved file using \code{\link{parseFileIntoModel}}. Once a Model is created,
#' it can be queried using \code{\link{Query}}.
#' @seealso View examples of creating models by viewing the \code{'redland_overview'} vignette: \code{'vignette("redland_overview")'}
#' @slot librdf_model A redland model object
#' @rdname Model-class
#' @aliases Model
#' @include redland.R
#' @include World.R
#' @include Storage.R
#' @include Statement.R
#' @keywords classes
#' @section Methods:
#' \itemize{
#'   \item{\code{\link{Model-initialize}}}{: Initialize a Model object}
#'   \item{\code{\link{addStatement}}}{: Add a Statement object to the Model}
#'   \item{\code{\link{freeModel}}}{: Free memory used by a librdf model object}
#' }
#' @seealso \code{\link{redland}}{: redland package}
#' @export
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
setClass("Model", slots=c(librdf_model = "_p_librdf_model_s"))

#' Constructor for a Model object.
#' @param .Object a Node object
#' @param world a World object
#' @param storage a Storage object
#' @param options extra options for model initialization
#' @rdname Model-initialize
#' @aliases Model-initialize
#' @return the World object
#' @export
setMethod("initialize", signature = "Model", definition = function(.Object, world, storage, options) {
    stopifnot(!is.null(world), !is.null(storage))
    .Object@librdf_model <- librdf_new_model(world@librdf_world, storage@librdf_storage, options);
    return(.Object)
})

#' Add a Statement object to the Model
#' @rdname addStatement
#' @param .Object a Model object
#' @param statement the Statement that will be added
#' @export 
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
setGeneric("addStatement", function(.Object, statement) {
  standardGeneric("addStatement")
})

#' @rdname addStatement
setMethod("addStatement", signature("Model", "Statement"), function(.Object, statement) {
  
  librdf_model_add_statement(.Object@librdf_model, statement@librdf_statement);
})

#' Free memory used by a librdf model.
#' @rdname freeModel
#' @details After this method is called, the Model object is no longer usable and should
#' be deleted \code{"rm(model)"} and a new object created.
#' @param .Object a Model object
#' @export
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' # At this point, some operations would be performed with the model.
#' # See '?redland' for a complete example.
#' # When the Model object is no longer needed, the resources it has allocated can be freed.
#' freeModel(model)
#' rm(model)
setGeneric("freeModel", function(.Object) {
  standardGeneric("freeModel")
})

#' @rdname freeModel
setMethod("freeModel", signature("Model"), function(.Object) {
  librdf_free_model(.Object@librdf_model)
})
