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
#' @title A Redland World object, used to initialize the Redland RDF library.
#' @description A World object is the top level object in the Redland RDF library
#' implementation, so it contains all other objects needed to process RDF
#' Models.
#' @slot librdf_world A redland world object
#' @aliases World
#' @include redland.R
#' @keywords classes
#' @useDynLib redland
#' @section Methods:
#' \itemize{
#'   \item{\code{\link{World-initialize}}}{: Initialize a World object}
#'   \item{\code{\link{freeWorld}}}{: Free memory used by a librdf world object}
#' }
#' @seealso \code{\link{redland}}{: redland package}
#' @examples
#' world <- new("World")
#' @import methods
#' @export
setClass("World", slots = c(librdf_world = "_p_librdf_world_s"))

#' Initialize the World object.
#' @rdname World-initialize
#' @aliases World-initialize
#' @param .Object the World object
#' @return the World object
#' @export
setMethod("initialize", signature = "World", definition = function(.Object) {
    .Object@librdf_world <- librdf_new_world();
    librdf_world_open(.Object@librdf_world)
    return(.Object)
})

#' Free memory used by a librdf world object
#' @details After this method is called, the World object is no longer usable and should
#' be deleted \code{"rm(world)"} and a new object created.
#' @rdname freeWorld
#' @param .Object a World object
#' @export
#' @examples
#' world <- new("World")
#' # At this point we would perform some operations using the world object.
#' # When the world object is no longer needed, we can free the resources it has allocated.
#' result <- freeWorld(world)
#' rm(world)
setGeneric("freeWorld", function(.Object) {
  standardGeneric("freeWorld")
})

#' @rdname freeWorld
setMethod("freeWorld", signature("World"), function(.Object) {
  librdf_free_world(.Object@librdf_world)
})
