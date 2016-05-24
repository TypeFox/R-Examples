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

#' @title A Redland QueryResults object is used to inspect query results from a Query object.
#' @description The QueryResults object contains the RDF statements that were returned from
#' a query on an RDF model.
#' @slot librdf_query_results A redland query object
#' @rdname QueryResults-class
#' @aliases QueryResults
#' @include redland.R
#' @keywords classes
#' @export
#' @section Methods:
#' \itemize{
#'   \item{\code{\link{QueryResults-initialize}}}{: Initialize a QueryResults object.}
#'   \item{\code{\link{getNextResult}}}{: Get the next query result.}
#'   \item{\code{\link{freeQueryResults}}}{: Free memory used by a librdf query result.}
#' }
#' @seealso \code{\link{redland}}{: redland package}
#' @examples
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' stmt <- new("Statement", world=world, 
#'   subject="https://cn.dataone.org/cn/v1/resolve/urn:uuid:274a0c5c-3082-4562-bbd3-2b1288768cac",
#'   predicate="http://www.w3.org/ns/prov#hadPlan",
#'   object="https://cn.dataone.org/cn/v1/resolve/urn:uuid:01305f45-f22b-40c8-8d27-00357d01e4a5")
#' status <- addStatement(model, stmt)
#' stmt <- new("Statement", world=world, 
#'   subject="https://orcid.org/0000-0002-2192-403X",
#'   predicate="http://www.w3.org/ns/prov#Agent",
#'   object="slaughter", 
#'   objectType="literal", datatype_uri="http://www.w3.org/2001/XMLSchema#string")
#' status <- addStatement(model, stmt)
#' queryString <- paste("PREFIX orcid: <https://orcid.org/>",
#'                      "PREFIX dataone: <https://cn.dataone.org/cn/v1/resolve/>",
#'                      "PREFIX prov: <http://www.w3.org/ns/prov#>",
#'                      "SELECT ?a ?c WHERE { ?a prov:Agent ?c . }", sep=" ")
#' query <- new("Query", world, queryString, base_uri=NULL, 
#'   query_language="sparql", query_uri=NULL)
#' queryResult <- executeQuery(query, model)
#' result <- getNextResult(queryResult)
setClass("QueryResults", slots = c(librdf_query_results = "_p_librdf_query_results"))

#' Initialize the QueryResults object.
#' @description The QueryResults object is initialized with the librdf query result from
#' return value of \code{'Query.execute()'}.
#' @details A QueryResults object is returned by the \code{Query.executeQuery()} method, so typically a user
#' does not initialize a QueryResult object by calling \code{new("QueryResult", ...)}
#' @rdname QueryResults-initialize
#' @aliases QueryResults-initialize
#' @param .Object the QueryResults object.
#' @param results a librdf query result
#' @return the QueryResults object
#' @export
setMethod("initialize", signature = "QueryResults", definition = function(.Object, results) {
  .Object@librdf_query_results <- results     
  return(.Object)
})

#' Get the next query result.
#' @description The next query result is returned. .
#' @rdname getNextResult
#' @param .Object a QueryResults object
#' @export
setGeneric("getNextResult", function(.Object) {
  standardGeneric("getNextResult")
})

#' @rdname getNextResult
setMethod("getNextResult", signature("QueryResults"), function(.Object) {

  nodeNames <- list()
  nodeValues <- list()
  
  # Process the next result, storing the bound values in a list
  if (!is.null(.Object@librdf_query_results) && librdf_query_results_finished(.Object@librdf_query_results) == 0) {
    num_nodes <- librdf_query_results_get_bindings_count(.Object@librdf_query_results)
    for (i in 1:num_nodes-1) {
      binding_name <- librdf_query_results_get_binding_name(.Object@librdf_query_results, i)
      val = librdf_query_results_get_binding_value(.Object@librdf_query_results, i)
      # If no value returned for this binding, set to "NA"
      if (!is.null.externalptr(val@ref)) {
        nval <- librdf_node_to_string(val)
      } else {
        nval = as.character(NA)
      }
      nodeNames <- c(nodeNames, binding_name)
      nodeValues <- c(nodeValues, nval)
    }
  } else {
    return(NULL)
  }
  
  names(nodeValues) <- nodeNames
  nextResult <-librdf_query_results_next(.Object@librdf_query_results)
  
  return(nodeValues)
})

#' Free memory used by a librdf query results
#' @description After this method is called, the QueryResults object is no longer usable and should
#' be deleted with \code{"rm(query)"}.
#' @rdname freeQueryResults
#' @param .Object a QueryResults object
#' @examples 
#' world <- new("World")
#' storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
#' model <- new("Model", world, storage, options="")
#' stmt <- new("Statement", world=world, 
#'   subject="https://orcid.org/0000-0002-2192-403X",
#'   predicate="http://www.w3.org/ns/prov#Agent",
#'   object="slaughter", 
#'   objectType="literal", datatype_uri="http://www.w3.org/2001/XMLSchema#string")
#' status <- addStatement(model, stmt)
#' queryString <- paste("PREFIX orcid: <https://orcid.org/>",
#'                      "PREFIX dataone: <https://cn.dataone.org/cn/v1/resolve/>",
#'                      "PREFIX prov: <http://www.w3.org/ns/prov#>",
#'                      "SELECT ?a ?c WHERE { ?a prov:Agent ?c . }", sep=" ")
#' query <- new("Query", world, queryString, base_uri=NULL, 
#'   query_language="sparql", query_uri=NULL)
#' queryResult <- executeQuery(query, model)
#' result <- getNextResult(queryResult)
#' 
#' # When the queryResult is no longer needed, the resources it had allocated can be freed.
#' freeQueryResults(queryResult)
#' rm(queryResult)
#' @export
setGeneric("freeQueryResults", function(.Object) {
  standardGeneric("freeQueryResults")
})

#' @rdname freeQueryResults
setMethod("freeQueryResults", signature("QueryResults"), function(.Object) {
  # Have to free all of the nodes that were created by the query result that
  # hold the bound node values
  if (!is.null(.Object@librdf_query_results) && librdf_query_results_finished(.Object@librdf_query_results) == 0) {
    num_nodes <- librdf_query_results_get_bindings_count(.Object@librdf_query_results)
    for (i in 1:num_nodes-1) {
      binding_name <- librdf_query_results_get_binding_name(.Object@librdf_query_results, i)
      val = librdf_query_results_get_binding_value(.Object@librdf_query_results, i)    
      if (!is.null.externalptr(val@ref)) {
        librdf_free_node(val)
      } 
    }
  }
  
  librdf_free_query_results(.Object@librdf_query_results)
})
