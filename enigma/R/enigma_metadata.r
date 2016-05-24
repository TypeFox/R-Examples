#' Search for metadata on a dataset from Enigma.
#'
#' @export
#'
#' @param dataset Dataset name. Required.
#' @param key (character) Required. An Enigma API key. Supply in the function call, or store in
#' your \code{.Rprofile} file, or do \code{options(enigmaKey = "<your key>")}. Obtain an API key
#' by creating an account with Enigma at \url{http://enigma.io}, then obtain an API key from
#' your account page.
#' @param ... Named options passed on to \code{\link[httr]{GET}}
#' @details Notice when you run the examples that the format of output is different for the
#' "parent nodes" vs. the "table nodes". Where the parent nodes have ouput$meta slots for
#' paths, immediate nodes and children tables, while the table nodes have ouput$meta slots for
#' info, table, ancestor datapaths, database boundary datapath, database boundary label, database
#' boundary tables, and paths, and an additional slot for description of table column attributes.
#' @examples \dontrun{
#' # After obtaining an API key from Enigma's website, pass in your key to the function call
#' # or set in your options (see above instructions for the key parameter)
#' # If you pass in your key to the function call use the key parameter
#'
#' # Parent node response attributes from Whitehouse dataset
#' ## US white house
#' enigma_metadata(dataset='us.gov.whitehouse')
#'
#' ## UCLA Ethnic power relations dataset
#' enigma_metadata(dataset='edu.ucla.epr')
#'
#' # Table node response attributes
#' ## US white house visitor list table
#' enigma_metadata(dataset='us.gov.whitehouse.visitor-list')
#'
#' ## UCLA Ethnic power relations dataset - armed conflict table
#' enigma_metadata('edu.ucla.epr.ethnic-armed-conflict')
#' }

enigma_metadata <- function(dataset=NULL, key=NULL, ...) {
  key <- check_key(key)
  check_dataset(dataset)

  url <- sprintf('%s/meta/%s/%s', en_base(), key, dataset)
  json <- enigma_GET(url, list(), ...)
  meta <- process_meta(json)
  result_names <- names(json$result)
  if(any(result_names %in% "columns")){
    out <- list(success = json$success, datapath = json$datapath, info = meta, columns = process_cols(json))
  } else {
    out <- list(success = json$success, datapath = json$datapath, info = meta)
  }
  structure(out, class="enigma_meta")
}

process_meta <- function(x){
  result_names <- names(x$result)
  if(!any(result_names %in% "columns")){
    list(paths=lapply(x$result$path, as.list),
         immediate_nodes=lapply(x$result$immediate_nodes, as.list),
         children_tables=lapply(x$result$children_tables, as.list))
  } else {
    tmp <- x$result$metadata
    tablemeta <- sapply(tmp, "[[", "value")
    names(tablemeta) <- c('total_rows','last_updated')
    list(info=as.list(x$info),
         table=as.list(tablemeta),
         ancestor_datapaths=x$result$ancestor_datapaths,
         db_boundary_datapath=x$result$db_boundary_datapath,
         db_boundary_label=x$result$db_boundary_label,
         db_boundary_tables=lapply(x$result$db_boundary_tables, as.list),
         paths=lapply(x$result$path, as.list))
  }
}

process_cols <- function(x){
  columns <- x$result$columns
  colsdat <- do.call(rbind.fill,
                     lapply(columns, function(x){
                       x[sapply(x, is.null)] <- NA; data.frame(x, stringsAsFactors = FALSE)
                     }))
  colsdat_table <- data.frame(id=colsdat$id,label=colsdat$label,type=colsdat$type,index=colsdat$index)
  list(table=colsdat_table, description=as.list(colsdat$description))
}
