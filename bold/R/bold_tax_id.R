#' Search BOLD for taxonomy data by BOLD ID.
#'
#' @export
#' @param id (integer) One or more BOLD taxonomic identifiers. required.
#' @param dataTypes (character) Specifies the datatypes that will be returned. 'all' returns all
#' data. 'basic' returns basic taxon information. 'images' returns specimen images.
#' @param includeTree (logical) If TRUE (default: FALSE), returns a list containing information
#' for parent taxa as well as the specified taxon.
#' @template otherargs
#' @references \url{http://boldsystems.org/index.php/resources/api?type=taxonomy#idParameters}
#' @seealso \code{bold_tax_name}
#' @examples \dontrun{
#' bold_tax_id(id=88899)
#' bold_tax_id(id=88899, includeTree=TRUE)
#' bold_tax_id(id=88899, includeTree=TRUE, dataTypes = "stats")
#' bold_tax_id(id=c(88899,125295))
#'
#' ## dataTypes parameter
#' bold_tax_id(id=88899, dataTypes = "basic")
#' bold_tax_id(id=88899, dataTypes = "stats")
#' bold_tax_id(id=88899, dataTypes = "images")
#' bold_tax_id(id=88899, dataTypes = "geo")
#' bold_tax_id(id=88899, dataTypes = "sequencinglabs")
#' bold_tax_id(id=88899, dataTypes = "depository")
#' bold_tax_id(id=88899, dataTypes = "thirdparty")
#' bold_tax_id(id=88899, dataTypes = "all")
#' bold_tax_id(id=c(88899,125295), dataTypes = "geo")
#' bold_tax_id(id=c(88899,125295), dataTypes = "images")
#'
#' ## Passing in NA
#' bold_tax_id(id = NA)
#' bold_tax_id(id = c(88899,125295,NA))
#'
#' ## get httr response object only
#' bold_tax_id(id=88899, response=TRUE)
#' bold_tax_id(id=c(88899,125295), response=TRUE)
#'
#' ## curl debugging
#' library('httr')
#' bold_tax_id(id=88899, config=verbose())
#' }

bold_tax_id <- function(id, dataTypes='basic', includeTree=FALSE, response=FALSE, ...) {
  
  tmp <- lapply(id, function(x)
    get_response(args = bc(list(taxId = x, dataTypes = dataTypes, includeTree = if (includeTree) TRUE else NULL)),
                 url = paste0(bbase(), "API_Tax/TaxonData"), ...)
  )
  if (response) { 
    tmp 
  } else {
    res <- do.call(rbind.fill, Map(process_response, x = tmp, y = id, z = includeTree, w = dataTypes))
    if (NCOL(res) == 1) { 
      res$noresults <- NA
      return(res)
    } else { 
      res 
    }
  }
}
