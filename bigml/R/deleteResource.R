#' Deleting BigML Resources
#' @export 
#' @param resource_id the resource to delete.
#' @template author
#' @details This function deletes bigml resources referenced by their resource
#'	id.
#' @template dots
#' @return TRUE if successful, FALSE otherwise.
#' @examples 
#' \dontrun{
#' # replace with your valid credentials:
#'	deleteResource("source/1")
#' }
deleteResource <-
function (resource_id, ...){
    res = ''
    if (class(resource_id) == "list") {
        if ("resource" %in% names(resource_id)) {
            res = resource_id$resource
        }
        else {
            stop("argument doesn't apear to be a valid BigML response")
        }
    }
    else if (class(resource_id) == "character") {
        res = resource_id
    } else {
        stop("wrong type for resource_id")
    }

    vals = strsplit(res, '/')
    vals = vals[[1]]
    resources = c("source","dataset","model","prediction")
    if (! vals[1] %in% resources || length(vals) != 2){
        stop("illegal resource id")
    }
    url = paste(.BIGML_URL,vals[1],sep='')
    response = .basic_api(url)$delete(res, ...)
    return(response =='')
}
