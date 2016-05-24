#' Get Dryad metrics from Total Impact (http://totalimpact.org/) API.
#' @import stringr RJSONIO plyr
#' @param id The Dryad identifier - see examples.
#' @param sleep Time (in seconds) before function sends API call - defaults to
#'    zero.  Set to higher number if you are using this function in a loop with
#'    many API calls.  Although, with TotalImpact, you can also just put in
#'    multiple id's in one call (much faster) instead of using a loop (slower).
#' @param url The base URL (do not change from default).
#' @return A data.frame, or list of data.frame's, of results.
#' @export
#' @examples \dontrun{
#' totimp_dryad(id = '10.5061/dryad.8671')
#' totimp_dryad(id = list('10.5061/dryad.8671', '10.5061/dryad.8908'))
#' }

totimp_dryad <- function(id = NA, sleep = 0, url = "http://total-impact.org/api/v1/items/") {
    
    Sys.sleep(sleep)
    
    unix2POSIXct <- function(time) structure(time, class = c("POSIXt", 
        "POSIXct"))
    
    id_ <- paste(laply(id, str_replace_all, pattern = "/", replacement = "%252F"), 
        collapse = ",")
    url2 <- paste(url, id_, ".json", "?fields=metrics", sep = "")
    out <- fromJSON(url2)
    getdf <- function(x) {
        out_ <- list()
        for (i in 1:length(x[[1]])) {
            df <- rbind(data.frame(x[[1]][[i]]$metrics[[1]][1:3]), data.frame(x[[1]][[i]]$metrics[[2]][1:3]), 
                data.frame(x[[1]][[i]]$metrics[[3]][1:3]), data.frame(x[[1]][[i]]$metrics[[4]][1:3]))
            # This is a comment
            df$last_update_converted <- ldply(df$last_update, unix2POSIXct)[, 
                1]
            out_[[x[[1]][[i]]$id]] <- df
        }
        out_
    }
    getdf(out)
} 
