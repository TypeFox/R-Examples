#' Initialize polidata
#'
#' @export

InitPolidata <- function() {

    for(apiSource in apiSources) {
        reply <- ""
        apiKey <- GetAPIKey(apiSource)
        while(!IsValidReply(reply)) {
            reply <- PromptAPIKey(apiSource, apiKey)
        }
    }

    ShowAPIKeys()
    cat("API setting complete.")
}
