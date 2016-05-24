getStatus <- function(jobId, api.key){
    # Internal function to get status based on single job IDs.
    # Not in NAMESPACE. Used by retrieveText() to check status.
    #    
    # Args:
    #   jobID: Job ID to GET.
    #   api.key: API key for the HP IDOL API.
    #    
    # Returns:
    #   Job status (character)
    url <- "http://api.idolondemand.com/1/job/status/"
    get.url <- paste(url, jobId, "?apikey=", api.key, sep = "")
    # Get status and return
    status = as.character(content(GET(get.url))$status)
    return(status)
}
