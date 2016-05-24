getRequestResults <- function(job.id, api.key){
    # Internal function to get results based on single job IDs.
    # Not in NAMESPACE. Used by retrieveText() to get transcription.
    #    
    # Args:
    #   job.id: Job ID to GET.
    #   api.key: API key for the HP IDOL API.
    #    
    # Returns:
    #   Transcribed text (character)
    # Check status, if queued, return status
    status <- getStatus(job.id, api.key)
    if(status != 'finished') {
        return(status)
    }
    # If finished, snag that text!
    url <- "http://api.idolondemand.com/1/job/result/"
    get.url <- paste(url, job.id, "?apikey=", api.key, sep = "")
    call.results = content(GET(get.url))
    # Get the text and clean it of tags (<NOISE/MUSIC>)
    out.text <- call.results$actions[[1]]$result$document[[1]]
    out.text <- cleanFun(out.text)
    return(out.text)
}
