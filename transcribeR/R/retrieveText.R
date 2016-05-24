retrieveText <- function(job.file, api.key) {
    # Function to retrieve asynchronous 
    #    
    # Args:
    #   job.file: job.file created from sendAudio()
    #   api.key: api key for HP IDOL API
    #    
    # Returns:
    #    
    jobs <- read.csv(job.file)
    # If jobs csv has already been through retrieveText,
    # it will have a 'transcribed.text' column. Only get the
    # jobs that were in queue when last called.
    untranscribed.inds <- which(jobs$TRANSCRIPT == 'queued' | jobs$TRANSCRIPT == '' | is.na(jobs$TRANSCRIPT) | is.null(jobs$TRANSCRIPT))
    
    # try to transcribe all job.IDs 
    for(ind in untranscribed.inds){
        ID <- jobs$JOBID[ind]
        text <- getRequestResults(ID, api.key = api.key)        
        jobs$TRANSCRIPT[ind] <- text
    }
    write.csv(jobs, file = job.file, row.names = FALSE)
    return(jobs)
}
