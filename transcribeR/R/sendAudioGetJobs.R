sendAudioGetJobs <- function(wav.dir, api.key, interval = "-1",
                             encode = "multipart", existing.csv = NULL, csv.location,
                             language = "en-US", verbose = FALSE){
    # Main user function to POST to HP IDOL Speech Recognition
    # API and write jobs to job.file (a filename)
    #    
    # Args:
    #   wav.dir: Directory that contains wav files. Will
    #            only try to upload wav files in dir.
    #   api.key: API key for the HP IDOL API.
    #   interval: HP API arg.
    #   encode: HP API arg.
    #   job.file: CSV file where jobs.csv will be written.
    #   language: HP API arg.
    #    
    # Returns:
    #   Message indicating success, automatically writes jobs
    #   csv to file.

    error.messages <- NULL
    url <- "https://api.idolondemand.com/1/api/async/recognizespeech/v1"
    # get all files in wav directory
    wav.dir <- gsub('/?$', '/', wav.dir) # add trailing '/' if missing
    wav.filenames <- Sys.glob(c(paste(wav.dir,'*.wav', sep = ''),paste(wav.dir,'*.mp4', sep = ''),
                                paste(wav.dir,'*.mp3', sep = ''),paste(wav.dir,'*.wma', sep = '')))
    total.number.of.files <- length(wav.filenames)
    # holder for content
    out.list <- list()
    ex.v = c(1:6)
    ex.v[1] <- "DATE"
    ex.v[2] <- "APIKEY"
    ex.v[3] <- "FILENAME"
    ex.v[4] <- "LANGUAGE"
    ex.v[5] <- "JOBID"
    ex.v[6] <- "TRANSCRIPT"
    
    is.file.created <- createJobCSV(existing.csv, csv.location) # Boolean, TRUE if a file is created
    print(paste("Is a new file created?",is.file.created))
    i <- 0 # used for counting uploaded files
    x <- 0 # used for verbose mode
    if(is.file.created == FALSE){
        existing.job.csv <- read.csv(existing.csv)
        if(any(colnames(existing.job.csv) != ex.v)){ # Check if the provided file is correctly formatted
            error.messages <- "incorrect csv type"
            stop("This doesn't appear to be a transcribeR jobs.csv, please provide another filename")
        }
        for(fpath in wav.filenames){
            fn <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", fpath)
            if(!(fn %in% existing.job.csv$FILENAME)){
                attempt <- POST(
                    url,
                    body = list(
                        file = upload_file(fpath),
                        apikey = api.key,
                        language = language, 
                        interval = interval
                        ),
                    encode = encode)
                stop_for_status(attempt)
                out.list[[fn]] <- content(attempt)
                i = i + 1
            } #
            # print out the status of the upload and current filename
            if (verbose) {
                x = x + 1
                print(paste("Current Upload: ", fn, sep = ""))
                print(paste("Current % of all audio uploaded: ", round(x / total.number.of.files * 100, 4), "%", sep = ""))
            }
        }
    }
    else { # a new file was created
        for(fpath in wav.filenames){
            fn <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", fpath)
            attempt <- POST(
                url,
                body = list(
                    file = upload_file(fpath),
                    apikey = api.key,
                    language = language,
                    interval = interval
                    ),
                encode = encode)
            stop_for_status(attempt)
            out.list[[fn]] <- content(attempt)
            i = i + 1
            if (verbose) {
                x = x + 1
                print(paste("Current Upload: ", fn, sep = ""))
                print(paste("Current % of all audio uploaded: ", round(x / total.number.of.files * 100, 4), "%", sep = ""))
            }
        }
    }
    DATE <- rep(as.character(Sys.Date()),length(out.list))
    APIKEY <- rep(api.key,length(out.list))
    FILENAME <- names(out.list)
    LANGUAGE <- rep(language, length(out.list))
    JOBID <- unname(unlist(lapply(out.list, function(x) x[['jobID']])))
    TRANSCRIPT <- rep("",length(out.list))
    
    df <- data.frame(DATE, APIKEY, FILENAME, LANGUAGE, JOBID, TRANSCRIPT, stringsAsFactors=FALSE)
    row.names(df) <- NULL
    appendToCSV(csv.location, df, append = TRUE, sep=",", row.names=FALSE, col.names=FALSE)
   
    if(is.null(error.messages)){
       print(paste("Jobs successfully uploaded, a transcribeR CSV was written to", csv.location))
       print(paste("The number of files uploaded was", i))
    }
    else {
       print("Error in uploading jobs and/or collecting job IDs.")
       print(error.messages)
    }
}
