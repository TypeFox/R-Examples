getTicket <- function(ticket.id){
        curl = getCurlHandle()
        stopPaging <- FALSE
        result <- list()
        i <- 1

	## Need to page through the results since only 100 are returned at a time
	while(stopPaging == FALSE){
            result[[i]]<-getURL(paste(.ZendeskEnv$data$url, gsub(".json", "", .ZendeskEnv$data$tickets), "/", ticket.id, ".json?page=" ,i, sep=""), curl=curl, 
                          ssl.verifypeer=FALSE, .opts=list(userpwd=(paste(.ZendeskEnv$data$username, .ZendeskEnv$data$password, sep=":"))))
            if(is.null(fromJSON(result[[i]])$next_page)){
                stopPaging <- TRUE
            }
            i <- i + 1
        }

	## Transform the JSON data to a data.frame
        json.data <- lapply(unlist(result), fromJSON)
        pre.result <- lapply(json.data, function(x) do.call("rbind", x))
        final.result<-do.call("rbind", pre.result)
        ticket.df <- data.frame(final.result)
	ticket.df <- unlistDataFrame(ticket.df)
        return(ticket.df)
}



