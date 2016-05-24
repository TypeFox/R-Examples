.ZendeskEnv <- new.env()
.ZendeskEnv$data <- list()

.onLoad <- function(libname, pkgname){
    if(is.null(.ZendeskEnv$data) == FALSE){
    	.ZendeskEnv$data <- list(
        	username <- NULL,
        	password <- NULL,
        	url <- NULL,
		      users = "/api/v2/users.json",
		      tickets = "/api/v2/tickets.json",
		      audits = "/audits.json",
		      organizations = "/api/v2/organizations.json",
          ticket_metrics = "/api/v2/ticket_metrics.json",
          satisfaction_ratings = "/api/v2/satisfaction_ratings.json"
		
            )
    }
}

unlistDataFrame <- 
function(dataframe){
	n <- dim(dataframe)[1]
	dim.df <- dim(dataframe)[2]

##	Not sure why vectorized version of this doesn't work
##	apply(dataframe, 2, function(x) { if (length(unlist(x)) == n) unlist(x)} )

	for(i in 1:dim.df){
		if(length(unlist(dataframe[,i])) == n){
			dataframe[,i] <- unlist(dataframe[,i])
		} 
		if(length(unlist(dataframe[,i])) < n){
			dataframe[,i] <- as.character(dataframe[,i])
		}
	}

	return (dataframe)
}

      
