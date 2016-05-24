tibbr.subject <-
function(connection, subject = NULL) {
	# Sanity check arguments
	if(!is(connection, "tibbrConnection")) {
		stop("connection must be a tibbrConnection")
	}
	subjvalid <- (is.character(subject) && length(subject) == 1 && !is.na(subject)) || is.null(subject)
	if(!subjvalid) {
		stop("subject must be a one-element character vector")
	}

	# Initial empty list
	results <- list()

	# Create the get request query structures
	query.names <- c("params[page]", "params[per_page]")
	query.vals <- c("1", "25")
	if(!is.null(subject)) {
		query.names <- c(query.names, "params[search_str]")
		query.vals <- c(query.vals, subject)
	}

	# Do the first request
	res.obj <- getSubjectPage(connection, query.names, query.vals)
	results <- append(results, res.obj$items)

	# Loop over the remaining pages, appending their results
	if(res.obj$total_pages > 1) {
		for(page in as.character(seq(2, res.obj$total_pages))) {
			query.vals[1] <- page
			res.obj <- getSubjectPage(connection, query.names, query.vals)
			results <- append(results, res.obj$items)
		}
	}

	# Correctly class and return the final results
	lapply(results, function(x) {
		class(x) <- c("tibbrSubject", "tibbrPostable")
		class(x$user) <- c("tibbrUser", "tibbrPostable")
		x
	})
}

getSubjectPage <-
function(connection, query.names, query.vals) {
	# Create the get request query string
	query <- paste(paste(sapply(query.names, URLencode), sapply(query.vals, URLencode), sep="="), collapse="&")

	# Perform the HTTP request
	res <- httpGet(paste("https://", connection$server, "/a/users/", connection$user$id, "/search_subjects.json?", query, sep=""),
	               c(`Client_key` = connection$clientkey,
			 `Auth_token` = connection$user$auth_token))

	# Analyze the result body
	res.obj <- try(json.parse(res), silent=TRUE)
	if(is(res.obj, "try-error")) {
		# json.parse could not parse the response as JSON
		stop("malformed response recieved from tibbr")
	}
	if(!is.null(res.obj$error)) {
		# Error reported by tibbr
		# {"error":"The authentication process failed."}
		stop(res.obj$error)
	} else {
		# Successful search; tibbr returns a collection of subjects
		res.obj
	}
}

