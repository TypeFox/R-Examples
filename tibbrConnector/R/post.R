tibbr.post <-
function(connection, content, subject = NULL, links = character(0), attachments = character(0)) {
	# Sanity check arguments
	if(!is(connection, "tibbrConnection")) {
		stop("connection must be a tibbrConnection")
	}
	if(missing(content)) {
		stop("must specify message content to post")
	}
	if(!is.character(links)) {
		stop("links argument must be a character vector")
	}
	if(!is.character(attachments)) {
		stop("attachments argument must be a character vector")
	}

	# Create the post request body
	boundary <- paste(sample(LETTERS, 32, replace=TRUE), collapse="")
	## Message content
	req <- c(charToRaw("--"), charToRaw(boundary), charToRaw("\r\n"))
	req <- c(req, charToRaw("Content-Disposition: form-data; name=\"message[content]\"\r\n\r\n"))
	subj.str <- if(missing(subject)) "" else {
		subs <- if(is(subject, "tibbrPostable")) postingName(subject) else sapply(subject, postingName)
		paste(subs, collapse=" ")
	}
	req <- c(req, charToRaw(paste(subj.str, content)))
	## Links
	for(i in seq_along(links)) {
		req <- c(req, charToRaw("\r\n--"), charToRaw(boundary), charToRaw("\r\n"))
		req <- c(req, charToRaw(paste("Content-Disposition: form-data; name=\"message[links[[", i - 1, "][url]]]\"\r\n\r\n", sep="")))
		req <- c(req, charToRaw(links[i]))
	}
	## Attachments
	for(i in seq_along(attachments)) {
		if(!file.exists(attachments[i])) {
			stop("attachment \"", attachments[i], "\" does not exist")
		}
		req <- c(req, charToRaw("\r\n--"), charToRaw(boundary), charToRaw("\r\n"))
		req <- c(req, charToRaw(paste("Content-Disposition: form-data; name=\"message[assets[[", i - 1, "][data]]]\"; filename=\"", basename(attachments[i]), "\"\r\n", sep="")))
		req <- c(req, charToRaw(paste("Content-Type: ", guessMIMEType(attachments[i]), "\r\n\r\n", sep="")))
		req <- c(req, readBin(attachments[i], raw(), file.info(attachments[i])$size))
	}
	req <- c(req, charToRaw("\r\n--"), charToRaw(boundary), charToRaw("--\r\n"))

	# Perform the HTTP request
	res <- httpPost(paste("https://", connection$server, "/a/messages.json", sep=""),
			c(`Content-Type` = paste("multipart/form-data; boundary=", boundary, sep=""),
			  `Client_key` = connection$clientkey,
			  `Auth_token` = connection$user$auth_token,
			  `Content-Length` = length(req)
			 ),
			req)

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
	} else if(is(res.obj, "JsonArray")
		  && length(res.obj) == 1
		  && length(res.obj[[1]]) == 2
		  && res.obj[[1]][[1]] == "base") {
		# Error reported by tibbr
		# [["base","tibbr cannot find the subject BadName."]]
		stop(res.obj[[1]][[2]])
	} else {
		# Successful posting
		invisible(TRUE)
	}
}

