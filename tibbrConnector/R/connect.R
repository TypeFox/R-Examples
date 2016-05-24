tibbr.connect <-
function(server, user, password) {
	# Sanity check arguments
	if(missing(server)) {
		stop("must specify tibbr server to connect to")
	}
	if(regexpr("^https?:", server) != -1) {
		stop("must not specify 'http://' or 'https://' prefix")
	}
	if(missing(user)) {
		stop("must specify user to connect to tibbr as")
	}
	if(missing(password)) {
		stop("must specify password for user")
	}

	# Create the authentication request body
	req.names <- c("params[login]", "params[password]")
	req.vals <- c(user, password)
	req <- paste(paste(sapply(req.names, FormEncode), sapply(req.vals, FormEncode), sep="="), collapse="&")

	# Create a new client key
	baseuuid <- paste(sample(c(letters[1:6],0:9), 30, replace=TRUE), collapse="")
	client.key <- paste(substr(baseuuid,1,8),
			    "-",
			    substr(baseuuid,9,12),
			    "-",
			    "4",
			    substr(baseuuid,13,15),
			    "-",
			    sample(c("8","9","a","b"),1),
			    substr(baseuuid,16,18),
			    "-",
			    substr(baseuuid,19,30),
			    sep="", collapse="")

	# Perform the HTTP request
	res <- httpPost(paste("https://", server, "/a/users/login.json", sep=""),
			c(`Content-Type` = "application/x-www-form-urlencoded", 
			  `Client_key` = client.key),
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
	} else {
		# Successful login; tibbr returns logged in user
		# Return the connection object
		class(res.obj) <- c("tibbrUser", "tibbrPostable")
		structure(list(server = server, user = res.obj, clientkey = client.key), class = "tibbrConnection")
	}
}

tibbr.disconnect <-
function(connection) {
	# Sanity check arguments
	if(!is(connection, "tibbrConnection"))
		stop("connection must be a tibbrConnection")

	# Perform the HTTP request
	res <- httpPut(paste("https://", connection$server, "/a/users/", connection$user$id, "/logout.json", sep=""),
		       c(`Content-Type` = "application/x-www-form-urlencoded",
			 `Client_key` = connection$clientkey,
			 `Auth_token` = connection$user$auth_token),
		       "")

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
		# Successful logout
		invisible(TRUE)
	}
}

