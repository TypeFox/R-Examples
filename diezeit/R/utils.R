zeit_get_url <- function(path, ..., key=zeit_key()) {
	#auth <- zeit_auth(key)
	#req <- GET("http://api.zeit.de/", path=path, auth, query=list(...))
	req <- GET("http://api.zeit.de/", path=path, query=list(api_key = key, ...))
	zeit_check(req)
	#message("Request: ", req$url) # for debugging
	return(req)
}


#zeit_auth <- function(key=zeit_key()) {
#	authenticate(key, "")
#}


zeit_check <- function(req) {
	if(req$status_code < 400) return(invisible())
	message <- zeit_parse(req)$message
	stop("HTTP failure: ", req$status_code, "\n", message, call.=FALSE)
}


zeit_parse <- function(req) {
	json <- content(req, as="text")
	if(identical(json, "")) stop("Not output to parse", call.=FALSE)
	if(length(grep("application/json", req$headers$'content-type', fixed=TRUE)) == 0) stop("No JSON to parse", call.=FALSE)
	fromJSON(json, simplifyVector=FALSE)
}


zeit_key <- function(force=FALSE) {
	env <- Sys.getenv('ZEIT_KEY')
	if(!identical(env, "") && !force) return(env)

	if(!interactive()) {
		stop("Please set env var 'ZEIT_KEY' to your personal ZEIT api key", call.=FALSE)
	}

	message("Couldn't find env var ZEIT_KEY.")
	message("Please enter your key and press enter")
	key <- readline(": ")

	if(identical(key, "")) {
		stop("Key entry failed", call.=FALSE)
	}

	message("Updating ZEIT_KEY env var to given key")
	Sys.setenv(ZEIT_KEY=key)

	return(key)
}


hasKey <- function() {
	!identical(zeit_key(), "")
}
