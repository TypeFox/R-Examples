httpLog <-
function(method, url, headers, reqbody, respbody) {
	if(!is.null(dbg <- getOption("tibbrConnector.debug"))) {
		f <- file(dbg, "ab")
		cat(method, "to", url, "\n", file=f)
		cat("REQUEST HEADERS:\n", file=f)
		cat(paste(names(headers), headers, sep=": ", collapse="\n"), file=f)
		if(!is.null(reqbody)) {
			cat("\nREQUEST BODY:\n", file=f)
			writeBin(reqbody, f)
			cat("\nEND REQUEST BODY", file=f)
		}
		cat("\nRESPONSE BODY:\n", file=f)
		writeBin(respbody, f)
		cat("\nEND RESPONSE BODY\n--------------------------------------------------------------------------\n", file=f)
		close(f)
	}
}

if(version$language == "R") {
## BEGIN R SPECIFIC IMPLEMENTATION

httpGet <-
function(url, headers) {
	reader <- basicTextGatherer()
	curlPerform(url = url, httpheader = headers, writefunction=reader$update)
	b <- reader$value()
	httpLog("GET", url, headers, NULL, b)
	b
}

httpPost <-
function(url, headers, body) {
	reader <- basicTextGatherer()
	curlPerform(url = url, httpheader = headers, postfields = body, writefunction=reader$update)
	b <- reader$value()
	httpLog("POST", url, headers, body, b)
	b
}

httpPut <-
function(url, headers, body) {
	bytes <- charToRaw(body)
	reader <- basicTextGatherer()
	curlPerform(url = url, httpheader = headers, infilesize = length(bytes), readfunction = bytes, writefunction = reader$update, customrequest = "PUT", upload = TRUE)
	b <- reader$value()
	httpLog("PUT", url, headers, body, b)
	b
}

## END R SPECIFIC IMPLEMENTATION
} else {
## BEGIN TERR SPECIFIC IMPLEMENTATION

httpGet <-
function(url, headers) {
	tryCatch(b <- http.get(url, headers)$body, error=function(e) {
		httpLog("GET", url, headers, NULL, paste("ERROR:", geterrmessage()))
		stop(e)
	})
	httpLog("GET", url, headers, NULL, b)
	b
}

httpPost <-
function(url, headers, body) {
	tryCatch(b <- http.post(url, body, headers)$body, error=function(e) {
		httpLog("POST", url, headers, body, paste("ERROR:", geterrmessage()))
		stop(e)
	})
	httpLog("POST", url, headers, body, b)
	b
}

httpPut <-
function(url, headers, body) {
	bytes <- charToRaw(body)
	tryCatch(b <- http.put(url, body, headers)$body, error=function(e) {
		httpLog("PUT", url, headers, body, paste("ERROR:", geterrmessage()))
		stop(e)
	})
	httpLog("PUT", url, headers, body, b)
	b
}

## END TERR SPECIFIC IMPLEMENTATION
}

FormEncode <-
function(x) {
	gsub("%20", "+", URLencode(x))
}
