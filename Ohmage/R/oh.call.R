oh.call <- function(xpath, serverurl=getOption("SERVERURL"), token=getOption("TOKEN"), 
	responseformat="json", style="post", verbose=FALSE, recycle=FALSE, ...){
	
	if(is.null(serverurl)){
		stop("Serverurl is missing. Either use oh.login() or specify 'serverurl' and 'token' arguments in plot call.")
	}
	
	posturl <- paste(serverurl, xpath, sep="");
	
	if(recycle){
		curl = getOption("OhCurlHandle");
		h = getOption("OhCurlReader");
		h$reset();
		increment.curlcount();
	} else {
		curl = getCurlHandle()
		h = dynCurlReader(curl, binary = TRUE)		
	}

	#some calls don't need a token (e.g. /auth_token)
	if(!is.null(token)){
		HTTPPARAMS <- list(client=getOption("CLIENTNAME"), auth_token=token, ...);	  		
	} else {
		if(verbose) message("Note: calling without auth_token...")		
		HTTPPARAMS <- list(client=getOption("CLIENTNAME"), ...);		
	}
	
	if(verbose){
		cat(" -------- Ohmage verbose output (", as.character(Sys.time(), usetz = TRUE), ") --------\n\n")
		cat("POST ", posturl, "\n\n")
		print.param.list(HTTPPARAMS);
	}

	#actual HTTP POST
	postForm(curl = curl, uri=posturl, style=style, binary=TRUE, .params=HTTPPARAMS,
		.opts = list(sslversion=3, ssl.verifyhost= FALSE, ssl.verifypeer=FALSE, headerfunction = h$update, verbose = verbose, connecttimeout=10));	
	
	
	#parse response
	headers <- parseHTTPHeader(h$header());
	httpstatus <- headers[["status"]];
	response <- h$value();
	
	if(verbose) {
		cat("-------- Ohmage response: --------\n\n")
		cat("HTTP", httpstatus, "\n");
	}

	if(httpstatus != 200){
		if(is.raw(response)){
			stop("Ohmage error: HTTP ", httpstatus, ".\n", rawToChar(response), "\n");	
		} else {
			stop("Tomcat error: HTTP ", httpstatus, ".\n", response, "\n");
		}
	}
	
	#this should never happen:
	if(length(response)==0){
		stop("server returned no content (check tomcat error log).");
	}
	
	if(!is.raw(response)){
		stop(response);
	}

	if(verbose){
		cat(rawToChar(response),"\n\n");
	}	
	
	if(responseformat == "file"){

		if(length(response) < 1000){
			warning("File less than 1KB: ", rawToChar(response));
		}

		tf <- tempfile(pattern="image", tmpdir="/tmp");
		if(file.create(tf)){
			con <- file(tf,"wb");
		} else {
			stop("Could not create temporary file.")
		}
		writeBin(as.vector(response), con);
		close(con);

		attr(tf,"Content-Type") <- attr(response,"Content-Type");
		return(tf);		
	}
	
	xhr <- fromJSON(rawToChar(response), simplifyWithNames=FALSE);
	
	if(xhr$result == "success"){
		return(xhr)
	} else {
		if(!is.null(xhr$errors[[1]]$text) && xhr$errors[[1]]$text == "The token is unknown."){
			oh.logout();
			stop("The token is unknown. Your session might have timed out. Please re-login.")
		}
		stop(xhr$errors[[1]]$text);
	}
}
