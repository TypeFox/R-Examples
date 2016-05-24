dnb_get_url <- function(path, query, limit, start, token=dnb_token()) {
	req <- GET("http://services.dnb.de/", path=path, query=list(version="1.1", operation="searchRetrieve", accessToken=token, query=query, maximumRecords=limit, startRecord=start, recordSchema="MARC21-xml"))
	dnb_check(req)
	message("Request: ", req$url) # for debugging
	return(req)
}


dnb_check <- function(req) {
	if(req$status_code<400) return(invisible())
	message <- dnb_parse(req)$message
	stop("HTTP failure: ", req$status_code, "\n", message, call.=FALSE)
}


dnb_parse <- function(req) {
	xml <- content(req, as="text")
	if(identical(xml, "")) stop("Not output to parse", call.=FALSE)
	if(length(grep("text/xml", req$headers$'content-type', fixed=TRUE))==0) stop("No XML to parse", call.=FALSE)
	parsed <- as_list(read_xml(gsub("\n +", "", xml)))
	return(parsed)
}


dnb_token <- function(force=FALSE) {
	env <- Sys.getenv('DNB_TOKEN')
	if(!identical(env, "") && !force) return(env)
	if(!interactive()) {
		stop("Please set env var 'DNB_TOKEN' to your personal access token", call.=FALSE)
	}
	message("Couldn't find env var DNB_TOKEN.")
	message("Please enter your token and press enter")
	token <- readline(": ")
	if(identical(token, "")) {
		stop("Token entry failed", call.=FALSE)
	}
	message("Updating DNB_TOKEN env var to given token")
	Sys.setenv(DNB_TOKEN=token)
	return(token)
}


dnb_to_df <- function(lst) {
	# prepare data.frame
	nrec <- length(lst$records)
	df <- data.frame(matrix(nrow=nrec, ncol=17))
	names(df) <- c("id", "link", "author", "title", "subtitle", "publisher", "year", "language", "isbn", "price", "pages", "format", "edition", "keyword", "toc", "description", "cover")
	
	# get data
	for(r in 1:nrec) {
		rec <- lst$records[[r]]$recordData$record
		rec <- setNames(rec, sapply(rec, function(x) attributes(x)$tag))
		rec <- lapply(rec, function(x) lapply(x, function(y) setNames(y, attributes(y)$code)))
		rec <- lapply(rec, function(x) unlist(x, recursive=FALSE))
		
		if(!is.null(rec[["001"]])) {	# id/link
			df$id[r] <- rec[["001"]]
			df$link[r] <- paste0("http://d-nb.info/", rec[["001"]])	
		}
		if(!is.null(rec[["100"]][["a"]])) {	# author
			aut <- rec[["100"]][["a"]]
			if(!is.null(rec[["100"]][["4"]])) aut <- paste0(aut, " (", rec[["100"]][["4"]], ")")
			df$author[r] <- aut
		}
		if(length(which(names(rec)=="700"))>0) {	# co-author
			for(ca in which(names(rec)=="700")) {
				if(!is.null(rec[[ca]][["a"]])) {
					coaut <- rec[[ca]][["a"]]
					if(!is.null(rec[[ca]][["4"]])) coaut <- paste0(coaut, " (", rec[[ca]][["4"]], ")")
					if(is.na(df$author[r])) {
						df$author[r] <- coaut
					} else {
						df$author[r] <- paste(df$author[r], coaut, sep="; ")
					}
				}
			}
		}
		if(!is.null(rec[["245"]][["a"]])) {	# title
			df$title[r] <- rec[["245"]][["a"]]
		}
		if(!is.null(rec[["245"]][["b"]])) {	# subtitle
			df$subtitle[r] <- rec[["245"]][["b"]]
		}
		if(!is.null(rec[["264"]][["b"]])) {	# publisher
			pub <- rec[["264"]][["b"]]
			if(!is.null(rec[["264"]][["a"]])) pub <- paste0(pub, ", ", rec[["264"]][["a"]])
			df$publisher[r] <- pub
		}
		if(!is.null(rec[["264"]][["c"]])) {	# year
			df$year[r] <- rec[["264"]][["c"]]
		}
		if(!is.null(rec[["041"]][["a"]])) {	# language
			df$language[r] <- rec[["041"]][["a"]]
		}
		if(!is.null(rec[["024"]][["a"]])) {	# isbn
			df$isbn[r] <- rec[["024"]][["a"]]
		} else if(!is.null(rec[["020"]][["a"]])) {
			df$isbn[r] <- rec[["020"]][["a"]]
		}
		if(!is.null(rec[["020"]][["c"]])) {	# price
			df$price[r] <- rec[["020"]][["c"]]
		}
		if(!is.null(rec[["300"]][["a"]])) {	# pages
			df$pages[r] <- rec[["300"]][["a"]]
		}
		if(!is.null(rec[["300"]][["c"]])) {	# format
			df$format[r] <- rec[["300"]][["c"]]
		}
		if(!is.null(rec[["250"]][["a"]])) {	# edition
			df$edition[r] <- rec[["250"]][["a"]]
		}
		if(length(which(names(rec)=="689"))>0) {	# keyword
			for(kw in which(names(rec)=="689")) {
				if(!is.null(rec[[kw]][["a"]])) {
					if(is.na(df$keyword[r])) {
						df$keyword[r] <- rec[[kw]][["a"]]
					} else {
						df$keyword[r] <- paste(df$keyword[r], rec[[kw]][["a"]], sep="; ")
					}
				}
			}
		}
		if(length(which(names(rec)=="856"))>0) {	# toc/description
			for(kw in which(names(rec)=="856")) {
				if(!is.null(rec[[kw]][["3"]]) && !is.null(rec[[kw]][["u"]])) {
					if(rec[[kw]][["3"]]=="Inhaltsverzeichnis") {
						df$toc[r] <- rec[[kw]][["u"]]
					} else if(rec[[kw]][["3"]]=="Inhaltstext") {
						df$description[r] <- rec[[kw]][["u"]]
					}
				}
			}
		}
		if(!is.null(rec[["020"]][["9"]])) {	# cover
			df$cover[r] <- paste0("https://portal.dnb.de/opac/mvb/cover.htm?isbn=", rec[["020"]][["9"]])
		}	
	}
	
	return(df)
}
