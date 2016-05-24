
# Returns TRUE if object with given ident exists in Opasnet Base (WARNING! returns FALSE on any error!!!)
opbase.obj.exists <- function(ident, username = NULL, password = NULL)
{
	q <- list('ident' = ident)
	ret <- tryCatch(opbase.query(q, username, password), error = function(e) return(NULL))
	if (is.null(ret)){
		return(FALSE)
	} else {
		return(TRUE)
	}
}

# Returns locations of index
opbase.locations <- function(ident, index_name, series_id = NULL, username = NULL, password = NULL)
{
	query = list()
	query[['ident']] <- ident
	query[['index_name']] <- URLencode(index_name, reserved = TRUE)
	if (! is.null(series_id)) query[['series_id']] <- series_id
	ret <- opbase.query(query, username, password)
	return(ret$locations)
}

# Fetch all series ids of object
opbase.series <- function(ident, username = NULL, password = NULL, verbose = FALSE)
{
	query = list()
	query[['ident']] <- ident
	query[['username']] <- username
	query[['password']] <- password
	object <- opbase.query(query, username, password)

	if(is.null(object$acts)) stop("Acts not found!")
	
	ret = c()
	
	if (verbose) print(object$acts)
	
	for (a in object$acts)
	{
		ret <- append(ret, as.numeric(a$act$series_id))
	}
	
	return(unique(ret))
}

opbase.indices <- function(ident, act = NULL, username = NULL, password = NULL, verbose = FALSE)
{
	query = list()
	query[['ident']] <- ident
	query[['username']] <- username
	query[['password']] <- password
	object <- opbase.query(query, username, password)
	
	if(is.null(object$acts)) stop("Acts not found!")
	
	ret = c()
	
	if (verbose) print(object$acts)
	
	if (is.null(act)) act <- length(object$acts)
	
	for (a in object$acts[[act]]$indices)
	{
		ret <- c(ret, as.character(a$name))
	}
	
	return(ret)
}

# Read data from opasnet base 2
opbase.data <- function(ident, series_id = NULL, subset = NULL, verbose = FALSE, username = NULL, password = NULL, samples = NULL, exclude = NULL, include = NULL, range = NULL, optim_test = TRUE, ...) {
	
	query = list()

	query[['ident']] <- ident
	if (! is.null(subset)) query[['ident']] <- paste(query[['ident']], opbase.sanitize_subset_name(subset), sep='.')
	
	if (is.null(series_id))
	{
		query[['act']] <- 0
	}
	else
	{
		query[['series']] <- series_id
	}
	
	if (! is.null(samples)) query[['samples']] <- samples
	if (! is.null(exclude)) query[['exclude']] <- opbase.parse_locations(exclude, query[['ident']], series_id, username, password)
	if (! is.null(include)) query[['include']] <- opbase.parse_locations(include, query[['ident']], series_id, username, password)
	if (! is.null(range)) query[['range']] <- opbase.parse_range(range, query[['ident']], series_id, username, password)
	
	query[['username']] <- username
	query[['password']] <- password
	
	#if (verbose) print(query)
	
	# Run query to get KEY for downloading the actual data
	object <- opbase.query(query, username, password)
	
	#if (verbose) print(object)
	if (verbose) print("Object info and download key loaded!")
	
	if(is.null(object$key) || object$key == '') stop("Invalid download key retrieved!")
	
	out <- data.frame()
	data <- NULL
	first <- TRUE
	
	query = list()
	query[['key']] <- object$key
	
	while ((!is.null(data) && data != '') | first) {
		first <- FALSE

		if (verbose) print(paste('Loading data chunk from server... ',format(Sys.time(), "%H:%M:%OS3"),sep=''))
		temp <- opbase.query(query, username, password)
		data <- temp$data
		if (verbose) print('Data loaded ok!')
		if (verbose) print(paste('Processing data... ',format(Sys.time(), "%H:%M:%OS3"),sep=''))
		if (!is.null(data) && data != '')
		{
			temp <- fromJSON(data)
			if (optim_test){
				tmp <- list()
				iterator <- list()
				for (rivi in temp) {
					for (sarake in names(rivi)) {
						tmp[[sarake]] <- c(tmp[[sarake]], rivi[[sarake]])
						iterator[[sarake]] <- c(iterator[[sarake]], length(rivi[[sarake]]))
					}
				}
				# Using iterator it would be possible to implement multiple results on any column, but since this is not possible
				# at the moment only res will be checked.
				iterate <- FALSE
				if (is.null(samples) || (! is.null(samples) && samples > 0)) {
					if (prod(iterator[["res"]])>1) {
						for (sarake in names(tmp[names(tmp)!="res"])) {
							tmp[[sarake]] <- rep(tmp[[sarake]], times = iterator[["res"]])
						}
						iterate <- TRUE
					}
				}
				tmp <- data.frame(tmp)
				if (iterate) {
					iterations <- lapply(iterator[["res"]], f.iter)
					tmp$Iteration <- unlist(iterations)
				}
				# This method appears to be slower than the original with heavily iterated data (~15% difference with 5000 samples).
				# As the number of rows per chunk gets smaller, the difference between per-row and per-column approaches 
				# diminishes while this method wastes more resources calculating cell data lengths.
				out <- rbind(out, tmp)
			} else {
				if (verbose) print(paste('JSON parsed',format(Sys.time(), "%H:%M:%OS3"),sep=''))
				temp <- lapply(temp, data.frame)	# THIS FUNCTION IS RELATIVELY SLOW!!!! Could this be done in any other way?
				if (verbose) print(paste('Converted to data frame',format(Sys.time(), "%H:%M:%OS3"),sep=''))
				lengths <- lapply(temp, nrow)	
				if (verbose) print(paste('Row lengths resolved',format(Sys.time(), "%H:%M:%OS3"),sep=''))
				temp <- do.call("rbind", temp)
				if (verbose) print(paste('Rbind done',format(Sys.time(), "%H:%M:%OS3"),sep=''))
				
				if (is.null(samples) || (! is.null(samples) && samples > 0)) {
					if  (sum(unlist(lengths) > 1) > 0) {
						iterations <- lapply(lengths, f.iter)
						temp$Iteration <- unlist(iterations)
					}
				}	
				out <- rbind(out, temp)
				if (verbose) print(paste('Concatenated chunk to output',format(Sys.time(), "%H:%M:%OS3"),sep=''))
			}
		}
		if (verbose) print('Data processed ok!')
	}
	
	if (nrow(out) == 0) stop('Empty result set!')
	
	#if (verbose) print(out)


	if (is.null(samples) || samples > 0) {
		out <- out[,!colnames(out) %in% c("sid", "aid", "mean", "sd")]
		if ('res' %in% colnames(out)) {
			colnames(out)[colnames(out) == "res"] <- "Result"	
			a <- suppressWarnings(as.numeric(as.character(out$Result)))
			if (sum(is.na(a)) == 0) out$Result <- a
		}
	} else {
		out <- out[,!colnames(out) %in% c("sid", "aid")]
		colnames(out)[colnames(out) == "mean"] <- "Mean"	
		colnames(out)[colnames(out) == "sd"] <- "Sd"	
	}

	for(i in 1:length(object$indices)) {
		ind <- object$indices[[i]]
		temp <- as.character(ind$name)
		if (verbose) print(paste("Index ",i," name is ",temp,sep=''))
		colnames(out)[i] <- temp
	}
	
	
	
	return(out)
}

# Write data to the new opasnet database
opbase.upload <- function(
	input, 
	ident = NULL, 
	name = NULL, 
	subset = NULL, 
	obj_type = 'variable', 
	act_type = 'replace', 
	language = 'eng', 
	unit = '', 
	who = NULL, 
	rescol = NULL, 
	chunk_size = NULL, 
	verbose = FALSE, 
	username = NULL, 
	password = NULL, 
	index_units = NULL, 
	index_types = NULL
) {
	
	args <- opbase.parse_args()
	
	if (is.null(ident))
	{
		if (! is.null(args$wiki_page_id)){
			ident <- args$wiki_page_id
		} else {
			stop('Ident missing!')
		}
	}
	
	server <- 'cl1.opasnet.org'
	path <- '/opasnet_base_2/index.php'
	
	# Coerce input into a data frame if it isn't one already; get rid of empty cells
	if (is.array(input)) dataframe <- as.data.frame(as.table(input)) else dataframe <- input
	if (is.null(rescol)) {
		rescol <- colnames(dataframe) == "Freq"
		if (sum(rescol) == 1) rescol <- "Freq" else {
			rescol <- colnames(dataframe) == "Result"
			if (sum(rescol) == 1) rescol <- "Result" else {
				rescol <- colnames(dataframe) == "result"
				if (sum(rescol) == 1) rescol <- "result"
			}
		}}

	if (is.null(rescol)) stop('No result column could be defined!')
	if (verbose) print(paste('Rescol:',rescol, sep= ' '))
	
	dataframe <- dataframe[is.na(dataframe[,rescol]) == FALSE,]
	ColNames <- colnames(dataframe)[!(colnames(dataframe)%in%c(rescol, "id", "obs"))]
	
	ident <- tolower(ident)
	
	# Wiki id
	if (substr(ident, 1,5)=="op_en") {wiki_id <- 1; page <- substr(ident, 6, nchar(ident))} else {
		if (substr(ident, 1,5)=="op_fi") {wiki_id <- 2; page <- substr(ident, 6, nchar(ident))} else {
			if (substr(ident, 1,6)=="heande") {wiki_id <- 3; page <- substr(ident, 7, nchar(ident))} else {
				if (substr(ident, 1,4)=="test") {wiki_id <- 4; page <- substr(ident, 5, nchar(ident))} else {
					stop(paste("No wiki id found in ident ",ident,sep=''))}}}}
	
	n <- length(dataframe[1,rescol])
	page <- as.numeric(page)
	if (is.na(page)) stop("Could not convert characters following the wiki ident into a page number!\n")
	if (is.null(who)==TRUE) stop("uploader name not given")
	
	# If trying to append, check if object exists first (this functionality should be updated into opbase.obj.exists)
	if (act_type == 'append') {
		temp_id <- ident
		if (!is.null(subset)) {
			temp_id <- paste(temp_id, subset, sep = ".")
		}
		series <- tryCatch(
			opbase.series(temp_id),
			error = function(e) {
				if(grepl("^Query response error: Load by columns: record not found!", e[[1]])) { 
					return(NULL)
				} else {
					stop(paste("Unexpected error:", e[[1]]))
				}
			}
		)
		if (is.null(series)) {
			act_type <- 'replace'
		}
	}
	
	if (is.null(name)==TRUE && act_type == 'replace') stop("object name not given")
	
	# Build index list
	indices = list()

	for (i in 1:length(ColNames)) {
		if (is.null(index_types))
		{
			t = 'entity'; 
		} else {
			t = index_types[i]
		}
		if (is.null(index_units))
		{
			u = ''; 
		} else {
			u = index_units[i]
		}	
		indices[[i]] = list(type=t,name=ColNames[i],page=page,wiki_id=wiki_id,order_index=i,hidden=0,unit=u) 
	}
	
	header <- list(
			object = list(
					name = opbase.ensure_utf8(name),
					ident = ident,
					type = obj_type,
					page = page,
					wiki_id = wiki_id
			),
			act = list(
					unit = opbase.ensure_utf8(unit),
					who = opbase.ensure_utf8(who),
					samples = n,
					comments = "R upload",
					language = language
			),
			indices = indices
	)
	
	if (! is.null(subset))
	{
		header[['object']][['subset_name']] <- opbase.ensure_utf8(subset)
		header[['object']][['ident']] <- paste(ident, opbase.sanitize_subset_name(subset), sep='.')
	}
		
	if (act_type == 'replace')
	{
		method <- 'POST';
	}
	if (act_type == 'append')
	{
		method <- 'PUT';
	}
	
	# Do some authentication!!!
	if (is.null(username))
	{
		if (! is.null(args$user))
		{	
			header[['username']] <- args$user
			header[['password']] <- opbase.hashed_password(opbase.read_auth(args$user),  ident=header[['object']][['ident']])
		}
	}
	else
	{
		if (! is.null(password))
		{	
			header[['username']] <- username
			header[['password']] <- opbase.hashed_password(password,  ident=header[['object']][['ident']])
		}
	}	
	
	json <- toJSON(header)
	
	if (verbose) print(json)
	
	data <- list('_method' = method, 'json' = json)
	
	response <- postToHost(server, path, data)
	
	if (is.null(response)) stop('Server is not responding!!!')

	if (verbose) print(response)
	
	# Parse JSON data from the server response
	response <- fromJSON(regmatches(response, regexpr('\\{.+\\}',response)))
	if (! is.null(response$error)) stop(response$error)
	if (is.null(response$key) || response$key == '') stop("Invalid upload key retrieved!")
	
	total_rows <- nrow(dataframe)
	
	# Automatic chunksize?
	if (is.null(chunk_size))
	{
		chunk_size <- 1000
	
		if (n > 1)
		{
			div <-  n/ 30
			if (div < 1) div <- 1
			chunk_size = round(chunk_size / div)	
			if (chunk_size < 1) chunk_size <- 1
		}
	}
	
	if (chunk_size > total_rows) chunk_size <- total_rows
	
	start <- 1
	end <- chunk_size
	
	rows <- 0
	
	raw_data <- list('key' = response$key, 'indices' =  indices)

	# Do some authentication!!!
	if (is.null(username))
	{
		if (! is.null(args$user))
		{	
			raw_data[['username']] <- args$user
			raw_data[['password']] <- opbase.hashed_password(opbase.read_auth(args$user),  key=response$key)
		}
	}
	else
	{
		if (! is.null(password))
		{	
			raw_data[['username']] <- username
			raw_data[['password']] <- opbase.hashed_password(password,  key=response$key)
		}
	}	
	
	
	# Write the data
	repeat
	{
		data_chunk <- data.frame(lapply(dataframe[start:end,], as.character), stringsAsFactors=FALSE, check.names=FALSE)
		chunk_rows <- nrow(data_chunk)
		
		if (verbose) print(data_chunk)
		
		data_rows = list()
		# Create data list for JSON
		for (r in 1:chunk_rows)
		{
			row <- data_chunk[r,]
			
			v <- as.character(row[rescol])
			record <- list('res' = v)
			
			for (i in 1:length(ColNames))
			{
				v <- as.character(row[ColNames[i]])
				record[[toString(i)]] <- v
			}
			
			if (verbose) print(record)

			data_rows[[r]] <- record
		}
		
		#if (verbose) print(data_rows)
	
		raw_data[['data']] <- data_rows
		json <- toJSON(raw_data)
		
		if (verbose) print(json)
		
		data <- list('json' = json)
		response <- postToHost(server, path, data)
		
		if (is.null(response)) stop('Server is not responding!!!')
		
		if (verbose) print(response)
		
		# Parse JSON data from the server response
		response <- fromJSON(regmatches(response, regexpr('\\{.+\\}',response)))
		if (! is.null(response$error)) stop(response$error)
		
		if (verbose) print(response)
		
		rr <- as.integer(response$rows)
		
		if (rr != chunk_rows) stop(paste('Invalid inserted rows count! ', rr, ' vs ', chunk_rows, sep=''))
	
		rows <- (rows + rr)
		
		if (end >= total_rows) break
		
		start <- (end + 1)
		end <- (end + chunk_size)
		if (end > total_rows) end <- total_rows
	}
	
	return(rows)
	
}

# Private function to parse include or exclude list and return vector containing corresponding index idents and location ids
opbase.parse_locations <- function(locs, ident, series_id = NULL, username = NULL, password = NULL) {
	ret = c()
	
	#print(locs)
	
	for(i in names(locs))
	{
		query = list()
		query[['ident']] <- ident
		query[['index_name']] <- URLencode(i, reserved = TRUE)
		if (! is.null(series_id)) query[['series']] <- series_id
		object <- opbase.query(query, username, password)
		#ind <- object$index
		loc_ids = c()
		for (loc in locs[[i]])
		{
			# Seek thru all locations of an index
			for (lid in names(object$locations))
			{
				if (object$locations[[lid]] == loc)
				{
					loc_ids <- c(loc_ids, lid) 
					break
				}
			}
		}
		ret <- c(ret, paste(object$index$ident,paste(loc_ids,collapse=','),sep=','))
	}
	#print(ret)
	return(ret)
}

# Private function to parse range list and return vector containing corresponding index idents and values
#
# range format: list(<index name> = c(<min>,<max>), ...)

opbase.parse_range <- function(range, ident, series_id = NULL, username = NULL, password = NULL) {
	ret = c()
	
	#print(locs)
	
	for(i in names(range))
	{
		query = list()
		query[['ident']] <- ident
		query[['index_name']] <- URLencode(i, reserved = TRUE)
		if (! is.null(series_id)) query[['series']] <- series_id
		object <- opbase.query(query, username, password)
		if (length(range[[i]]) != 2) stop('Invalid range array!')
		if (! is.na(range[[i]][[1]])) {
			min = range[[i]][[1]] } else {
			min = ''
		}
		if (! is.na(range[[i]][[2]])) {
			max = range[[i]][[2]] } else {
			max = ''
		}
		ret <- c(ret, paste(object$index$ident,min,max,sep=';'))
	}
	#print(ret)
	return(ret)
}

# Private function to get authentication
opbase.read_auth <- function(user) {
	auth <- fromJSON(paste(readLines("/var/www/html/rtools_server/offline/opasnet.json"),  collapse = ""))
	return(auth[[user]])
}

# Private function to create hashed password
opbase.hashed_password <- function(password, index = NULL, ident = NULL, key = NULL) {
	str <- ''
	if (! is.null(index)) str <- paste(str, index, sep='')
	if (! is.null(ident)) str <- paste(str, ident, sep='')
	if (! is.null(key)) str <- paste(str, key, sep='')
	str <- paste(str, password, sep='')
	return(digest(str, algo="md5", serialize=FALSE))
}

# Private function to make queries to server
opbase.query <- function(data, username = NULL, password = NULL) {

	args <- opbase.parse_args() 
	
	index <- NULL
	ident <- NULL
	key <- NULL
	
	if (! is.null(data[['index']])) index <- data[['index']]
	if (! is.null(data[['ident']])) ident <- data[['ident']]
	if (! is.null(data[['key']])) key <- data[['key']]
	
	# Do some authentication!!!
	if (is.null(username))
	{
		if (! is.null(args$user))
		{	
			data[['username']] <- args$user
			data[['password']] <- opbase.hashed_password(opbase.read_auth(args$user), index = index, ident = ident, key = key)
		}
	}
	else
	{
		if (! is.null(password))
		{	
			data[['username']] <- username
			data[['password']] <- opbase.hashed_password(password, index = index, ident = ident, key = key)
		}
	}
	
	# Build http-query key / value pairs
	tmp = c(1:length(data))
	i <-1
	for (k in names(data)){
		if (length(data[[k]]) > 1)
		{
			sepi <- paste(k,'[]=',sep='')
			tmp[i] <- paste(sepi, paste(data[[k]], collapse=paste("&", sepi, sep='')), sep='')
		} else {
			tmp[i] <- paste(k, '=', data[[k]], sep= '')
		}
		i <- i + 1
	}
	
	#print( paste(tmp, collapse='&') )
		
	url <- paste("http://cl1.opasnet.org/opasnet_base_2/index.php?", paste(tmp, collapse='&'), sep = "")
	
	response <- fromJSON(
			paste(
					readLines(url),  
					collapse = ""
			)
	)
	
	if (is.null(response))
	{
		stop("Opasnet server is not responding! Unable to query!")
	}
	
	if (! is.null(response$error))
	{
		stop(paste("Query response error: ",response$error, sep= ''))
	}
	
	return(response)
}

f.iter <- function(x) {
	1:x
}

#opbase.objects <- fromJSON(
#	paste(
#		readLines(
#			"http://cl1.opasnet.org/opasnet_base_2/index.php", 
#		), 
#		collapse = ""
#	)
#)

# Private function to parse arguments
opbase.parse_args <- function()
{
	# Parse arguments
	targs <- strsplit(commandArgs(trailingOnly = TRUE),",")
	args = list()
	
	if (length(targs) > 0)
		for(i in targs[[1]])
		{
			tmp = strsplit(i,"=")
			key <- tmp[[1]][1]
			value <- tmp[[1]][2]
			args[[key]] <- value
		}
	return(args)
}

# Private function to sanitize object subset data name 
opbase.sanitize_subset_name <- function(name)
{
	enc <- Encoding(name)
	# Make lowercase
	name <- tolower(name)
	# Remove all punctuation marks: ! " # $ % & ' ( ) * + , - . / : ; < = > ? @ [ \ ] ^ _ ` { | } ~.
	name <- gsub('[[:punct:] ]','_',name)
	# Convert to ASCII?
	if (enc != 'unknown') name <- iconv(name, enc, "ASCII//TRANSLIT","_")
	# Truncate multiple underscores
	name <- gsub('_+','_',name)
	# Remove trailing or leading underscores
	name <- gsub('^_|_$','',name)
	return(name)
}

opbase.ensure_utf8 <- function(str)
{
	enc <- Encoding(str)
	if (enc == 'UTF8' || enc == 'unknown') return(str)
	return(iconv(str, enc, "UTF8","?"))
}
