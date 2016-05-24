# Encode single R object with given key using ECB encryption, returns ciphertext as raw vector

objects.encode <- function(obj, key){
	
	# Check key
	key <- charToRaw(key)
	if (length(key) %% 16 != 0) stop('Invalid key length! Must be 16, 32 or 64 ASCII chars!')
	
	sobj <- serialize(obj, NULL, TRUE)
	m <- length(sobj) %% 16
	add <- 0

	if (m != 0){
		add <- 16 - m
		sobj <- c(sobj, rep(as.raw(00),each=add))
	}
	
	if (length(sobj) %% 16 != 0) stop('Something went wrong adding extra chars!')
	
	aes <- AES(key, mode="ECB")
	eobj <- aes$encrypt(sobj)
	
# Add character (as HEX) to tell us the amount of added zeros
	eobj <- c(as.raw(sprintf("%x", add)), eobj)
	return(eobj)
}

# Decode raw ciphertext encoded with objects.encode, use the same key obviously!

objects.decode <- function(eobj, key){

	# Check key
	key <- charToRaw(key)
	if (length(key) %% 16 != 0) stop('Invalid key length! Must be 16, 32 or 64 ASCII chars!')
	
	xtra <- as.integer(eobj[[1]])
	aes <- AES(key, mode="ECB")
	sobj <- aes$decrypt(eobj[-1], raw=TRUE)
	
	# Remove extra chars?
	if (xtra > 0) {
		l <- length(sobj)
		tmp <- c((l-xtra+1):l)
		sobj <- sobj[-(tmp)]
	}
	
	obj <- unserialize(sobj)
	return(obj)
}

# Wrapper for save-method, writes desired objects to run path as rdata

objects.put <- function(..., list = character()){
	
	# Parse arguments
	targs <- strsplit(commandArgs(trailingOnly = TRUE),",")
	args = list()
	
	if (length(targs) == 0) stop('This function can be used within Opasnet only!!!')
	
	for(i in targs[[1]])
	{
		tmp = strsplit(i,"=")
		key <- tmp[[1]][1]
		value <- tmp[[1]][2]
		args[[key]] <- value
	}
	
	fname <- paste(args$token,'_objs.RData.gz',sep='')
	
	save(..., list = list,
			file = fname,
			ascii = FALSE, version = NULL, envir = parent.frame(),
			compress = 'gzip', compression_level = 6,
			eval.promises = TRUE, precheck = TRUE)
}

# Wrapper for load-method, reads object for given run token

objects.get <- function(token){

	# Try locally first
	fname <- paste(token,'_objs.RData.gz',sep='')	
	
	if (file.exists(fname)) {
		load(fname, .GlobalEnv)
	} else {
		# And then via web server
		fname <- paste('http://cl1.opasnet.org/rtools_server/runs/',token,'_objs.RData.gz',sep='')	
		load(url(fname), .GlobalEnv)
	}
}

# New method for storing objects, writes key to the opasnet base as well
# Returns the key

objects.store <- function(..., list = character(), verbose = FALSE){
	
	# Parse arguments
	targs <- strsplit(commandArgs(trailingOnly = TRUE),",")
	args = list()
	
	if (length(targs) == 0) stop('This function can be used within Opasnet only!!!')
	
	for(i in targs[[1]])
	{
		tmp = strsplit(i,"=")
		key <- tmp[[1]][1]
		value <- tmp[[1]][2]
		args[[key]] <- value
	}
	
	now <- Sys.time()
	okey <- gsub("\\.","",as.character(as.numeric(now)))
	okey <- substr(okey,0,12)
	
	if (is.null(args$code_name)) stop('R-code block must have NAME to save objects!')
	
	# Write to base
	data <- matrix(c(args$wiki_page_id, args$code_name, format(now,"%Y-%m-%dT%I:%M:%OS2Z",tz='GMT'), okey), ncol=4, byrow=TRUE)
	colnames(data) <- c("Page ident","Code name","Time","result")
	data <- as.data.frame(data)
	
	index_types <- c("entity","entity","time")
	
	obj_name <- "Saved R objects"
	unit <- "#"
	who <- 'RTools'
	ident <- objects.page_ident(args$wiki_page_id)
	
	if (verbose) paste('Objects page ident:',print(ident),sep=' ')
	if (verbose) paste('Data to insert:',print(data),sep=' ')
	
	if (opbase.obj.exists(ident)){
		opbase.upload(input = data, ident = ident, name = obj_name, act_type = 'append', unit = unit, who = who, verbose = verbose)
	} else {
		opbase.upload(input = data, ident = ident, name = obj_name, act_type = 'replace', unit = unit, who = who, index_types = index_types, verbose = verbose)
	}
	
	# Now finally write objects
	fname <- paste(okey,'_objs.RData.gz',sep='')
	
	save(..., list = list,
			file = fname,
			ascii = FALSE, version = NULL, envir = parent.frame(),
			compress = 'gzip', compression_level = 6,
			eval.promises = TRUE, precheck = TRUE)	
	
	return(okey)
}


objects.latest <- function(page_ident, code_name, verbose = FALSE){
	
	ident <- objects.page_ident(page_ident)
	
	if (verbose) print(paste('Saved R objects page ident is ', ident, sep=''))
	
	series <- opbase.series(ident)
	
	if (verbose) print(paste('Series ids: ',paste(series, collapse=','), sep=''))
	
	res <- NULL
	
	for (s in series)
	{
		tmp <- tryCatch(opbase.data(ident, series_id = s, include = list('Page ident' = page_ident, 'Code name' = code_name), verbose = verbose), error = function(e) return(NULL))

		if (verbose) print(tmp)
		
		if (! is.null(tmp))
		{				
			if (is.null(res))
			{
				res <- tmp
			} else {
				res <- rbind(res, tmp)
			}
		}
	}
	
	if (is.null(res)) stop(paste("No stored objects found! Run initiation code first? Page ident: ",page_ident, " Code name: ", code_name, sep=''))
	
	k <- max(res$Result)
	
	if (verbose) print(paste('Object key is ', k, sep=''))
	
	objects.get(k)
}

# Private function for getting the ident for page holding the key data
objects.page_ident <- function(ident){
	
	ident <- tolower(ident)
	
	# Wiki id
	if (substr(ident, 1,5)=="op_en") return('Op_en5897')
	if (substr(ident, 1,5)=="op_fi") return('Op_fi3382')
	if (substr(ident, 1,6)=="heande")  return('Heande3827')
	if (substr(ident, 1,4)=="test")  return('test4228')
	
	stop(paste("Wiki for ident not determined: ",ident,sep=''))
}

