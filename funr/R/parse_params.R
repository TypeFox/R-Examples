

parse_params <- function(func, paramPairs, verbose = FALSE){
	#func <- as.character(paramPairs[1])
	#if(length(func) == 0) return(help())

	args <- formals(func)
	paramPairs <- paramPairs[grep("=", paramPairs)] ## get those with =

	if(verbose){message("args:");print(args)}

	#print("paramPairs:");print(paramPairs)
	#args_supplied = sapply(strsplit(paramPairs, "="), "[[", 1)


	if(verbose)
		message("\nget_params: we have ",
						length(paramPairs), " parameters\n",
						paste(paramPairs, collapse = "\n"))

	for(param in paramPairs){
		if(verbose)
			message("\nstarting process with: ", param)

		if(!grepl("=", param))
			stop("seems you are missing a = sign between variable and value ?")

		## split using = signs
		splt <- unlist(strsplit(param, "="));
		nm = splt[1]
		value = splt[2]

		## check if this is stdin
		if(value %in% c("/dev/stdin", "-"))
			suppressWarnings({value = readLines(file("stdin"))})

		if(verbose)
			print(value)

		value <- strsplit(value,",")[[1]] #handling those with , in value.. for multiple R values

		## --- if function supports ... need to pass ALL arguments
		if(sum( names(args) %in% "...") & !nm %in% names(args) ){
			## -- remove the dots argument
			if(verbose)
				message("Adding ", nm, ":", value)
			l = list(nm = value);names(l) = nm
			args <- c(args, l)
		}

		# 		if(verbose)
		# 			message("processing param: ", nm, " value ", args[[nm]])

		if(class(args[[nm]]) == "numeric" ){
			args[[nm]] = as.numeric(value)
		}else if(class(args[[nm]]) %in% c("logical") | (value[1] %in% c("TRUE", "FALSE") & length(value) == 1)){
			args[[nm]] = as.logical(value)
		}else if(class(args[[nm]]) %in% c("character", "name" )){
			args[[nm]] = as.character(value)
		}else if(class(args[[nm]]) %in% c("list")){
			args[[nm]] = as.list(value)
		}else if(class(args[[nm]]) %in% c("call")){ ## example call to getOption
			args[[nm]] = as.character(value)
		}
		if(verbose)
			message("processed param: ", nm, " value ", args[[nm]])
	}

	## remove dots
	dots = which( names(args) == "..." )
	if(length(dots) > 0 )
		args = args[-dots]

	## check values if NULL, remove.
	rm = which(sapply(args, is.name))
	if(length(rm) > 0)
		args = args[-rm]

	##print(do.call(rbind, as.list(args)))
	if(verbose) print(args)
	return(as.list(args))

}

