enve.cliopts <- function(
      ### Generates nicely formatted command-line interfaces for
      ### functions (_closures_ only).
      fx,
      ### Function for which the interface should be generated.
      rd_file,
      ### (Optional) .Rd file with the standard documentation of the function.
      positional_arguments,
      ### (Optional) Number of _positional_ arguments passed to parse_args
      ### (package:optparse).
      usage,
      ### (Optional) Usage passed to OptionParser (package:optparse).
      mandatory=c(),
      ### Mandatory arguments.
      vectorize=c(),
      ### Arguments of the function to vectorize (comma-delimited). If numeric,
      ### use also `number`.
      ignore=c(),
      ### Arguments of the function to ignore.
      number=c(),
      ### Force these arguments as numerics. Useful for numeric
      ### vectors (see `vectorize`) or arguments with no defaults.
      defaults=list(),
      ### Defaults to use instead of the ones provided by the formals.
      o_desc=list(),
      ### Descriptions of the options. Help from `rd` is ignored for arguments
      ### present in this list.
      p_desc=""
      ### Description of the function. Help from `rd` is ignored for the
      ### function description unless this value is an empty string.
      ){
   
   #= Load stuff
   if(!suppressPackageStartupMessages(
      requireNamespace("optparse", quietly=TRUE)))
	 stop("Package 'optparse' is required.")
   requireNamespace("tools", quietly=TRUE)
   if(missing(positional_arguments)) positional_arguments <- FALSE
   if(missing(usage)) usage <- "usage: %prog [options]"
   
   #= Get help (if any)
   if(!missing(rd_file)){
      rd <- tools::parse_Rd(rd_file)
      for(i in 1:length(rd)){
         tag <- attr(rd[[i]],'Rd_tag')
	 if(tag=="\\description" && p_desc==""){
	    p_desc <- paste("\n\t",as.character(rd[[i]]),sep='')
	 }else if(tag=="\\arguments"){
	    for(j in 1:length(rd[[i]])){
	       if(length(rd[[i]][[j]])==2){
		  name <- as.character(rd[[i]][[j]][[1]])
		  if(length(o_desc[[name]])==1) next
		  desc <- as.character(rd[[i]][[j]][[2]])
		  o_desc[[name]] <- paste(gsub("\n","\n\t\t",desc), collapse='')
	       }
	    }
	 }
      }
   }

   #= Set options
   o_i <- 0
   opts <- list()
   f <- formals(fx)
   if(length(defaults)>0){
      for(i in 1:length(defaults)) f[[names(defaults)[i]]] <- defaults[[i]]
   }
   for(i in names(f)){
      if(i=="..." || i %in% ignore) next
      o_i <- o_i + 1
      flag <- gsub("\\.","-",i)

      optopt <- list(help="")
      if(length(o_desc[[i]])==1) optopt$help <- o_desc[[i]]
      if(!is.null(f[[i]]) && !suppressWarnings(is.na(f[[i]])) && is.logical(f[[i]])){
	 optopt$opt_str <- paste(ifelse(f[[i]], "--no-", "--"), flag, sep='')
	 optopt$action  <- ifelse(f[[i]], "store_false", "store_true")
      }else{
	 optopt$opt_str <- paste("--", flag, sep='')
	 optopt$action  <- "store"
	 optopt$help <- paste(optopt$help, "\n\t\t[",
	    ifelse(i %in% mandatory, "** MANDATORY", "default %default"),
	    ifelse(i %in% vectorize, ", separate values by commas", ""),
	    "].", sep="")
      }
      if(!is.name(f[[i]])){
	 optopt$default <- f[[i]]
	 optopt$metavar <- class(f[[i]])
      }
      if(i %in% number) optopt$metavar <- "NUMERIC"
      optopt$dest <- i
      
      opts[[o_i]] <- do.call(optparse::make_option, optopt)
   }
   opt <- optparse::parse_args(
      optparse::OptionParser(option_list=opts, description=p_desc, usage=usage),
      positional_arguments=positional_arguments)

   #= Post-hoc checks
   if(length(opt[['options']])==0) opt <- list(options=opt, args=c())
   for(i in mandatory){
      if(length(opt$options[[i]])==0) stop('Missing mandatory argument: ',i)
   }
   for(i in vectorize){
      if(length(opt$options[[i]])==1)
	 opt$options[[i]] <- strsplit(opt$options[[i]],",")[[1]]
   }
   for(i in number){
      if(length(opt$options[[i]])>0)
	 opt$options[[i]] <- as.numeric(opt$options[[i]])
   }
   opt$options$help <- NULL
   
   return(opt)
   ### Returns a `list` with keys: `options`, a named list with the values for
   ### the function's arguments; and `args`, a vector with zero or more strings
   ### containing the positional arguments.
}

