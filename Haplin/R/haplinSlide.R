haplinSlide <- function(filename, data, pedIndex, markers = "ALL", winlength = 1, strata = NULL, table.output = TRUE, cpus = 1, slaveOutfile = "", printout = FALSE, verbose = FALSE, ...)
{
##
## RUN HAPLIN ON SLIDING WINDOWS
##
## MERK: BRUKER "#i#" TIL AA KOMMENTERE VEKK DET SOM GIR ADVARSEL MED R CMD check, SIDEN DE RETTE PAKKENE IKKE ER INSTALLERT
#
## SELECT PARALLEL MECHANISM (SOME MAY BE DISABLED)
#.para <- "snow"
#.para <- "snowfall"
#.para <- "doSMP"
.para <- "parallel"
#
## RUNS IN PARALLEL ONLY IF cpus IS NUMERIC > 1
if(!missing(cpus)){
	if(!is.numeric(cpus)) stop('The number of cpu-s "cpus" must be numeric!', call. = F)
	.run.para <- (cpus > 1)
} else .run.para <- F
#
## GET HAPLIN DEFAULTS, OVERRIDE DEFAULTS FOR SPECIFIED ARGUMENTS
.info <- f.catch(match.call(), formals())
#
##
if(winlength > 1) cat("\nImportant: Remember that SNPs must be in correct physical ordering\nfor a sliding window approach to be meaningful!\n")
#
## CHECK IF DATA DERIVE FROM AN R OBJECT OR FROM FILE
.missdata <- missing(data)
.misspedIndex <- missing(pedIndex)
#
## IMPORTANT: DO NOT REMOVE MISSING YET, EVEN IF use.missing = F (I.E. REMOVE WINDOW-WISE, NOT LISTWISE)
.use.missing <- .info$model$use.missing ## SAVE OLD VALUE
.info$model$use.missing <- T ## ENFORCE
#
## WITH A gwaa.data OBJECT, DELAY CONVERSION, ELSE READ AND CONVERT
.is.gwaa <- !.missdata && (class(data) == "gwaa.data")
if(.is.gwaa){
	## ONLY REDUCE NUMBER OF MARKERS, BUT DO NOT CONVERT YET
	.data.read <- data[, .info$filespecs$markers]
}else{
	## READ AND CONVERT DATA
	.data.read <- f.get.data(data, pedIndex, .info)
	.info <- attr(.data.read, "info")
}
.marker.names <- f.get.marker.names(.data.read, n.vars = .info$filespecs$n.vars)
if(length(.marker.names) != length(.info$filespecs$markers)) warning("Something's wrong with the marker names...", call. = F)
#
##
.info$model$use.missing <- .use.missing ## REVERT TO ORIGINAL VALUE
#
## FIND MARKERS INCLUDED IN EACH WINDOW. NB: THEY NOW REFER TO MARKERS IN THE _REDUCED_ FILE, NOT THE ORIGINAL ONE
.slides <- f.windows(markers = seq(along = .info$filespecs$markers), winlength = winlength)
## USE ORIGINAL MARKER SPECIFICATION AS NAMES
###.names <- .info$filespecs$markers[.slides]
.names <- .marker.names[.slides]
.names <- matrix(.names, ncol = ncol(.slides))
.names <- f.create.tag(.names, sep = "-")
#
## PREPARE TO RUN ON EACH WINDOW
.nres <- dim(.slides)[1]
#
## REPRODUCE HAPLIN ARGUMENTS FROM .info
.args <- f.args.from.info(.info)
#
##
cat("\n")
#
## FUNCTION TO RUN HAPLIN ON A SINGLE WINDOW
.f.run <- function(i){
	cat("Running Haplin on Window '", .names[i], "' (", i, "/", .nres, ")...:  ", sep = "")
	args_ <- .args
	#
	if(.run.para){
		## HAPLIN MUST BE MADE AVAILABLE IN EACH RUN
		suppressPackageStartupMessages(require(Haplin, quietly = T))
	}
	#
	##
	# if(.info$filespecs$database){
	# 	## HAPLIN DATABASE
	# 	args_$markers <- .info$filespecs$markers[.slides[i,]]
	# }else{
		## FILE, DATA MATRIX, AND GWAA
		args_$filename <- NULL
		args_$markers <- .slides[i,]
		args_$data <- .data.read
		if(.is.gwaa & !.misspedIndex){
			args_$pedIndex <- pedIndex
		}
	# }
		#
		## RUN HAPLIN
		if(is.null(strata)){
			.res <- try(do.call("haplin", args_), silent = T)
			#
			## CHECK IF ERRORS
			if(class(.res) == "try-error") cat("RUN FAILED\n")
			else{
				cat("done\n")
				if(table.output) .res <- haptable(.res)
			}		
		}else{
			.res <- try(do.call("haplinStrat", args_), silent = T)
			## CHECK IF ERRORS
			if(class(.res) == "try-error") cat("RUN FAILED\n")
			else{
				cat("done\n")
				#.res <- posttest(.res)
			}		
		}
		return(.res)
}
#
## DO THE RUN
## EITHER IN SEQUENCE (SINGLE CPU) 
## OR IN PARALLEL (MULTIPLE CPUs)
if(!.run.para){
    if(!missing(slaveOutfile)) sink(file = slaveOutfile)
	.res.list <- lapply(seq(length.out = .nres), .f.run)
    if(!missing(slaveOutfile)) sink()
	names(.res.list) <- .names
	#
	## CHECK FOR HAPLIN FAILURES
	.errs <- sapply(.res.list, class) == "try-error"
}else{
	#
	#
	#	## INITIALIZE CPUS
	if(.para == "snowfall"){
		## SNOWFALL
#i#		sfInit(parallel = T, slaveOutfile = slaveOutfile, cpus = cpus)
#i#		on.exit(sfStop())
	}
	if(.para == "parallel"){
		## parallel
		w <- parallel::makeCluster(spec = cpus, outfile = slaveOutfile)
		#on.exit(stopCluster(w))
		#
		## shut down nodes explicitly
		.plat <- .Platform$OS.type
		if(.plat == "windows") .cmds <- "taskkill /pid "
		if(.plat == "unix") .cmds <- "kill "
		# get process id's
		.pids <- sapply(w, function(x){
			snow::sendCall(x, eval, list(quote(Sys.getpid())))
			snow::recvResult(x)
		})
		.cmds <- paste(.cmds, .pids, sep = "")
		on.exit({
			if(!missing(slaveOutfile)) sink(file = slaveOutfile, append = T)
				sapply(.cmds, system)
				for(.w in w){
					try(snow::postNode(.w, "DONE"), silent = T)
					try(snow::closeNode(.w), silent = T)
				}
			if(!missing(slaveOutfile)) sink()
		})
	}
	if(.para == "doSMP"){
		## doSMP BACKEND TO foreach
#i#		w <- startWorkers(workerCount = cpus)
#i#		on.exit(stopWorkers(w))
#i#		registerDoSMP(w)
	}
	if(.para == "snow"){
		## SNOW
		# w <- makeCluster(cpus, type = type)
#i#		w <- makeCluster(cpus)
#i#		on.exit(stopCluster(w))
	}
	#
	## RUN BY SPLITTING ON CPUS
	if(!missing(slaveOutfile)) sink(file = slaveOutfile, append = T)
		cat("\n--- Running haplinSlide using ", cpus, " cpu's ---\n\n", sep = "")
	if(!missing(slaveOutfile)) sink()
#i#	if(.para == "snowfall")	.res.list <- sfLapply(seq(length.out = .nres), .f.run)
	if(.para == "parallel") .res.list <- parLapply(w, seq(length.out = .nres), .f.run)
#i#	if(.para == "snow") .res.list <- parLapply(w, seq(length.out = .nres), .f.run)
#i#	if(.para == "doSMP") .res.list <- foreach(i = seq(length.out = .nres)) %dopar% .f.run(i)
	names(.res.list) <- .names
	#
	## CHECK FOR HAPLIN FAILURES
	.ftest <- function(x){
		return(identical(class(x), "try-error"))
	}
#i#	if(.para == "snowfall") .errs <- sfSapply(.res.list, .ftest)
	if(.para == "parallel") .errs <- unlist(parLapply(w, .res.list, .ftest))
#i#	if(.para == "snow") .errs <- parLapply(w, .res.list, class) == "try-error"
#i#	if(.para == "doSMP") .errs <- foreach(i = seq(along = .res.list)) %dopar% class(.res.list[[i]]) == "try-error"
}
#
## GO THROUGH POSSIBLE ERRORS
if(any(.errs)){
	## COLLECT ERROR MESSAGES, CHANGE OUTPUT TO NA WITH ERR. MESS. AS ATTRIBUTES
	.mess <- .res.list[.errs]
	.err.res <- lapply(.mess, function(x) {
		.tmpres <- NA
		attributes(x) <- NULL # REMOVE "try-error"
		attr(.tmpres, "error.message") <- x
		return(.tmpres)
	})
	.res.list[.errs] <- .err.res
}
if(!missing(slaveOutfile)) sink(file = slaveOutfile, append = T)
	cat("\n --- haplinSlide has completed ---\n")
if(!missing(slaveOutfile)) sink()
#
## SET CLASS
###if(!table.output) 
class(.res.list) <- "haplinSlide"
#
## RETURN RESULT LIST
return(.res.list)
}
