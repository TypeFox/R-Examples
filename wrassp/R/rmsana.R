##' rmsana function adapted from libassp
##'
##' Analysis of short-term Root Mean Square amplitude of
##' the signal in <listOfFiles>. Per default, the RMS values are
##' expressed in decibel (dB) so that they correspond to
##' the short-term power of the signal.
##' Analysis results will be written to a file with the
##' base name of the input file and extension '.rms'.
##' Default output is in SSFF binary format (track 'rms').
##' @title rmsana
##' @param listOfFiles vector of file paths to be processed by function
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>:  set begin of analysis interval to <time> seconds (default = 0: begin of file)
##' @param centerTime = <time>: set single-frame analysis with the analysis window centred at <time> seconds; 
##' overrules beginTime, endTime and windowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds (default: end of file)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param windowSize = <dur>: set analysis window size to <dur> ms; overrules effectiveLength option
##' @param effectiveLength make window size effective rather than exact
##' @param linear calculate linear RMS values (default: values in dB)
##' @param window = <type>: set analysis window function to <type> (default: HAMMING)
##' @param toFile write results to file (default extension is .rms)
##' @param explicitExt set if you wish to overwride the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param forceToLog is set by the global package variable useWrasspLogger. This is set
##' to FALSE by default and should be set to TRUE is logging is desired.
##' @return nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @useDynLib wrassp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"), 
##'                        pattern = glob2rx("*.wav"), 
##'                        full.names = TRUE)[1]
##' 
##' # calculate rms values
##' res <- rmsana(path2wav, toFile=FALSE)
##' 
##' # plot rms values
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res$rms, 
##'      type='l', 
##'      xlab='time (s)', 
##'      ylab='RMS energy (dB)')
##' 
##' @export
'rmsana' <- function(listOfFiles = NULL, optLogFilePath = NULL,
                     beginTime = 0.0, centerTime = FALSE, 
                     endTime = 0.0, windowShift = 5.0, 
                     windowSize = 20.0, effectiveLength = TRUE, 
                     linear = FALSE, window = 'HAMMING', 
                     toFile = TRUE, explicitExt = NULL,
                     outputDirectory = NULL, forceToLog = useWrasspLogger){


	###########################
	# a few parameter checks and expand paths
	
	if (is.null(listOfFiles)) {
		stop(paste("listOfFiles is NULL! It has to be a string or vector of file",
		           "paths (min length = 1) pointing to valid file(s) to perform",
		           "the given analysis function."))
	}

	if (is.null(optLogFilePath) && forceToLog){
	  stop("optLogFilePath is NULL! -> not logging!")
	}else{
	  if(forceToLog){
	    optLogFilePath = path.expand(optLogFilePath)  
	  }
	}
	
	if(!isAsspWindowType(window)){
		stop("WindowFunction of type '", window,"' is not supported!")
	}

	if (!is.null(outputDirectory)) {
	  outputDirectory = normalizePath(path.expand(outputDirectory))
	  finfo  <- file.info(outputDirectory)
	  if (is.na(finfo$isdir))
	    if (!dir.create(outputDirectory, recursive=TRUE))
	      stop('Unable to create output directory.')
	  else if (!finfo$isdir)
	    stop(paste(outputDirectory, 'exists but is not a directory.'))
	}
  
    ###########################
    # Pre-process file list
    listOfFiles <- prepareFiles(listOfFiles)
  
	###########################
	# perform analysis

	if(length(listOfFiles)==1){
          pb <- NULL
	}else{
	  if(toFile==FALSE){
	    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
	  }
    cat('\n  INFO: applying rmsana to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
	}	
	
	externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "rmsana", beginTime = beginTime, 
                                    centerTime = centerTime, endTime = endTime, 
                                    windowShift = windowShift, windowSize = windowSize, 
                                    effectiveLength = effectiveLength, linear = linear, 
                                    window = window, toFile = toFile, 
                                    explicitExt = explicitExt, 
                                    progressBar = pb, outputDirectory = outputDirectory,
                                    PACKAGE = "wrassp"))
	
  ############################
  # write options to options log file  
	if (forceToLog){
	  optionsGivenAsArgs = as.list(match.call(expand.dots = TRUE))
	  wrassp.logger(optionsGivenAsArgs[[1]], optionsGivenAsArgs[-1],
	                optLogFilePath, listOfFiles)
        
	}        
  #############################
  # return dataObj if length only one file
        
	if(!(length(listOfFiles)==1)){
    close(pb)
  }else{
    return(externalRes)
  }      
}

