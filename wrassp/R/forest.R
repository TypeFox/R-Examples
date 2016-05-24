##' forest function adapted from libassp
##'
##' Formant estimation of the signal(s) in <listOfFiles>.
##' Raw resonance frequency and bandwidth values are
##' obtained by root-solving of the Linear Prediction
##' polynomial from the autocorrelation method and the
##' Split-Levinson-Algorithm (SLA). Resonances are then
##' classified as formants using the so-called Pisarenko
##' frequencies (by-product of the SLA) and a formant
##' frequeny range table derived from the nominal F1
##' frequency. The latter may have to be increased by
##' about 12\% for female voices (see NominalF1 and Gender options).
##' Formant estimates will be written to a file with the
##' base name of the input file and extension '.fms'.
##' Default output is in SSFF binary format (tracks 'fm'
##' and 'bw')
##' @title forest
##' @param listOfFiles vector of file paths to be processed by function
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds (default = 0: begin of data)
##' @param endTime = <time>:  set end of analysis interval to <time> seconds (default = 0: end of data)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param windowSize  = <dur>: set analysis window size to <dur> ms (default: 30.0)
##' @param effectiveLength make window size effective rather than exact
##' @param nominalF1 = <freq>: set nominal F1 frequency to <freq> Hz (default: 500.0 Hz)
##' @param gender = <code>: set gender specific parameters where 
##' <code> = f[emale], m[ale] or u[nknown] (when <code>=f: eff. window length = 12.5 ms nominal F1 = 560.0 Hz)
##' @param estimate insert rough frequency estimates of missing formants (default: frequency set to zero)
##' @param order decrease default order by 2 (one resonance less)
##' @param incrOrder increase default order by 2 (one resonance more)
##' @param numFormants = <num>: set number of formants to <num> (default: 4;  maximum: 8 or half the LP order)
##' @param window = <type>: set analysis window function to <type> (default: BLACKMAN)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (-1 <= val <= 0) 
##' (default: dependent on sample rate and nominal F1)
##' @param toFile write results to file (default extension is .fms)
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
##' # calculate formant values
##' res <- forest(path2wav, toFile=FALSE)
##' 
##' # plot formant values
##' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) + 
##'           attr(res, 'startTime'), 
##'         res$fm, 
##'         type='l', 
##'         xlab='time (s)', 
##'         ylab='Formant frequency (Hz)')
##' 
##' @export
'forest' <- function(listOfFiles = NULL, optLogFilePath = NULL,
                     beginTime = 0.0, endTime = 0.0, 
                     windowShift = 5.0, windowSize = 20.0, 
                     effectiveLength = TRUE, nominalF1 = 500, 
                     gender = 'm', estimate = FALSE, 
                     order = 0, incrOrder = 0, 
                     numFormants = 4, window = 'BLACKMAN', 
                     preemphasis = -0.8, toFile = TRUE, 
                     explicitExt = NULL, outputDirectory = NULL, 
                     forceToLog = useWrasspLogger){
	
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
	  OutputDirectory = normalizePath(path.expand(outputDirectory))
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
	#perform analysis
	
	if(length(listOfFiles)==1){
    pb <- NULL
	}else{
	  if(toFile==FALSE){
	    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
	  }
    cat('\n  INFO: applying forest to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
	}

	externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "forest", beginTime =  beginTime, 
                                    endTime = endTime, windowShift = windowShift, 
                                    windowSize = windowSize, effectiveLength = effectiveLength, 
                                    nominalF1 = nominalF1, gender = gender, 
                                    estimate = estimate, order = as.integer(order), 
                                    incrOrder = as.integer(incrOrder), numFormants = as.integer(numFormants), 
                                    window = window, preemphasis = preemphasis, 
                                    toFile = toFile, explicitExt = explicitExt, 
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
