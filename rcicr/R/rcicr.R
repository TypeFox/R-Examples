# Script by Ron Dotsch, based on Matlab code by Oliver Langner and Python code by Ron Dotsch
# r.dotsch@psych.ru.nl

#' Generate single sinusoid patch
#'
#' @export
#' @import aspace
#' @import matlab
#' @param img_size Integer specifying size of sinusoid patch in number of pixels
#' @param cycles Integer specifying number of cycles sinusoid should span
#' @param angle Value specifying the angle (rotation) of the sinusoid
#' @param phase Value specifying phase of sinusoid
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid
#' @return The sinusoid image with size \code{img_size}.
#' @examples
#' generateSinusoid(512, 2, 90, pi/2, 1.0)
generateSinusoid <- function(img_size, cycles, angle, phase, contrast) {
  
  # Generates an image matrix containing a sinusoid, angle (in degrees) of 0 will give vertical, 90 horizontally oriented sinusoid  
  angle <- aspace::as_radians(angle)
  sinepatch = matlab::repmat(matlab::linspace(0, cycles, img_size), img_size, 1)
  sinusoid <- (sinepatch * cos(angle) + t(sinepatch) * sin(angle)) * 2 * pi
  sinusoid <- contrast * sin(sinusoid + phase)
  return(sinusoid)
}


#' Generate sinusoid noise pattern
#' 
#' @export
#' @import matlab
#' @param img_size Integer specifying size of the noise pattern in number of pixels
#' @param pre_0.3.0 Boolean specifying whether the noise pattern should be created in a way compatible with older versions of rcicr (< 0.3.0). If you are starting a new project, you should keep this at the default setting (FALSE). There is no reason to set this to TRUE, with the sole exception to recreate behavior of rcicr prior to version 0.3.0.  
#' @return List with two elements: the 3D sinusoid matrix with size \code{img_size}, and and indexing
#' matrix with the same size to easily change contrasts.
#' @examples
#' generateNoisePattern(256)
generateNoisePattern <- function(img_size=512, pre_0.3.0=FALSE) {
  # Settings of sinusoids
  scales <- c(1, 2, 4, 8, 16)
  orientations <- c(0, 30, 60, 90, 120, 150)
  phases <- c(0, pi/2)
  
  # Size of sinusoids per scale
  mg <- matlab::meshgrid(1:img_size, 1:img_size,1:length(scales))
  x <- mg$x
  y <- mg$y
  rm(mg)
  sinSize = x / y
  
  # Number of sinsoids needed
  nrSin = length(scales) * length(orientations) * length(phases)
  
  # Pre allocate memory
  sinusoids = matlab::zeros(c(img_size, img_size, nrSin))
  sinIdx = matlab::zeros(c(img_size, img_size, nrSin))
  
  # counters
  
  if (pre_0.3.0) {
    co = 0 # sinusoid layer counter
    idx = 0 # contrast index counter
  } else {
    co = 1 # sinusoid layer counter
    idx = 1 # contrast index counter
  }
  
  for (scale in scales) {
    for (orientation in orientations) {
      for (phase in phases) {
        # Generate single sinusoid
        size <- sinSize[scale, img_size]
        s <- generateSinusoid(size, 2, orientation, phase, 1)
        
        # Repeat to fill scale
        sinusoids[,,co] <- matlab::repmat(s, scale, scale)
        
        # Create index matrix
        for (col in 1:scale) {
          for (row in 1:scale) {
            
            # insert absolute index for later contrast weighting
            sinIdx[(size * (row-1) + 1) : (size * row), (size * (col-1) + 1) : (size * col), co] = idx
            
            # Update contrast counter
            idx = idx + 1
            
          }
        } 
        
        # Update layer counter
        co = co + 1
        
      }
    }
  }
  return(list(sinusoids=sinusoids, sinIdx=sinIdx))
}


#' Generate single noise image based on parameter vector
#' 
#' @export
#' @param params Vector with 4092 values specifying the contrast of each sinusoid in noise
#' @param s 3D sinusoid matrix (generated using \code{generateNoisePattern()})
#' @return The noise pattern as pixel matrix
#' @examples
#' #params <- rnorm(4092) # generates 4092 normally distributed random values
#' #s <- generateNoisePattern(img_size=256)
#' #noise <- generateNoiseImage(params, s)
generateNoiseImage <- function(params, s) {
  noise <- apply(s$sinusoids * array(params[s$sinIdx], dim(s$sinusoids)), 1:2, mean)
  return(noise)
}

#' Generate classification image noise pattern based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid
#' 
#' @export
#' @param stimuli Matrix with one row per trial, each row containing the 4092 parameters for the original stimulus
#' @param responses Vector containing the response to each trial (1 if participant selected original , -1 if participant selected inverted;
#' this can be changed into a scale)
#' @param s 3D sinusoid matrix (generated using \code{generateNoisePattern()})
#' @return The classification image as pixel matrix
generateCINoise <- function(stimuli, responses, s) {
  
  weighted <- responses * stimuli
  
  # Only aggregate if more than one stimulus/response row
  if(is.null(dim(weighted))) {
    params <- weighted
  } else{
    params <- colMeans(weighted)
  }
  
  return(generateNoiseImage(params, s))
}

#' Determines optimal scaling constant for a list of ci's
#' 
#' @export
#' @import matlab
#' @import jpeg
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to jpegs, als the base image (\code{$base})
#' @param saveasjpegs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as jpegs (using the key of the list as name)
#' @param targetpath Optional string specifying path to save jpegs to (default: ./cis)
#' @return List of scaled noise patterns and determind scaling factor
autoscale <- function(cis, saveasjpegs=TRUE, targetpath='./cis') {
  
  # Get range of each ci
  ranges <- matlab::zeros(length(names(cis)), 2)
  for (ciname in names(cis)) {
    ranges[which(ciname==names(cis)), ] <- range(cis[[ciname]]$ci)
  }  

  # Determine the lowest possible scaling factor constant
  if (abs(min(ranges[,1])) > max(ranges[,2])) {
    constant <- abs(min(ranges[,1]))
  }  else {
    constant <- max(ranges[,2])
  }

  print(paste0("Using scaling factor constant:", constant))
  
  # Scale all noise patterns
  for (ciname in names(cis)) {
    cis[[ciname]]$scaled <-  (cis[[ciname]]$ci + constant) / (2*constant)
    
    # Combine and save to jpeg if necessary
    if (saveasjpegs) {
      ci <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2
      
      dir.create(targetpath, recursive=T, showWarnings = F)
      
      jpeg::writeJPEG(ci, paste0(targetpath, '/', ciname, '_autoscaled.jpg'), quality=1.0)
    }
  
  }

  cis[['autoscaling.constant']] <- constant
  return(cis)
}

#' Generates classification image 
#' 
#' Generate classification image for for any reverse correlation task. 
#' 
#' This funcions saves the classification image as jpeg to a folder and returns the CI. Your choice of scaling
#' matters. The default is \code{'matched'}, and will match the range of the intensity of the pixels to
#' the range of the base image pixels. This scaling is non linear and depends on the range of both base image
#' and noise pattern. It is truly suboptimal, because it shifts the 0 point of the noise (that is, pixels that would
#' have not changed base image at all before scaling may change the base image after scaling and vice versa). It is
#' however the quick and dirty way to see how the CI noise affects the base image.
#' 
#' For more control, use \code{'constant'} scaling, where the scaling is independent of 
#' the base image and noise range, but where the choice of constant is arbitrary (provided by the user with t
#' the \code{constant} parameter). The noise is then scale as follows: \code{scaled <- (ci + constant) / (2*constant)}.
#' Note that pixels can take intensity values between 0 and 1 If your scaled noise exceeds those values,
#' a warning will be given. You should pick a higher constant (but do so consistently for different classification images
#' that you want to compare). The higher the constant, the less visible the noise will be in the resulting image.
#' 
#' When creating multiple classification images a good strategy is to find the lowest constant that works for all 
#' classification images. This can be automatized using the \code{autoscale} function.
#' 
#' @export
#' @import jpeg
#' @importFrom stats aggregate
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveasjpeg Boolean stating whether to additionally save the CI as jpeg image
#' @param targetpath Optional string specifying path to save jpegs to (default: ./cis)
#' @param filename Optional string to specify a file name for the jpeg image
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, or \code{matched} (default)
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'})
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined 
generateCI <- function(stimuli, responses, baseimage, rdata, saveasjpeg=TRUE, filename='', targetpath='./cis', antiCI=FALSE, scaling='constant', constant=0.1) {
  
  # Load parameter file (created when generating stimuli)
  load(rdata)
  
  # Check whether critical variables have been loaded
  if (!exists('s', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain s variable.')
  }
  
  if (!exists('base_faces', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain base_faces variable.')
  }
  
  if (!exists('stimuli_params', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain stimuli_params variable.')
  }
  
  # Get base image
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    stop(paste0('File specified in rdata argument did not contain any reference to base image label: ', baseimage, ' (NOTE: file contains references to the following base image label(s): ', paste(names(base_faces), collapse=', '), ')'))
  }
  
  
  # Average responses for each presented stimulus (in case stimuli have been presented multiple times,
  # or group-wise classification images are being calculated, in order to reduce memory and processing
  # load)
  aggregated <- aggregate(responses, by=list(stimuli=stimuli), FUN=mean)
  responses <- aggregated$x
  stimuli <- aggregated$stimuli
  
  # Retrieve parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]
  
  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }
  
  # Compute classification image
  if (antiCI) {
    params = -params
  } 
  ci <- generateCINoise(params, responses, s)
  
  # Scale 
  if (scaling == 'none') {
    scaled <- ci
  } else if (scaling == 'constant') {
    scaled <- (ci + constant) / (2*constant)
    if (max(scaled) > 1.0 | min(scaled) < 0) {
      warning('Chosen constant value for constant scaling made noise of classification image exceed possible intensity range of pixels (<0 or >1). Choose a lower value, or clipping will occur.')
    } 
  } else if (scaling == 'matched') {
    scaled <- min(base) + ((max(base) - min(base)) * (ci - min(ci)) / (max(ci) - min(ci)))
  } else {
    warning(paste0('Scaling method \'', scaling, '\' not found. Using none.'))
    scaled <- ci
  }
  
  # Combine with base image
  combined <- (scaled + base) / 2
  
  # Save to file
  if (saveasjpeg) {
    
    if (filename == '') {
      filename <- paste0(baseimage, '.jpg')
    }
    
    if (antiCI) {
      filename <- paste0('antici_', filename)
    } else {
      filename <- paste0('ci_', filename)
    }
    
    dir.create(targetpath, recursive=T, showWarnings = F)
    
    jpeg::writeJPEG(combined, paste0(targetpath, '/', filename))
    
  }
  
  # Return list
  return(list(ci=ci, scaled=scaled, base=base, combined=combined))
}

#' Generates multiple classification images by participant or condition 
#' 
#' Generate classification image for any reverse correlation task that displays independently generated alternatives. 
#' 
#' This funcions saves the classification images by participant or condition as jpeg to a folder and returns the CIs.
#' 
#' @export
#' @import dplyr
#' @param data Data frame 
#' @param by String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for
#' @param stimuli String specifying column name in data frame that contains the stimulus numbers of the presented stimuli
#' @param responses String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveasjpeg Boolean stating whether to additionally save the CI as jpeg image
#' @param saveunscaledjpeg Optional boolean specifying whether unscaled versions of classification images should be saved as jpeg
#' @param targetpath Optional string specifying path to save jpegs to (default: ./cis)
#' @param label Optional string to insert in file names of jepgs to make them easier to identify 
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant},  \code{matched} or \code{autoscale} (default)
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'})
#' @return List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined) 
batchGenerateCI <- function(data, by, stimuli, responses, baseimage, rdata, saveasjpeg=TRUE, saveunscaledjpeg=FALSE, targetpath='./cis', label='', antiCI=FALSE, scaling='autoscale', constant=0.1) {
  
  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }
  
  pb <- dplyr::progress_estimated(length(unique(data[,by])))
  cis <- list()

  for (unit in unique(data[,by])) {
    
    # Update progress bar
    pb$tick()$print()
    
    # Get subset of data 
    unitdata <- data[data[,by] == unit, ]
    
    # Specify filename for CI jpeg
    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[1,by])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[1,by])
    }

    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI(unitdata[,stimuli], unitdata[,responses], baseimage, rdata, saveunscaledjpeg, paste0(filename, '.jpg'), targetpath, antiCI, scaling, constant)
  }
  
  if (doAutoscale) {
    cis <- autoscale(cis, saveasjpegs=saveasjpeg, targetpath=targetpath)
  }
  
  pb$stop()
  return(cis)
  
}




#' Computes cumulative trial CIs correlations with final/target CI
#' 
#' Computes cumulative trial CIs correlations with final/target CI.
#' 
#' Use for instance for plotting curves of trial-final/target CI correlations to estimate how many trials are necessary in your task
#' 
#' @export
#' @import dplyr
#' @importFrom stats cor
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param targetci List Target CI object generated with rcicr functions to correlate cumulative CIs with
#' @param step Step size in sequence of trials to compute correlations with
#' @return Vector containing correlation between cumulative CI and final/target CI 
computeCumulativeCICorrelation <- function(stimuli, responses, baseimage, rdata, targetci=list(), step=1) {
  
  # Load parameter file (created when generating stimuli)
  load(rdata)
  
  # Check whether critical variables have been loaded
  if (!exists('s', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain s variable.')
  }
  
  if (!exists('base_faces', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain base_faces variable.')
  }
  
  if (!exists('stimuli_params', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain stimuli_params variable.')
  }
  
  # Get base image
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    stop(paste0('File specified in rdata argument did not contain any reference to base image label: ', baseimage, ' (NOTE: file contains references to the following base image label(s): ', paste(names(base_faces), collapse=', '), ')'))
  }
  
  
  # Retrieve parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]
  
  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }
  
  # Compute final classification image if necessary
  if (length(targetci) == 0) {
    finalCI <- generateCINoise(params, responses, s)
  } else {
    finalCI <- targetci$ci
  }
  
  # Compute correlations with final CI with cumulative CI
  pb <- dplyr::progress_estimated(length(responses))
  
  correlations <- vector()
  corcounter <- 1
  for (trial in seq(1,length(responses), step)) {
    pb$tick()$print()
    
    cumCI <- generateCINoise(params[1:trial,], responses[1:trial], s)
    correlations[corcounter] <- cor(as.vector(cumCI), as.vector(finalCI))
    corcounter <- corcounter + 1
  }
  pb$stop()
 
  # Return correlations
  return(correlations)
}