##' read.AsspDataObj creates an object of class dobj from a signal or parameter 
##' file readable by the ASSP Library (WAVE, SSFF, AU, ...)
##'
##' @title read.AsspDataObj from a signal/parameter file
##' @param fname filename of the signal or parameter file
##' @param begin begin time (default is in ms) of segment to retrieve
##' @param end end time (default is in ms) of segment to retrieve
##' @param samples (BOOL) if set to false ms values of begin/end are sample numbers
##' @return list object containing file data
##' @author Lasse Bombien
##' @aliases getAsspDataObj
##' @useDynLib wrassp
##' @export
'read.AsspDataObj' <- 'getAsspDataObj' <- function(fname, begin=0, end=0, samples=FALSE) {
  fname <- path.expand(fname)
  .External("getDObj2", fname, begin=begin, end=end, samples=samples, PACKAGE="wrassp")
}

##' Prints an overview of ASSP Data Objects
##'
##' @title print a summary of an AsspDataObj
##' @param x an object of class AsspDataObj
##'
##' @param ... other arguments that might be passed on to other functions 
##' @author Lasse Bombien
##' @method print AsspDataObj
##' @seealso \code{\link{read.AsspDataObj}}
##' @useDynLib wrassp
##' @aliases summary.AsspDataObj
##' @export
"print.AsspDataObj" <- summary.AsspDataObj<- function(x, ...)
{
    temp <- attr(x, "filePath")
    if (is.null(temp)) {
        cat("In-memory Assp Data Object\n")
    }
    else {
        cat(paste("Assp Data Object of file ", temp, ".\n", sep=""))
    }
    cat(sprintf("Format: %s (%s)\n", AsspFileFormat(x), AsspDataFormat(x)))
    cat(paste(as.integer(numRecs.AsspDataObj(x)),
              "records at", attr(x, 'sampleRate'), "Hz\n"))
    cat(sprintf("Duration: %f s\n", dur.AsspDataObj(x)))
    cat(paste("Number of tracks:", length(names(x)), "\n"))
    for (track in names(x)) {
        cat('\t', track)
        cat(paste(" (", ncol(x[[track]]), " fields)\n", sep=''))
    }
    genVars <- attr(x, 'genericVars')
    if (!is.null(genVars)) {
        cat("\nGeneric variables:\n")
        for (var in names(genVars)) {
            cat(sprintf("  %s:", var))
            if (genVars[[var]]$Type %in% c("CHAR", "BYTE")) {
                cat(sprintf("\t%s\n", genVars[[var]]$Value))
            } else {
                cat(sprintf("\t%f\n", genVars[[var]]$Value))
            }
            cat(sprintf("    (%s)\n", genVars[[var]]$Type))
        }
    }
}


##' Writes an object of class AsspDataObj to a file given the meta information
##' contained in the object.
##'
##' @title write.AsspDataObj to file
##' @param dobj an object of class AsspDataObj
##' @param file file name as a character string, defaults to the
##' \code{filePath} attribute of the AsspDataObj
##' @return NULL
##' @author Lasse Bombien
##' @useDynLib wrassp
##' @export
"write.AsspDataObj" <- function (dobj, file=attr(dobj, 'filePath'))
  {
    if (is.null(file))
      stop('File path not set internally. Please specify!')
    file <- path.expand(file)
    .Call("writeDObj", dobj, file, PACKAGE="wrassp")
  }

##' Checks whether x is a valid AsspDataObj
##'
##' @title Checks whether x is a valid AsspDataObj
##' @param x an object of class AsspDataObj
##' @param ... optional other arguments passed to further functions
##' @return TRUE or FALSE
##' @author Lasse Bombien
##' @useDynLib wrassp
##' @export
is.AsspDataObj <- function (x, ...)
  {
    if (class (x) != "AsspDataObj")
      return (FALSE)
    return (TRUE)
  }


##' Remove a track from an
##' AsspDataObj object
##'
##' @title Remove track from an AsspDataObj
##' @param dobj An object of class AsspDataObj
##' @param trackname the name of a track in this object
##' @return The object without the track named trackname
##' @author Lasse Bombien
##' @useDynLib wrassp
##' @export
delTrack <- function (dobj, trackname)
  {
    if (!is.AsspDataObj (dobj))
      stop ('First argument must be a AsspDataObj.')

    w <- which (names (dobj) == trackname)
    if (length (w) != 1)
      stop ('Invalid trackname')

    ## remove track
    dobj[[trackname]] <- NULL
    ## remove
    attr(dobj, 'trackFormats') <- attr(dobj, 'trackFormats')[-w]
    
    return (dobj)
  }

##' Add a track to an AsspDataObj
##'
##' The specified data object is extended by a new track named \code{trackname}.
##' If there already is a track with the same name and \code{deleteExisiting}
##' is \code{FALSE} the function does nothing but returns with an error. If
##' \code{deleteExisting} is \code{TRUE} the existing track will be removed
##' (see \code{\link{delTrack}}.
##' \code{data} to be added is a numeric matrix (or will be coerced to one).
##' It must have
##' the same number of rows as the tracks that already exist in the object
##' (if any). TODO add \code{format} information.
##' @title Add a track to an AsspDataObj
##' @param dobj The data object to which the data is to be added
##' @param trackname The name of the new track
##' @param data a matrix with values
##' @param format format for binary writing to file (defaults to 'INT16') 
##' @param deleteExisting Delete existing track with the same (default: FALSE)
##' @return the object including the new track
##' @author Lasse Bombien
##' @seealso \code{\link{delTrack}}
##' @useDynLib wrassp
##' @export
addTrack <- function (dobj, trackname, data, format = 'INT16',
                      deleteExisting=FALSE) {
  if (!is.AsspDataObj(dobj))
    stop('dobj must be an AsspDataObj.')
  
  if (!is.numeric(data))
    stop('data must be a numeric matrix')
  
  if (!is.character(trackname) | length(trackname) != 1)
    stop('trackname must be an atomic string.')
  
  data <- as.matrix(data)
  
  tracks <- names(dobj)
  w <- tracks  == trackname
  if (any(w) & !deleteExisting)
    stop(paste('Track', trackname,
                'exists and will not be deleted',
                '("deleteExisting" argument)'))
  if (length(tracks) == 1 & any(w)) {
      ## this is fine: the only track will be replaced
  } else if (length(tracks) > 0) {
    if (nrow(data) != nrow(dobj[[1]]))
      stop(paste("number of rows in data must match number of rows in",
                  "existing tracks."))
  }

  dobj[[trackname]] <- data
  if (any(w))
    attr(dobj, 'trackFormats')[w] <- format
  else
    append(attr(dobj, 'trackFormats'), format)

  return(dobj)
}

##' List the tracks of an AsspDataObj
##'
##' AsspDataObj contain tracks (at least one). This function lists the names
##' of these tracks. This function is equvalent to calling \code{names(x)}.
##' @title tracks.AsspDataObj
##' @param x an object of class AsspDataObj
##' @return a character vector containing the names of the tracks
##' @author Lasse Bombien
##' @export
##' @useDynLib wrassp
tracks.AsspDataObj <- function(x) {
  names(x)
}

##' Function to get or set the file format of an AsspDataObj.
##' 
##' \code{libassp} handles a number of file formats common in speech research. 
##' This function enables the user to determine the file format of an object 
##' read from file and to set it for subsequent writing. This allows for file 
##' format conversion to some degree. Note, that many conversions are not 
##' reasonable/possible: conversions are therefore discouraged unless the user 
##' knows what they are doing. Format specifiers can be found in
##' \code{\link{AsspFileFormats}} and exist in two forms: a code name and a
##' code number. Both are suitable for setting the format.
##' @title Get and set AsspFileFormat
##' @param x an object of class AsspDataObj
##' @return for \code{AsspFileFormat} the code name of the object's 
##'   currently set file format
##' @author Lasse Bombien
##' @seealso \code{\link{AsspFileFormats}}, \code{\link{AsspDataFormat}}
##' @examples
##' \dontrun{
##' obj  <- read.AsspDataObj('/path/to/file.wav')
##' AsspFileFormat(obj)
##' AsspFileFormat(obj) <- 'SSFF' ## or
##' AsspFileFormat(obj) <- 20
##' }
##' @useDynLib wrassp
##' @export
AsspFileFormat <- function(x) {
  ## file format is in the first element (of two) in the fileInfo attribute
  xx <- x
  if (!is.AsspDataObj(xx))
    stop('Argument must be an object of class AsspDataObj')
  curFormat <- attr(xx, 'fileInfo')[1]
  ind <- match(curFormat, AsspFileFormats)
  if (is.na(ind))
    stop('Invalid file format. This AsspDataObj has been messed with!')
  return(names(AsspFileFormats)[ind])
}

##' @rdname AsspFileFormat
##' @param value an integer or a string indicating the new file format
##' @usage AsspFileFormat(x)  <- value
##' @return for \code{AsspFileFormat<-}, the updated object
##' @export
"AsspFileFormat<-" <- function(x, value) {
  value <- value[1]
  if (!is.AsspDataObj(x))
    stop('Argument must be an object of class AsspDataObj')
  fi  <- attr(x, 'fileInfo')
  if (is.numeric(value)) {
    ind <- match(value, AsspFileFormats)
  } else if (is.character(value)) {
    ind <- match(value, names(AsspFileFormats))
  } else {
    stop ('format must be an integer or a string.')
  }
  if (is.na(ind))
    stop('format does not specify a valid file format.')
  fi[1]  <- AsspFileFormats[ind]
  attr(x, 'fileInfo')  <- as.integer(fi)
  x
}

##' Function to get or set the data format of an AsspDataObj.
##'
##' \code{libassp} can store data in binary and ASCII format. 
##' This function enables the user to determine the data format of an object 
##' read from file and to set it for subsequent writing.
##' Valid values are 
##' \code{'ascii'} (or \code{1}) for ASCII format or \code{'binary'} (or \code{2}) for binary IO.
##' Use is discouraged unless the user knows what they are doing.
##' @title Get/set data format of an AsspDataObj
##' @param x an object of class AsspDataObj
##' @return a string representing the current data format
##' @useDynLib wrassp
##' @seealso \code{\link{AsspFileFormat}}
##' @export
##' @author Lasse Bombien
AsspDataFormat <- function(x) {
  f <- attr(x, 'fileInfo')[2]
  if (f==1) 
    return('ascii')
  else if (f==2)
    return('binary')
  else
    stop('Invalid data format. This AsspDataObj has been messed with!')
}

##' @rdname AsspDataFormat
##' @param value an integer or a string indicating the new data format
##' @usage AsspDataFormat(x)  <- value
##' @return for \code{AsspDataFormat<-}, the updated object
##' @export
##' 
"AsspDataFormat<-" <- function(x, value) {
  value <- value[1]
  fi <- attr(x, 'fileInfo')
  if (is.numeric(value)) {
    if (value %in% c(1,2))
      fi[2] <- value
    else
      stop('Invalid data format specified')
  } else if (is.character(value)) {
    formats <- c('ascii', 'binary')
    ind <- charmatch(tolower(value), formats)
    if (is.na(ind))
      stop('Invalid data format specified')
    fi[2] <- ind
  } else 
    stop('New value must be an integer or a string.')
  attr(x, 'fileInfo') <- as.integer(fi)
  x
}

##' Timing information on AsspDataObj
##'
##' Some utility function to retrieve duration, number of records, sample rate and so on.
##' @title Timing information on AsspDataObj
##' @param x an object of class AsspDataObj
##' @return dur: the duration of the AsspDataObj in ms
##' @author Lasse Bombien
##' @export
##' @useDynLib wrassp
dur.AsspDataObj <- function(x) {
  if (!is.AsspDataObj(x))
    stop('Argument must be of class AsspDataObj.')
  numRecs.AsspDataObj(x) / attr(x, 'sampleRate')
}

##' @rdname dur.AsspDataObj
##' @return numRecs: the number of records stored in the AsspDataObj
##' @export
numRecs.AsspDataObj <- function(x) {
  attr(x, 'endRecord') - attr(x, 'startRecord') + 1
}

##' @rdname dur.AsspDataObj
##' @return rate: the data/sample rate of the AsspDataObj in Hz
##' @export
rate.AsspDataObj <- function(x) {
  attr(x, 'sampleRate')
}
