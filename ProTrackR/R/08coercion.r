.rowValid <- function(row)
{
  if (length(row) != 1) stop("only single integer value allowed for 'row'")
  if (!(row %in% 1:maximumPatternTableRowCount)) stop("Row not valid!")
}

.trackValid <- function(track)
{
  if (length(track) != 1) stop("only single integer value allowed for 'track'")
  if (!(track %in% 1:maximumTrackCount)) stop("Track not valid!")
}

.patternValid <- function(pattern)
{
  if (length(pattern) != 1) stop("only single integer value allowed for 'pattern'")
  if (!(pattern %in% 1:100)) stop("Pattern not valid!")
}

setGeneric("PTCell", function(x, row, track, pattern) standardGeneric("PTCell"))
setGeneric("PTCell<-", function(x, row, track, pattern, value) standardGeneric("PTCell<-"))

#' Coerce to or replace PTCell
#'
#' This method will coerce a set of objects to a \code{PTCell} object. It can also
#' be used to select specific cells from \code{PTModule},
#' \code{PTPattern} and \code{PTTrack} objects and replace
#' the selected \code{PTCell}.
#'
#' Method to coerce \code{x} to class \code{\link{PTCell}}.
#'
#' When \code{x} is \code{raw} data, it should consist of a \code{vector} of
#' 4 elements, formatted as specified in the \link{PTCell-class}.
#'
#' When \code{x} is a \code{character} string, it should be formatted as follows:
#' "\code{NNO SS EEE}", where \code{NN} is the note (for instance
#' "C-" or "A#", where the dash has no particular meaning and may be omitted,
#' the hash sign indicates a sharp note). Use a dash if the cell holds no note.
#' \code{O} is the octave (with a value of 0, or a dash, if a note
#' is missing, otherwise any of 1, 2 or 3). \code{SS} is the sample index
#' number, formatted as two hexadecimal digits (for example `1A'). \code{EEE} is
#' a three hexadecimal digit \code{\link{effect}} or trigger code (for more details see the
#' \code{\link{PTCell-class}}). The method is not case sensitive, so
#' you can use both upper and lower case. White spaces are ignored, you can use
#' as many as you would like. A correct \code{character} input for \code{x}
#' would be for example: "\code{A#2 01 A0F}". An `blank' \code{character}
#' representation would look like this: "\code{--- 00 000}".
#'
#' When \code{x} is of class \code{\link{PTTrack}}, \code{\link{PTPattern}}, or
#' \code{\link{PTModule}}, the \code{PTCell} at the specified indices (\code{row},
#' \code{track} and \code{pattern}) is returned, or can be replaced.
#'
#' @docType methods
#' @name PTCell-method
#' @rdname PTCell-method
#' @aliases PTCell,raw,missing,missing,missing-method
#' @param x Object (any of \code{raw} data, a \code{character} string, a \code{PTTrack},
#' a \code{PTPattern} or a \code{PTModule})
#' to coerce to a \code{\link{PTCell}}. See details below for the
#' required format of \code{x}.
#' @param row When \code{x} is a \code{PTTrack}, a \code{PTPattern},
#' or a \code{PTModule}, provide an index [1,64] of the row that needs
#' to be coerced to a \code{PTCell}.
#' @param track When \code{x} is a \code{PTPattern},
#' or a \code{PTModule}, provide an index [1,4] of the track that needs
#' to be coerced to a \code{PTCell}.
#' @param pattern When \code{x} is a \code{PTModule}, provide an index
#' of the pattern that needs to be coerced to a \code{PTCell}. Note that
#' ProTracker uses indices for patterns that start at zero, whereas R uses indices
#' that start at one. Hence add one to an index obtained from a \code{PTModule}
#' object (e.g., \code{x$pattern.order})
#' @param value An object of \code{\link{PTCell}} with which the \code{\link{PTCell}}
#' object at the specified indices in object \code{x} needs to be replaced.
#' @return When \code{PTCell} is used, a \code{PTCell} object
#' based on \code{x} is returned.
#'
#' When \code{PTCell<-} is used, object \code{x} is returned in which
#' the selected \code{PTCell} is replaced with \code{value}.
#' @examples
#' ## This will create an empty PTCell (equivalent
#' ## to new("PTCell"):
#' PTCell(raw(4))
#'
#' ## Use a character representation to create
#' ## a new PTCell object. A cell with note
#' ## B in octave 2, sample number 10 and with
#' ## effect '105':
#' cell <- PTCell("B-2 0A 105")
#'
#' data("mod.intro")
#'
#' ## replace PTCell at pattern number 1, track
#' ## number 2, and row number 3:
#' PTCell(mod.intro, 3, 2, 1) <- cell
#'
#' @author Pepijn de Vries
#' @family cell.operations
#' @export
setMethod("PTCell", c("raw", "missing", "missing", "missing"), function(x){
  x <- new("PTCell", data = x)
  return (x)
})

#' @export
#' @rdname PTCell-method
#' @aliases PTCell,character,missing,missing,missing-method
setMethod("PTCell", c("character", "missing", "missing", "missing"), function(x){
  #remove all white spaces:
  x <- gsub("[[:space:]]", "", as.character(x))
  if (any(nchar(x) < 7)) stop("x contains character strings that cannot be coerced to PTCells")
  x[nchar(x) < 8] <- paste(substr(x[nchar(x) < 8], 1, 1),
                           substr(x[nchar(x) < 8], 2, nchar(x[nchar(x) < 8])), sep = "-")
  if (any(nchar(x) < 7)) stop("x contains character strings that cannot be coerced to PTCells")
  x <- toupper(x)

  if (any(!(substr(x, 3, 3) %in% c("1", "2", "3", "-")))) stop ("invalid octave used in string")
  period <- noteToPeriod(substr(x, 1, 3))
  if (any(is.na(period))) stop ("invalid notes used in provided character string")

  cell_data <- as.list(as.integer(paste("0x", substr(x, 4, 4),
                                        sprintf("%03X", period),
                                        substr(x, 5, 8), sep = "")))

  x <- lapply(cell_data, function(x) new("PTCell", data = unsignedIntToRaw(x, 4)))
  if (length(x) == 1) x <- x[[1]]
  return (x)
})

#' @export
#' @rdname PTCell-method
#' @aliases PTCell,PTModule,numeric,numeric,numeric-method
setMethod("PTCell", c("PTModule", "numeric", "numeric", "numeric"), function(x, row, track, pattern){
  return(PTCell(PTPattern(x, pattern), row, track))
})

#' @export
#' @name PTCell<-
#' @rdname PTCell-method
#' @aliases PTCell<-,PTModule,numeric,numeric,numeric,PTCell-method
setReplaceMethod("PTCell", c("PTModule", "numeric", "numeric", "numeric", "PTCell"), function(x, row, track, pattern, value){
  .rowValid(row)
  .trackValid(track)
  .patternValid(pattern)
  x@patterns[[pattern]]@data[row, (-3 + track*4):(track*4)] <-
    value@data
  return(x)
})

#' @export
#' @rdname PTCell-method
#' @aliases PTCell,PTPattern,numeric,numeric,missing-method
setMethod("PTCell", c("PTPattern", "numeric", "numeric", "missing"), function(x, row, track){
  .rowValid(row)
  .trackValid(track)
  return(new("PTCell", data = x@data[row, (-3 + track*4):(track*4)]))
})

#' @export
#' @rdname PTCell-method
#' @aliases PTCell<-,PTPattern,numeric,numeric,missing,PTCell-method
setReplaceMethod("PTCell", c("PTPattern", "numeric", "numeric", "missing", "PTCell"), function(x, row, track, value){
  if (is.null(value)) stop("PTCell cannot be set to NULL")
  .rowValid(row)
  .trackValid(track)
  x@data[row, (-3 + track*4):(track*4)] <- value@data
  return(x)
})

#' @export
#' @rdname PTCell-method
#' @aliases PTCell,PTTrack,numeric,missing,missing-method
setMethod("PTCell", c("PTTrack", "numeric", "missing", "missing"), function(x, row){
  .rowValid(row)
  return(new("PTCell", data = x@data[row,]))
})

#' @export
#' @rdname PTCell-method
#' @aliases PTCell<-,PTTrack,numeric,missing,missing-method
setReplaceMethod("PTCell", c("PTTrack", "numeric", "missing", "missing", "PTCell"), function(x, row, value){
  if (is.null(value)) stop("not allowed to set a PTTrack to NULL")
  .rowValid(row)
  x@data[row,] <- value@data
  return(x)
})

setGeneric("PTTrack", function(x, track, pattern) standardGeneric("PTTrack"))
setGeneric("PTTrack<-", function(x, track, pattern, value) standardGeneric("PTTrack<-"))

#' Coerce to or replace PTTrack
#'
#' This method will coerce a set of objects to a \code{PTTrack} object. It can also
#' be used to select specific tracks from \code{PTModule} and
#' \code{PTPattern} objects and replace the selected \code{PTTrack}.
#'
#' Method to coerce \code{x} to class \code{\link{PTTrack}}.
#'
#' When \code{x} is a 64 by 4 \code{matrix} of \code{raw} data, each row
#' implicitly represents a \code{\link{PTCell}} object and should
#' be formatted accordingly. See \code{\link{PTCell-class}} documentation for
#' more details.
#'
#' When \code{x} is a 64 element \code{vector} of \code{character} representation
#' of \code{\link{PTCell}} objects, the \code{character} representation must be
#' conform the specifications as documented at the \code{\link{PTCell-class}}.
#'
#' When \code{x} is of class \code{\link{PTPattern}}, or
#' \code{\link{PTModule}}, the \code{PTTrack} at the specified indices
#' (\code{track} and \code{pattern}) is returned, or can be replaced.
#'
#' @docType methods
#' @name PTTrack-method
#' @rdname PTTrack-method
#' @aliases PTTrack,raw,missing,missing-method
#' @param x Object (any of \code{raw} data, a 64 by 4 \code{matrix} of \code{raw}
#' data, a vector of \code{character} strings,
#' a \code{PTPattern} or a \code{PTModule})
#' to coerce to a \code{\link{PTTrack}}. See details below for the
#' required format of \code{x}
#' @param track When \code{x} is a \code{PTPattern},
#' or a \code{PTModule}, provide an index [1,4] of the track that needs
#' to be coerced to a \code{PTTrack}.
#' @param pattern When \code{x} is a \code{PTModule}, provide an index
#' of the pattern that needs to be coerced to a \code{PTTrack}. Note that
#' ProTracker uses indices for patterns that start at zero, whereas R uses indices
#' that start at one. Hence add one to an index obtained from a \code{PTModule}
#' object (e.g., \code{x$pattern.order})
#' @param value An object of \code{\link{PTTrack}} with which the \code{\link{PTTrack}}
#' object at the specified indices in object \code{x} needs to be replaced.
#' @return When \code{PTTrack} is used, a \code{PTTrack} object
#' based on \code{x} is returned.
#'
#' When \code{PTTrack<-} is used, object \code{x} is returned in which
#' the selected \code{PTTrack} is replaced with \code{value}.
#' @examples
#' ## This will create an 'empty' PTTrack with all nul
#' ## values, which is equivalent to new("PTTrack"):
#' PTTrack(as.raw(0x00))
#'
#' ## This will generate a PTTrack from a repeated
#' ## character representation of a PTCell:
#' chan <- PTTrack(rep("C-3 01 C20", 64))
#'
#' data("mod.intro")
#'
#' ## This will replace the PTTrack at pattern
#' ## number 1, track number 2 of mod.intro with chan:
#' PTTrack(mod.intro, 2, 1) <- chan
#'
#' @family track.operations
#' @author Pepijn de Vries
#' @export
setMethod("PTTrack", c("raw", "missing", "missing"), function(x){
  return (new("PTTrack", data = matrix(x, ncol = 4, nrow = maximumPatternTableRowCount)))
})

#' @rdname PTTrack-method
#' @aliases PTTrack,matrix,missing,missing-method
#' @export
setMethod("PTTrack", c("matrix", "missing", "missing"), function(x){
  x <- apply(x, 2, as.raw)
  return (new("PTTrack", data = x))
})

#' @rdname PTTrack-method
#' @aliases PTTrack,character,missing,missing-method
#' @export
setMethod("PTTrack", c("character", "missing", "missing"), function(x){
  if (length(x) != maximumPatternTableRowCount) stop ("x should be a vector of 64 character strings")
  cells <- PTCell(x)
  raw.dat <- matrix(unlist(lapply(cells, as.raw)), ncol = 4, byrow = T)
  return (PTTrack(raw.dat))
})


#' @rdname PTTrack-method
#' @aliases PTTrack,PTModule,numeric,numeric-method
#' @export
setMethod("PTTrack", c("PTModule", "numeric", "numeric"), function(x, track, pattern){
  .trackValid(track)
  .patternValid(pattern)
  return (PTTrack(PTPattern(x, pattern), track))
})

#' @rdname PTTrack-method
#' @name PTTrack<-
#' @aliases PTTrack<-,PTModule,numeric,numeric,PTTrack-method
#' @export
setReplaceMethod("PTTrack", c("PTModule", "numeric", "numeric", "PTTrack"), function(x, track, pattern, value){
  .trackValid(track)
  .patternValid(pattern)
  x@patterns[[pattern]]@data[,(-3 + track*4):(track*4)] <- value@data
  return (x)
})

#' @rdname PTTrack-method
#' @aliases PTTrack,numeric,missing-method
#' @export
setMethod("PTTrack", c("PTPattern", "numeric", "missing"), function(x, track){
  .trackValid(track)
  return(new("PTTrack", data = x@data[,(-3 + track*4):(track*4)]))
})

#' @rdname PTTrack-method
#' @aliases PTTrack<-,numeric,missing,PTTrack-method
#' @export
setReplaceMethod("PTTrack", c("PTPattern", "numeric", "missing", "PTTrack"), function(x, track, value){
  if (is.null(value)) stop("not allowed to set a PTTrack to NULL.")
  .trackValid(track)
  x@data[,(-3 + track*4):(track*4)] <- value@data
  return(x)
})

setGeneric("PTPattern", function(x, pattern) standardGeneric("PTPattern"))
setGeneric("PTPattern<-", function(x, pattern, value) standardGeneric("PTPattern<-"))

#' Coerce to or replace PTPattern
#'
#' This method will coerce a set of objects to a \code{PTPattern} object. It can also
#' be used to select specific patterns from \code{PTModule} objects and replace
#' the selected \code{PTPattern}.
#'
#' Method to coerce \code{x} to class \code{\link{PTPattern}}.
#'
#' When \code{x} is a 64 by 16 \code{matrix} of \code{raw} data, each row
#' implicitly represents the \code{\link{PTCell}} objects of each of the
#' four tracks. Each \code{\link{PTCell}} consists of four \code{raw}
#' values. The values in each row are formatted accordingly, where the values of the
#' cells of each track are concatenated. See \code{\link{PTCell-class}} documentation for
#' more details on the \code{raw} format of a \code{\link{PTCell}} object.
#'
#' When \code{x} is a 64 by 16 \code{matrix} of \code{character} representations
#' of \code{\link{PTCell}} objects, the \code{character} representation must be
#' conform the specifications as documented at the \code{\link{PTCell-class}}.
#'
#' When \code{x} is of class \code{\link{PTModule}}, the \code{PTPattern} at the
#' specified index (\code{pattern}) is returned, or can be replaced.
#'
#' @docType methods
#' @name PTPattern-method
#' @rdname PTPattern-method
#' @aliases PTPattern,raw,missing-method
#' @param x Object (any of \code{raw} data, a 64 by 16 \code{matrix} of \code{raw} data, a 64 by 4 \code{matrix} of \code{character} strings,
#' or a \code{PTModule})
#' to coerce to a \code{\link{PTPattern}}. See details below for the
#' required format of \code{x}.
#' @param pattern When \code{x} is a \code{PTModule}, provide an index
#' of the pattern that needs to be coerced to a \code{PTPattern}. Note that
#' ProTracker uses indices for patterns that start at zero, whereas R uses indices
#' that start at one. Hence add one to an index obtained from a \code{PTModule}
#' object (e.g., \code{x$pattern.order}).
#' @param value An object of \code{\link{PTPattern}} with which the \code{\link{PTPattern}}
#' object at the specified \code{index} in object \code{x} needs to be replaced.
#' @return When \code{PTPattern} is used, a \code{PTPattern} object
#' based on \code{x} is returned.
#'
#' When \code{PTPattern<-} is used, object \code{x} is returned in which
#' the selected \code{PTPattern} is replaced with \code{value}.
#' @examples
#' ## This will create an 'empty' PTPattern with
#' ## all 0x00 values, which is equivalent to
#' ## new("PTPattern"):
#' PTPattern(as.raw(0x00))
#'
#' ## Create a PTPattern based on repeated
#' ## PTCell character representations:
#' pat <- PTPattern(matrix("F#2 1A 20A", 64, 4))
#'
#' data("mod.intro")
#'
#' ## Replace the first pattern in the patternOrder
#' ## table in mod.intro with 'pat' (don't forget to
#' ## add one (+1) to the index):
#' PTPattern(mod.intro,
#'           patternOrder(mod.intro)[1] + 1) <- pat
#' @family pattern.operations
#' @author Pepijn de Vries
#' @export
setMethod("PTPattern", c("raw", "missing"), function(x){
  return (new("PTPattern", data = matrix(x, ncol = 4*maximumTrackCount, nrow = maximumPatternTableRowCount)))
})

#' @rdname PTPattern-method
#' @aliases PTPattern,matrix,missing-method
#' @export
setMethod("PTPattern", c("matrix", "missing"), function(x){
  if (typeof(x) == "character")
  {
    if (any(dim(x) != c(maximumPatternTableRowCount, maximumTrackCount))) stop ("x should be a matrix of 64 by 4 holding character strings")
    cells <- PTCell(as.character(t(x)))
    raw.dat <- matrix(unlist(lapply(cells, as.raw)), ncol = 4*maximumTrackCount, byrow = T)
    return (PTPattern(raw.dat))
  }
  if (typeof(x) == "raw")
  {
    return (new("PTPattern", data = x))
  }
  stop("Matrix of x should be of type raw or character")
})

#' @rdname PTPattern-method
#' @aliases PTPattern,PTModule,numeric-method
#' @export
setMethod("PTPattern", c("PTModule", "numeric"), function(x, pattern){
  .patternValid(pattern)
  return (x@patterns[[pattern]])
})

#' @rdname PTPattern-method
#' @name PTPattern<-
#' @aliases PTPattern<-,PTModule,numeric,PTPattern-method
#' @export
setReplaceMethod("PTPattern", c("PTModule", "numeric", "PTPattern"), function(x, pattern, value){
  if (is.null(value)) stop ("PTPattern cannot be set to NULL! Use deletePattern to remove patterns.")
  .patternValid(pattern)
  if (is.null(x@patterns[[pattern]])) stop("Pattern does not yet exist and cannot be replaced. Use insertPattern to insert a pattern")
  x@patterns[[pattern]] <- value
  return(x)
})

setGeneric("PTSample", function(x, index) standardGeneric("PTSample"))
setGeneric("PTSample<-", function(x, index, value) standardGeneric("PTSample<-"))

#' Coerce to or replace PTSample
#'
#' This method will coerce a set of objects to a \code{PTSample} object. It can also
#' be used to select specific samples from \code{PTModule} objects and replace
#' the selected \code{PTSample}.
#'
#' Method to coerce \code{x} to class \code{\link{PTSample}}.
#'
#' When \code{x} is a \code{\link[tuneR]{Wave}} object, this method will not
#' resample it. However, the sample rate will be adjusted and samples exeeding
#' the maximum length of \code{2*0xffff} = \code{131070} will be clipped to this
#' maximum length. When \code{x} is a stereo sample, it will be converted to
#' mono, by averaging the left and right channel.
#'
#' When \code{x} is a \code{vector} of \code{raw} data, it will be truncated
#' if the maximum length of \code{2*0xffff} = \code{131070} is exceeded.
#' The raw will be converted with \code{\link{rawToSignedInt}} in order
#' to represent an 8 bit mono \code{\link{waveform}}.
#'
#' As samples must have an even length (as per ProTracker specifications),
#' a 0x00 value is appended if the length is odd.
#'
#' When \code{x} is of class \code{\link{PTModule}}, the \code{PTSample} at the
#' specified \code{index} is returned, or will be replaced.
#' @docType methods
#' @rdname PTSample-method
#' @name PTSample-method
#' @aliases PTSample,Wave,missing-method
#' @param x Object (any of class \code{\link[tuneR]{Wave}}, a \code{vector}
#' of \code{raw} data, or of class \code{\link{PTModule}}) that needs to
#' be coerced to a \code{\link{PTSample}} object. In the latter case, the
#' object can also be replaced.
#' @param index A positive \code{integer} index of the sample in \code{\link{PTModule}}
#' \code{x} that needs to be returned or replaced.
#' @param value An object of \code{\link{PTSample}} with which the \code{\link{PTSample}}
#' object at the specified \code{index} in object \code{x} needs to be replaced.
#' @return When \code{PTSample} is used, a \code{PTSample} object
#' based on \code{x} is returned.
#'
#' When \code{PTSample<-} is used, object \code{x} is returned in which
#' the selected \code{PTSample} is replaced with \code{value}.
#' @examples
#' ## Create a raw data sine wave:
#' raw_sine <- signedIntToRaw(round(sin(2*pi*(0:275)/276)*127))
#'
#' data("mod.intro")
#'
#' ## Replace sample number 1 from mod.intro
#' ## with the sine wave:
#' PTSample(mod.intro, 1) <-
#'   PTSample(raw_sine)
#'
#' ## Note that the replacement above
#' ## could also (maybe more efficiently)
#' ## be done with the 'waveform' method
#'
#' ## Restore the loop in sample number 1:
#' loopLength(PTSample(mod.intro, 1)) <- 276
#'
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("PTSample", c("Wave", "missing"), function(x){
  temp <- as.matrix(cbind(x@left, x@right))
  if (nrow(temp) > 2*0xffff)
  {
    warning("Wave object x is too long. Object length will be clipped to maximum length.")
    temp <- as.matrix(temp[1:(2*0xffff),])
  }
  delta <- integer(1)
  if(any(temp < 0)) delta <- as.integer(128)
  # length of a sample should be even:
  if ((nrow(temp) %% 2) == 1)
  {
    temp <- rbind(temp, rep(delta*(2^x@bit)/(256), 2))
  }
  x    <- new("PTSample",
              left = delta + as.integer(apply(temp, 1, function(x) mean(x, na.rm = T))/(2^(x@bit - 8))),
              wlooplen = as.raw(c(0, 1)))
  rm(temp)

  return(x)
})

#' @rdname PTSample-method
#' @aliases PTSample,raw,missing-method
#' @export
setMethod("PTSample", c("raw", "missing"), function(x){
  if (length(x) > 2*0xffff)
  {
    warning("Raw data too long. Clipped to maximum length")
    x <- x[1:(2*0xffff)]
  }
  # length should be even!
  if ((length(x) %% 2) == 1)
  {
    x <- c(x, raw(1))
  }
  return(new("PTSample",
             left = 128 + rawToSignedInt(x),
             wlooplen = as.raw(c(0, as.numeric(length(x) > 0)))))
})

#' @rdname PTSample-method
#' @aliases PTSample,PTModule,numeric-method
#' @export
setMethod("PTSample", c("PTModule", "numeric"), function(x, index){
  index <- as.integer(index)
  if (length(index) != 1) stop ("Index should have a length of 1.")
  if (!(index %in% 1:31)) stop ("Index out of range")
  return(x@samples[[index]])
})

#' @rdname PTSample-method
#' @name PTSample<-
#' @aliases PTSample<-,PTModule,numeric,PTSample-method
#' @export
setReplaceMethod("PTSample", c("PTModule", "numeric", "PTSample"), function(x, index, value){
  index <- as.integer(index)
  if (length(index) != 1) stop ("Index should have a length of 1.")
  if (!(index %in% 1:31)) stop ("Index out of range")
  x@samples[[index]] <- value
  return(x)
})
