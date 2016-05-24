#' Convert raw vectors into a character string
#'
#' A function that converts \code{raw} data into a \code{character} string.
#'
#' The function \code{\link{rawToChar}} will fail on vectors of \code{raw} data
#' with embedded \code{0x00} data. This function will not fail on embedded \code{0x00} values.
#' Instead, it will replace embedded \code{0x00} data with white spaces. Note that
#' leading and trailing \code{0x00} data will be omitted from the result.
#'
#' @param raw_dat A vector of class \code{raw} to be converted into a \code{character}.
#' @return A \code{character} string based on the \code{raw} data
#' @examples
#' ## generate some raw data with an embedded 0x00:
#' some.raw.data <- as.raw(c(0x68, 0x65, 0x6c, 0x6c, 0x6f, 0x00,
#'                           0x77, 0x6f, 0x72, 0x6c, 0x64, 0x21))
#' \dontrun{
#' ## this will fail:
#' try(rawToChar(some.raw.data))
#' }
#'
#' ## this will succeed:
#' rawToCharNull(some.raw.data)
#'
#' @family raw.operations
#' @family character.operations
#' @author Pepijn de Vries
#' @export
rawToCharNull <-
  function(raw_dat)
    ## function that converts vectors of raw data into character strings
    ## where null-characters are converted into spaces.
  {
    result <- ""
    if (length(raw_dat) < 3) try(result <- (rawToChar(raw_dat)), silent = T) else
    {
      result    <- raw_dat
      runlength <- rle(result)$lengths
      if (length(runlength) > 2)
      {
        rel_range <- (runlength[1] + 1):(length(result) - runlength[length(runlength)])
        result[rel_range][result[rel_range] == as.raw(0x00)] <- as.raw(0x20)
      }
      try(result <- rawToChar(result), silent = T)
      if (class(result) == "raw") result <- ""
    }
    return(result)
  }

#' Convert raw vector into a single unsigned integer value
#'
#' This function converts \code{raw} data into an unsigned integer
#'
#' This function converts a vector of raw data into a single unsigned integer.
#' for conversion of raw data into a vector of unsigned integers [0,255] use
#' `\code{\link{as.integer}(x)}'. For an inverse of this function
#' see \code{\link{unsignedIntToRaw}}.
#'
#' @param raw_dat A vector of class \code{raw} to be converted into an unsigned integer
#' @return A single unsigned integer value based on the provided \code{raw} data
#' @examples
#' ## generate some raw data:
#' some.raw.data <- as.raw(c(0x01, 0x1e, 0x3f))
#'
#' ## convert raw data into an unsigned integer:
#' rawToUnsignedInt(some.raw.data)
#'
#' ## note the difference with
#' as.integer(some.raw.data)
#'
#' @family raw.operations
#' @family integer.operations
#' @author Pepijn de Vries
#' @export
rawToUnsignedInt <-
  function(raw_dat)
  {
    if (class(raw_dat) != "raw") stop ("This function requires raw data as input")
    result <- as.integer(raw_dat)
    return(as.integer(sum(result*(256^((length(result):1) - 1)))))
  }

#' Convert unsigned integer into a raw vector
#'
#' This function converts an unsigned integer into a vector of \code{raw} data.
#'
#' This function converts an unsigned integer value into a vector (with
#' a specified length, namely \code{length.out}) of \code{raw} data. For the
#' inverse of this function use \code{\link{rawToUnsignedInt}(raw_dat)}
#'
#' @param int_dat A single integer value. If a list or vector of values.
#' is provided, only the first element is evaluated. Input data are converted
#' to absolute integer values.
#' @param length.out Required length of the vector that will hold the resulting.
#' \code{raw} data. Defaults to 1. If the value of \code{int_dat} is to large to convert into
#' \code{raw} data of length \code{length.out}, data will be clipped.
#' @return A vector of length \code{length.out}, holding \code{raw} data.
#' @examples
#' ## generate some unsigned integer:
#' some.integer <- 43251
#'
#' ## convert the unsigned integer into raw data:
#' unsignedIntToRaw(some.integer, length.out = 4)
#'
#' \dontrun{
#' ## note that the integer is too large to store as raw with length.out = 1:
#' unsignedIntToRaw(some.raw.data, length.out = 1)
#' }
#'
#' @family raw.operations
#' @family integer.operations
#' @author Pepijn de Vries
#' @export
unsignedIntToRaw <-
  function(int_dat, length.out = 1)
  {
    if (length(int_dat) > 1) warning ("Argument 'int_dat' has more than 1 element. Only first element converted")
    if (int_dat < 0) warning ("Argument 'int_dat' is signed, taking absolute value")
    int_dat <- abs(as.integer(int_dat[[1]]))
    if (int_dat >= 256^length.out)
    {
      int_dat <- 256^length.out - 1
      warning ("Argument 'in_dat' is out of range, proceeding with clipped value")
    }
    result  <- NULL
    remaining <- int_dat
    repeat
    {
      result <- c(remaining%%256, result)
      remaining <- floor(remaining/256)
      if(remaining <= 0) break
    }
    result <- as.raw(result)
    if (length(result) < length.out) result <- c(raw(length.out - length(result)), result)
    return(result)
  }

#' Convert signed integers (short) into a raw vector
#'
#' This function converts signed integer values into a vector of \code{raw} data.
#'
#' This function converts signed integer values [-128,127] into a vector of
#' \code{raw} data. The function
#' will fail on values that are out of range (< -128 or > 127). To convert
#' raw data into a vector of unsigned integers use \code{\link{as.integer}(x)}.
#' For the inverse of this function see \code{\link{rawToSignedInt}(raw_dat)}.
#'
#' @param int_dat A vector of integer values, ranging from -128 up to 127.
#' @return A vector of the same length as \code{int_dat}, holding \code{raw} data.
#' @examples
#' ## generate some signed integers:
#' some.integers <- c(-100, 40, 0, 30, -123)
#'
#' ## convert the signed integers into a vector of raw data:
#' signedIntToRaw(some.integers)
#'
#' @family raw.operations
#' @family integer.operations
#' @author Pepijn de Vries
#' @export
signedIntToRaw <-
  function(int_dat)
  {
    int_dat <- as.integer(int_dat)
    if(any(int_dat < -128)) stop("Some or all values out of range")
    if(any(int_dat >  127)) stop("Some or all values out of range")

    result <- int_dat
    result[result < 0] <- result[result < 0] + 256
    return(as.raw(result))
  }

#' Convert a raw vector into signed integers (short)
#'
#' This function converts a vector of \code{raw} data into signed integer values.
#'
#' This function converts a vector of \code{raw} data into signed integer values
#' [-128,127]. To convert unsigned integers into raw data use \code{\link{as.raw}(x)}.
#' For the inverse of this function see \code{\link{signedIntToRaw}(int_dat)}.
#'
#' @param raw_dat A vector of \code{raw} data.
#' @return A vector of the same length as \code{raw_dat}, holding signed integer values.
#' @examples
#' ## generate some raw data:
#' some.raw.data <- as.raw(c(0x68, 0x65, 0x6c, 0x6c, 0x6f, 0x90))
#'
#' ## convert the raw data into a vector of signed intgers:
#' rawToSignedInt(some.raw.data)
#'
#' @family raw.operations
#' @family integer.operations
#' @author Pepijn de Vries
#' @export
rawToSignedInt <-
  function(raw_dat)
    ## function that converts raw data into a signed (short) integers
  {
    if (class(raw_dat) != "raw") stop("Argument 'raw_dat' is not of class 'raw'")
    result <- as.integer(raw_dat)
    result[result > 127] <- result[result > 127] - 256
    return(as.integer(result))
  }

.getPeriodIndex <-
  function(period)
    ## use internally only to convert period numbers to note and octave
    ## get the row index of the period table based on a period value
    ## it doesn't work perfectly can rely on this 100%!
  {
    return(128 - round(logb(period, 1.060199)))
  }

#' Get the note and octave from period table
#'
#' These functions return the note and octave that is closest to the provided period value.
#'
#' ProTracker uses a \link{period_table} to link period values to certain
#' octaves and notes. This function serves to look up corresponding
#' notes and octaves for specific period values.
#'
#' @rdname periodToChar
#' @name periodToChar
#' @param period \code{integer} value of a period value.
#' @return \code{periodToChar} returns a \code{character} representing the combination
#' of octave and note that is closest to
#' \code{period} in the ProTracker period table.
#' @examples
#' ## Note C# in octave 3 is closest to a period of 200 in the table:
#' periodToChar(200)

#' ## try with a range of period values:
#' periodToChar(200:400)
#'
#' @family character.operations
#' @family period.operations
#' @family note.and.octave.operations
#' @author Pepijn de Vries
#' @export
periodToChar <-
  function(period)
    ## get note and octave character that corresponds with a specific period value
  {
    oct <- as.character(octave(period))
    oct[oct == "0"] <- "-"
    paste(note(period), oct, sep = "")
  }

#' Extract period value for a specific note
#'
#' Extracts the ProTracker period value for a specific note.
#'
#' ProTracker uses a \link{period_table} to link period values to certain
#' octaves and notes. This function serves to look up corresponding
#' period values for specific notes and octaves.
#'
#' @rdname noteToPeriod
#' @name noteToPeriod
#' @param note \code{character} string representing a note and octave for which the
#' ProTracker period value needs to be determined
#' @param finetune \code{integer} value ranging from -8 up to 7. A value used to
#' tune an audio sample.
#' @return Returns the \code{numeric} ProTracker period value for a corresponding
#' note, octave and \code{\link{fineTune}}. Returns 0 if a note could not be found in the
#' table.
#' @examples
#' ## Determine the period value corresponding with note 'A-3':
#' noteToPeriod("A-3")
#'
#' ## get the period values for notes 'A-3' and 'A#3' with finetune at -1:
#' noteToPeriod(c("A-3", "A#3"), -1)
#'
#' ## get the period values for note 'A-3' with finetune at 0 and 1:
#' noteToPeriod("A-3", 0:1)
#'
#' @family period.operations
#' @family note.and.octave.operations
#' @author Pepijn de Vries
#' @export
noteToPeriod <-
  function(note = "C-3", finetune = 0)
  {
    if (length(note) != length(finetune) &&
        (length(note) > 1 && length(finetune) > 1))
          stop("note and finetune should have the same length,
                or either should have length 1!")
    note <- toupper(as.character(note))
    note[nchar(note) == 2] <- paste(substr(note[nchar(note) == 2], 1, 1),
                                    substr(note[nchar(note) == 2], 2, 2), sep = "-")
    finetune <- as.integer((finetune))
    if (length(note) > length(finetune)) finetune <- rep(finetune, length(note))
    if (length(note) < length(finetune)) note <- rep(note, length(finetune))
    if (any(finetune > 7))
    {
      warning("finetune is out of range. Value clipped")
      fintune[finetune > 7] <- 7
    }
    if (any(finetune < -8))
    {
      warning("finetune is out of range. Value clipped")
      fintune[finetune < -8] <- 8
    }
    r_index <- suppressWarnings(apply(outer(ProTrackR::period_table$octave, as.numeric(substr(note, 3,3)), "==") &
                                      outer(ProTrackR::period_table$tuning, finetune, "=="),
                                      2, which))
    r_index <- unlist(lapply(as.list(r_index), function(x) ifelse(length(x) == 0, NA, x)))
    index <- cbind(r_index, match(substr(note, 1, 2), names(ProTrackR::period_table)))
    if (ncol(index) != 2) return (0)
    period <- ProTrackR::period_table[index]
    period[is.na(period)] <- 0
    return(period)
  }

#' Calculate the sample rate for a note or period value
#'
#' Calculate the sample rate for a note or a ProTracker period value.
#'
#' The timing on a Commodore Amiga depends on the video mode, which could be
#' either `\href{https://en.wikipedia.org/wiki/PAL}{PAL}'
#' or `\href{https://en.wikipedia.org/wiki/NTSC}{NTSC}'. Therefore sample
#' rates also depend on these modes. As the PAL is mostly used in Europe, and
#' the Amiga was most popular in Europe, PAL is used by default.
#'
#' @rdname sampleRate
#' @name sampleRate
#' @aliases periodToSampleRate
#' @aliases noteToSampleRate
#' @param period A ProTracker \code{integer} value of a period value for which the sample rate
#' is to be calculated.
#' @param note A \code{character} string representing a note for which the sample
#' rate is to be calculated.
#' @param finetune An \code{integer} value ranging from -8 up to 7. A value used to
#' tune an audio sample.
#' @param video The video mode used to calculate the sample rate. A \code{character}
#' string that can have either the value `\href{https://en.wikipedia.org/wiki/PAL}{PAL}'
#' or `\href{https://en.wikipedia.org/wiki/NTSC}{NTSC}'. PAL is used by default.
#' @return Returns the sample rate in samples per seconds.
#' @examples
#' ## calculate the sample rate for a ProTracker period value of 200
#' periodToSampleRate(200)
#'
#' ## calculate the sample rate for a sample at note 'A-3'
#' noteToSampleRate("A-3")
#'
#' ## note that the NTSC video system gives a slightly different rate:
#' noteToSampleRate("A-3", video = "NTSC")
#'
#' ## fine tuning a sample will also give a slightly different rate:
#' noteToSampleRate("A-3", finetune = -1)
#'
#' @family character.operations
#' @family period.operations
#' @family sample.rate.operations
#' @family note.and.octave.operations
#' @author Pepijn de Vries
#' @export
noteToSampleRate <-
  function(note = "C-3", finetune = 0, video = c("PAL", "NTSC"))
    ## get the sample rate that corresponds with a specific note
  {
    rate <- periodToSampleRate(noteToPeriod(note, finetune), match.arg(video))
    if (is.infinite(rate)) stop ("Note not in ProTracker period table.")
    return (rate)
  }

#' @rdname sampleRate
#' @name sampleRate
#' @export
periodToSampleRate <-
  function(period, video = c("PAL", "NTSC"))
    ## get the sample rate that corresponds with a specific period value
    ## based on the most common Amiga clock speed
  {
    ProTrackR::paula_clock$frequency[ProTrackR::paula_clock$video == match.arg(video)]/period
  }

#' Get the high or low nybble of a raw value
#'
#' Get the high or low nybble of a raw value and return as integer value [0,15].
#'
#' A \code{raw} is basically a byte, composed of 8 bits (zeros and ones).
#' A nybble is a 4 bit value. Hence, a raw value (or byte) is composed of
#' two nybbles. The leftmost nybble of a raw value is refered to as the
#' high nybble, the rightmost nybble is referred to as the low nybble.
#' These functions return either the high or low nybbles of raw data as integer
#' values [0,15].
#' As ProTracker stores some information as nybbles this function can be
#' used to retrieve this info.
#'
#' @rdname nybble
#' @name nybble
#' @param raw_dat A vector of class \code{raw} from which the high or low nybble value
#' needs to be extracted.
#' @param which A character string indicating whether the high or low nybble should
#' be returnd. It should either be "\code{low}" (default) or "\code{high}".
#' @return A vector of the same length as \code{raw_dat} holding integer values.
#' @examples
#' ## this will return 0x0f:
#' hiNybble(as.raw(0xf3))
#'
#' ## which is the same as:
#' nybble(as.raw(0xf3), "high")
#'
#' ## this will return 0x03:
#' loNybble(as.raw(0xf3))
#'
#' ## which is the same as:
#' nybble(as.raw(0xf3), "low")
#' @author Pepijn de Vries
#' @family nybble.functions
#' @family raw.operations
#' @family integer.operations
#' @export
nybble <-
  function(raw_dat, which = c("low", "high"))
  {
    if (match.arg(which) == "low") return  (loNybble(raw_dat))
    if (match.arg(which) == "high") return (hiNybble(raw_dat))
#    stop ("No valid argument 'which' provided. Should be 'low' or 'high'.")
  }

#' @rdname nybble
#' @export
loNybble <-
  function(raw_dat)
    ## function that gets the value [0,16] of the 4 low bits of a raw byte
  {
    if (class(raw_dat) != "raw") stop ("Only raw data is accepted as input")
    return(as.integer(raw_dat)%%16)
  }

#' @rdname nybble
#' @export
hiNybble <-
  function(raw_dat)
    ## function that gets the value [0,16] of the 4 high bits of a raw byte
  {
    if (class(raw_dat) != "raw") stop ("Only raw data is accepted as input")
    return(as.integer(as.integer(raw_dat)/16))
  }

#' Get signed integer values from nybbles
#'
#' Get signed integer values from one or more nybble.
#'
#' Nybbles are 4 bit values, where each byte (8 bits) holds two nybbles.
#' A high nybble (left-hand side of a byte) and a low nybble (right-hand
#' side of a byte). This function extracts a nybble from \code{raw} data
#' and converts it into a signed \code{integer} value ranging from -8 up to 7.
#' @rdname nybbleToSignedInt
#' @name nybbleToSignedInt
#' @param raw_dat \code{raw} data (either a single value or a \code{vector}),
#' from which a nybble will be extracted and converted.
#' @param which A \code{character} string indicating whether the "\code{low}" (default)
#' or "\code{high}" nybble of \code{raw_dat} needs to be converted into a signed
#' \code{integer}.
#' @return Returns \code{integer} values of the same length as \code{raw_dat},
#' ranging from -8 up to 7.
#' @examples
#' ## generate some raw data:
#'
#' rdat <- as.raw(255*runif(100))
#'
#' ## get signed integers of low nybbles:
#'
#' sintl <- nybbleToSignedInt(rdat)
#'
#' ## get signed integers of high nybbles:
#'
#' sinth <- nybbleToSignedInt(rdat, "high")
#'
#' @family nybble.functions
#' @family raw.operations
#' @family integer.operations
#' @author Pepijn de Vries
#' @export
nybbleToSignedInt <- function(raw_dat, which = c("low", "high"))
{
  nb <- nybble(raw_dat, which)
  nb[nb > 7] <- nb[nb > 7] - 16
  return(nb)
}

#' Convert a signed integer to a nybble in raw data.
#'
#' This function converts a signed integer ranging from -8 up to 7 into
#' either the high or low nybble of a byte, represented by \code{raw} data.
#'
#' Nybbles are 4 bit values, where each byte (8 bits) holds two nybbles.
#' A high nybble (left-hand side of a byte) and a low nybble (right-hand
#' side of a byte). This function converts a signed \code{integer} value
#' ranging from -8 up to 7 to a nybble and sets it as either a high or a low
#' nybble in \code{raw} data.
#' @rdname signedIntToNybble
#' @name signedIntToNybble
#' @param int_dat A single\code{intger} value or a \code{vector} of
#' \code{integer} data ranging from -8 up to 7.
#' @param which A character string indicating whether the nybble should
#' be set to the "\code{low}" (default) or "\code{high}" position of the
#' raw data that is returned.
#' @return Returns \code{raw} data of the same length as \code{int_dat}.
#' The returned raw data holds either low or high nybbles (as specified
#' by \code{which}) based on the provided signed \code{integer}s.
#' @examples
#' ## generate some integers in the right range:
#'
#' dati <- sample(-8:7, 100, replace = TRUE)
#'
#' ## Set the low nybbles of rawl based on dati:
#'
#' rawl <- signedIntToNybble(dati)
#'
#' ## Set the high nybbles of rawl based on dati:
#'
#' rawh <- signedIntToNybble(dati, "high")
#' @family nybble.functions
#' @family raw.operations
#' @family integer.operations
#' @author Pepijn de Vries
#' @export
signedIntToNybble <- function(int_dat, which = c("low", "high"))
{
  int_dat <- as.integer(int_dat)
  if (any(int_dat < -8) || any(int_dat > 7)) stop("int_dat out of range [-8,7]!")
  int_dat[int_dat < 0] <- int_dat[int_dat < 0] + 16
  fact <- ifelse(match.arg(which) == "low", 1, 16)
  return(as.raw(fact*int_dat))
}

#' Get the vibrato table used by ProTracker
#'
#' Gets the vibrato table as used by ProTracker in vibrato effects.
#'
#' As the old Commodore Amiga computer didn't have built-in mathematical functions,
#' many programs on that machine used their own data tables. As did ProTracker
#' for vibrato effects for which a sine function was used. As there was no sine
#' function that could be called, sine values were stored in a table.
#'
#' This function returns the \code{integer} sine values (ranging from 0 up
#' to 255) as a function of the table index (ranging from 0 up to 31).
#'
#' @rdname proTrackerVibrato
#' @name proTrackerVibrato
#' @param x \code{integer} representing the table index ranging from 0
#' up to 31. Values outside this range can be used, but will produce
#' results that are not valid in the context of ProTracker.
#' @return Returns an \code{integer} sine value ranging from 0 up to 255
#' when a valid table index (\code{x}) is provided. It will otherwise return
#' a sine value ranging from -255 up to 255.
#' @examples
#' ## this will return the table as used in ProTracker
#' proTrackerVibrato(0:31)
#' @author Pepijn de Vries
#' @export
proTrackerVibrato <- function(x)
{
  return(as.integer(255*sin(pi*(as.integer(x))/32)))
}

.unpackFibonacciDelta <-
  function(raw_data)
    ## function to decompress audio data (used by read.sample())
    ## based on following specs:
    ## http://amigadev.elowar.com/read/ADCD_2.1/Devices_Manual_guide/node02D6.html
  {
    if (class(raw_data) != "raw") stop ("Only raw data is accepted as input")
    fibonacci <- c(-34, -21, -13, -8, -5, -3, -2, -1, 0, 1, 2, 3, 5, 8, 13, 21)
    result    <- as.raw(rep(0, 2*length(raw_data) - 2))
    result[1] <- raw_data[2]
    for (i_byte in 2:length(result))
    {
      nybble     <- ifelse((i_byte%%2) == 0, c(hiNybble), c(loNybble))
      work_byte  <- raw_data[1 + (i_byte/2)]
      result_val <- rawToSignedInt(result[i_byte - 1]) + fibonacci[1 + nybble[[1]](work_byte)]
      ## Not entirely sure if I handle values out of range correctly here.
      ## I choose to cut values of, while it is also possible to use a modulus
      result_val <- ifelse(result_val >  127,  127, result_val)
      result_val <- ifelse(result_val < -128, -128, result_val)
      result[i_byte] <- signedIntToRaw(result_val)
    }
    return(result)
  }
