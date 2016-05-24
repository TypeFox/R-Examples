validity.PTSample <- function(object)
{
  if (length(object@name)       != 22)                 return(F)
  if (length(object@finetune)   != 1)                  return(F)
  if (length(object@volume)     != 1)                  return(F)
  if (length(object@wloopstart) != 2)                  return(F)
  if (length(object@wlooplen)   != 2)                  return(F)
  if (2*rawToUnsignedInt(object@wloopstart) > length(object@left))   return(F)
  if ((2*(rawToUnsignedInt(object@wlooplen) +
          rawToUnsignedInt(object@wloopstart)) > length(object@left)) &&
      rawToUnsignedInt(object@wlooplen) != 1)       return(F)
  # loop length can only be zero when the sample is empty
  if (rawToUnsignedInt(object@wlooplen) == 0 &&
      length(object@left) > 0)                         return(F)
  if (object@bit != 8)                                 return(F)
  if (!object@pcm)                                     return(F)
  if (as.integer(object@volume)   > 0x40)              return(F)
  if (as.integer(object@finetune) > 0x0F)              return(F)
  if (length(object@right) > 0)                        return(F)
  if (length(object@left)  > 2*0xFFFF)                 return(F)
  # sample length should be even!
  if ((length(object@left)%%2) == 1)                   return(F)
  return (T)
}

#' The PTSample class
#'
#' This class holds audio fragments with meta-information, to be used in
#' \code{\link{PTModule}} objects.
#'
#' This class holds audio fragments with meta-information (so-called samples),
#' to be used in \code{\link{PTModule}} objects. This class extends
#' the \code{\link[tuneR]{Wave}} class from \code{\link[tuneR]{tuneR}}. It therewith inherits
#' all properties and cool methods available from the \code{\link[tuneR]{tuneR}} package.
#' This allows you, for instance, to generate power spectra (\code{\link[tuneR]{powspec}})
#' of them. You can also plot the waveform with the \code{\link[tuneR]{plot-Wave}} method.
#' See \code{\link[tuneR]{tuneR}} for all possibilities with \code{\link[tuneR]{Wave}}
#' objects.
#' If you want you can also explicitly coerce \code{\link{PTSample}} to
#' \code{\link[tuneR]{Wave}} objects like this: \code{as(new("PTSample"), "Wave")}.
#'
#' The \code{\link{PTSample}} class has some slots that are additional to the
#' \code{\link[tuneR]{Wave}} class, as ProTracker requires additional information on
#' the sample with respect to its name, fine tune, volume and loop positions.
#' The \code{\link{PTSample}} class restricts the enherited \code{\link[tuneR]{Wave}}
#' class such that it will only hold 8 bit, mono, pcm waves with a maximum of
#' \code{2*0xffff = 131070} samples, as per ProTracker standards. The length should
#' always be even.
#'
#' \code{PTSample}s can be imported and exported using the
#' \code{\link{read.sample}} ans \code{\link{write.sample}} methods respectively.
#' \code{\link[tuneR]{Wave}} objects and \code{raw} data can be coerced to
#' \code{PTSample}s with the \code{\link{PTSample-method}}.
#'
#' @slot name A \code{vector} of length 22 of class "\code{raw}", representing
#' the name of the \code{PTSample}. It is often used to include
#' descriptive information in a \code{\link{PTModule}}. The name
#' of a sample can be extracted or replaced with the \code{\link{name}} method.
#' @slot finetune Single value of class "\code{raw}". The \code{\link{loNybble}}
#' of the \code{raw} value, represents the sample fine tune value ranging from -8 up to
#' 7. This value is used to tweak the playback sample rate, in order to tune it.
#' Negative values will lower the sample rate of notes, positive values will
#' increase the sample rate of notes. Period values corresponding to specific
#' notes and fine tune values are stored in the \code{\link{period_table}}.
#' The fine tune value can be extracted or replace with the \code{\link{fineTune}}
#' method.
#' @slot volume Single value of class "\code{raw}". The raw data corresponds with
#' the default playback volume of the sample. It ranges from 0 (silent) up to
#' 64 (maximum volume). The volume value can be extracted or replaced with the
#' \code{\link{volume}} method.
#' @slot wloopstart A \code{vector} of length 2 of class "\code{raw}". The \code{raw}
#' data represent a single unsigned number representing the starting position of
#' a loop in the sample. It should have a value of \code{0} when there is no loop.
#' Its value could range from \code{0} up to \code{0xffff}. To get the actual position
#' in bytes the value needs to be multiplied with 2 and can therefore only be
#' can only be even. The sum of the loop start position and the loop length should
#' not exceed the \code{\link{sampleLength}}. Its value can be extracted or
#' replaced with the \code{\link{loopStart}} method.
#' @slot wlooplen A \code{vector} of length 2 of class "\code{raw}". The \code{raw}
#' data represent a single unsigned number representing the length of
#' a loop in the sample. To get the actual length in bytes, this value needs to
#' be multiplied by 2 and can therefore only be even. It should have a value of
#' \code{2} when there is no loop.
#' Its value could range from \code{2} up to \code{2*0xffff} (= \code{131070}) and
#' can only be even (it can be \code{0} when the sample is empty). The sum of the
#' loop start position and the loop length should
#' not exceed the \code{\link{sampleLength}}. Its value can be extracted or
#' replaced with the \code{\link{loopLength}} method.
#' @slot left Object of class "\code{numeric}" representing the waveform of the
#' left channel. Should be \code{integer} values ranging from 0 up to 255.
#' It can be extracted or replaced with the \code{\link{waveform}} method.
#' @slot right Object of class "\code{numeric}" representing the right channel.
#' This slot is inherited from the \code{\link[tuneR]{Wave}} class and should be
#' \code{numeric(0)} for all \code{PTSample}s, as they are all mono.
#' @slot stereo Object of class "\code{logical}" whether this is a stereo representation.
#' This slot is inherited from the \code{\link[tuneR]{Wave}} class. As
#' \code{PTSample}s are always mono, this slot should have the value \code{FALSE}.
#' @slot samp.rate Object of class "\code{numeric}" representing the sampling rate.
#' @slot bit Object of class "\code{numeric}" representing the bit-wise quality.
#' This slot is inherited from the \code{\link[tuneR]{Wave}} class. As
#' \code{PTSample}s are always of 8 bit quality, the value of this slot
#' should always be 8.
#' @slot pcm Object of class "\code{logical}" indicating whether wave format is PCM.
#' This slot is inherited from the \code{\link[tuneR]{Wave}} class, for
#' \code{PTSample}s its value should always be \code{TRUE}.
#'
#' @name PTSample-class
#' @rdname PTSample-class
#' @aliases PTSample
#' @family sample.operations
#' @author Pepijn de Vries
#' @exportClass PTSample
setClass("PTSample",
         representation(name       = "raw",
                        finetune   = "raw",
                        volume     = "raw",
                        wloopstart = "raw",
                        wlooplen   = "raw"),
         prototype(name       = raw(22),
                   finetune   = raw(1),
                   volume     = as.raw(0x40),
                   wloopstart = raw(2),
                   wlooplen   = raw(2),
                   samp.rate  = 16574.28, #can't seem to be able to call a function (noteToSampleRate) from the constructor
                   bit        = 8,
                   stereo     = F),
         contains = "Wave",
         validity = validity.PTSample)

#' Resample data
#'
#' Resample \code{numeric} data to a different rate.
#'
#' This function resamples \code{numeric} data (i.e., audio data) from a
#' source sample rate to a target sample rate. At the core it uses
#' the \code{\link[stats]{approx}} function.
#' @rdname resample
#' @name resample
#' @param x A \code{numeric} \code{vector} that needs to be resampled.
#' @param source.rate The rate at which \code{x} was sampled in Hz (or
#' another unit, as long as it is in the same unit as \code{target.rate}).
#' @param target.rate The desired target sampling rate in Hz (or
#' another unit, as long as it is in the same unit as \code{source.rate}).
#' @param ... Arguments passed on to \code{\link[stats]{approx}}.
#' To simulate the Commodore Amiga hardware, it's best to
#' use '\code{method = "constant"} for resampling 8 bit samples.
#' @return Returns a resampled \code{numeric} \code{vector} of length
#' \code{round(length(x) * target.rate / source.rate)} based on \code{x}.
#' @examples
#' some.data <- 1:100
#'
#' ## assume that the current (sample) rate
#' ## of 'some.data' is 100, and we want to
#' ## resample this data to a rate of 200:
#' resamp.data <- resample(some.data, 100, 200, method = "constant")
#' @author Pepijn de Vries
#' @export
resample <- function(x, source.rate, target.rate, ...)
{
  x <- as.numeric(x)
  source.rate <- as.numeric(source.rate[[1]])
  target.rate <- as.numeric(target.rate[[1]])
  if (source.rate <= 0) stop ("Source rate should be greater than 1.")
  if (target.rate <= 0) stop ("Target rate should be greater than 1.")
  xout <- seq(1, length(x) + 1, length.out = round(length(x)*target.rate/source.rate))
  return(stats::approx(x, xout = xout, rule = 2, ...)[[2]])
}

setGeneric("fineTune", def = function(sample) standardGeneric("fineTune"))
setGeneric("fineTune<-", def = function(sample, value) standardGeneric("fineTune<-"))

#' Fine tune a PTSample
#'
#' Extract or replace the fine tune value of a \code{\link{PTSample}}.
#'
#' \code{\link{PTSample}}s can be tuned with their fine tune values.
#' The values range from -8 up to 7 and affect the playback sample rate of
#' specific notes (see \code{\link{period_table}}). This method can be used
#' to extract this value, or to safely replace it.
#'
#' @docType methods
#' @rdname fineTune
#' @name fineTune
#' @aliases fineTune,PTSample-method
#' @param sample A \code{\link{PTSample}} for which the fine tune value
#' needs to be extracted or replace.
#' @param value A \code{numeric} value ranging from -8 up to 7, representing
#' the fine tune.
#' @return For \code{fineTune} the fine tune value, represented by an
#' \code{integer} value ranging from -8 up to 7, is returned.
#'
#' For \code{fineTune<-} A \code{\link{PTSample}} \code{sample}, updated
#' with the fine tune \code{value}, is returned.
#' @examples
#' data("mod.intro")
#'
#' ## get the finetune of the first sample of mod.intro:
#'
#' fineTune(PTSample(mod.intro, 1))
#'
#' ## Let's tweak the finetune of the first sample of
#' ## mod.intro to -1:
#'
#' fineTune(PTSample(mod.intro, 1)) <- -1
#'
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("fineTune", "PTSample", function(sample){
  return(nybbleToSignedInt(sample@finetune))
})

#' @rdname fineTune
#' @name fineTune<-
#' @aliases fineTune<-,PTSample,numeric-method
#' @export
setReplaceMethod("fineTune", c("PTSample", "numeric"), function(sample, value){
  sample@finetune <- signedIntToNybble(value[[1]])
  return(sample)
})

setGeneric("volume", function(sample) standardGeneric("volume"))
setGeneric("volume<-", function(sample, value) standardGeneric("volume<-"))

#' Default playback volume of PTSample
#'
#' Extract or replace the default volume of a \code{\link{PTSample}}.
#'
#' \code{\link{PTSample}}s have a default playback volume, ranging from
#' \code{0} (silent) up to 64 (maximum volume). This method can be used
#' to extract this value, or to safely replace it.
#'
#' @docType methods
#' @rdname volume
#' @name volume
#' @aliases volume,PTSample-method
#' @param sample A \code{\link{PTSample}} for which the default volume
#' needs to be extracted or replace.
#' @param value A \code{numeric} value ranging from 0 up to 64, representing
#' the volume level.
#' @return For \code{volume} the volume value, represented by an
#' \code{integer} value ranging from 0 up to 64, is returned.
#'
#' For \code{volume<-} A \code{\link{PTSample}} \code{sample}, updated
#' with the volume \code{value}, is returned.
#' @examples
#' data("mod.intro")
#'
#' ## get the volume of the first sample of mod.intro:
#'
#' volume(PTSample(mod.intro, 1))
#'
#' ## Let's lower the volume of this sample to 32
#' ## (or as a hexadecimal: 0x20):
#'
#' volume(PTSample(mod.intro, 1)) <- 0x20
#'
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("volume", "PTSample", function(sample){
  return(rawToUnsignedInt(sample@volume))
})

#' @rdname volume
#' @name volume<-
#' @aliases volume<-,PTSample,numeric-method
#' @export
setReplaceMethod("volume", c("PTSample", "numeric"), function(sample, value){
  value <- as.integer(value[[1]])
  if (value < 0 || value > 63) stop("Volume out of range [0-64]!")
  sample@volume <- as.raw(value)
  return (sample)
})

setGeneric("loopStart", function(sample) standardGeneric("loopStart"))
setGeneric("loopStart<-", function(sample, value) standardGeneric("loopStart<-"))

#' The loop start position of a PTSample
#'
#' Extract or replace the loop start position of a \code{\link{PTSample}}.
#'
#' \code{\link{PTSample}}s can have loops, marked by a starting position
#' and length of the loop (in samples), for more details see the
#' \code{\link{PTSample-class}}. This method can be used to extract
#' the loop starting position or safely replace its value.
#'
#' @docType methods
#' @rdname loopStart
#' @name loopStart
#' @aliases loopStart,PTSample-method
#' @param sample A \code{\link{PTSample}} for which the loop start position
#' needs to be extracted or replace.
#' @param value An even \code{numeric} value giving the loop starting position in
#' samples ranging from 0 up to 131070. The sum of the \code{\link{loopStart}} and
#' \code{\link{loopLength}} should not exceed the \code{\link{sampleLength}}.
#'
#' Use a \code{value} of either \code{character} "\code{off}" or \code{logical}
#' "\code{FALSE}", in order to turn off the loop all together.
#' @return For \code{loopStart} the loop start position (in samples), represented by
#' an even \code{integer} value ranging from 0 up to 131070, is returned.
#'
#' For \code{loopStart<-} A \code{\link{PTSample}} \code{sample}, updated
#' with the loop start position `\code{value}', is returned.
#' @examples
#' data("mod.intro")
#'
#' ## get the loop start position of the
#' ## first sample of mod.intro:
#'
#' loopStart(PTSample(mod.intro, 1))
#'
#' ## Let's change the starting position of
#' ## the loop to 500
#'
#' loopStart(PTSample(mod.intro, 1)) <- 500
#'
#' ## Let's turn off the loop all together:
#'
#' loopStart(PTSample(mod.intro, 1)) <- FALSE
#'
#' @family sample.operations
#' @family loop.methods
#' @author Pepijn de Vries
#' @export
setMethod("loopStart", "PTSample", function(sample){
  return(2*rawToUnsignedInt(sample@wloopstart))
})

#' @rdname loopStart
#' @name loopStart<-
#' @aliases loopStart<-,PTSample-method
#' @export
setReplaceMethod("loopStart", c("PTSample", "ANY"), function(sample, value){
  value <- value[[1]]
  if (is.na(value) || value == "off" || value == F)
  {
    sample@wloopstart <- unsignedIntToRaw(0, 2)
    sample@wlooplen <- unsignedIntToRaw(1, 2)
  } else
  {
    value <- as.integer(round(value/2))
    if (value < 0 || value > (0xffff)) stop("Loop start out of range [0-(2*0xffff)]!")
    if ((value + rawToUnsignedInt(sample@wlooplen))*2 > length(sample@left)) stop("Loop start plus length is greater than sample length")
    sample@wloopstart <- unsignedIntToRaw(value, 2)
  }
  return (sample)
})

setGeneric("loopLength", function(sample) standardGeneric("loopLength"))
setGeneric("loopLength<-", function(sample, value) standardGeneric("loopLength<-"))

#' The loop length of a PTSample
#'
#' Extract or replace the loop length of a \code{\link{PTSample}}.
#'
#' \code{\link{PTSample}}s can have loops, marked by a starting position
#' and length of the loop (in samples), for more details see the
#' \code{\link{PTSample-class}}. This method can be used to extract
#' the loop length or safely replace its value.
#'
#' @docType methods
#' @rdname loopLength
#' @name loopLength
#' @aliases loopLength,PTSample-method
#' @param sample A \code{\link{PTSample}} for which the loop length
#' needs to be extracted or replace.
#' @param value An even \code{numeric} value giving the loop length in
#' samples ranging from 2 up to 131070 (It can be 0 when the sample is
#' empty). The sum of the \code{\link{loopStart}} and
#' \code{\link{loopLength}} should not exceed the \code{\link{sampleLength}}.
#'
#' Use a \code{value} of either \code{character} "\code{off}" or \code{logical}
#' "\code{FALSE}", in order to turn off the loop all together.
#' @return For \code{loopLength} the loop length (in samples), represented by
#' an even \code{integer} value ranging from 0 up to 131070, is returned.
#'
#' For \code{loopLength<-} A \code{\link{PTSample}} \code{sample}, updated
#' with the loop length `\code{value}', is returned.
#' @examples
#' data("mod.intro")
#'
#' ## get the loop length of the
#' ## first sample of mod.intro:
#'
#' loopLength(PTSample(mod.intro, 1))
#'
#' ## Let's change the length of
#' ## the loop to 200
#'
#' loopLength(PTSample(mod.intro, 1)) <- 200
#'
#' ## Let's turn off the loop all together:
#'
#' loopLength(PTSample(mod.intro, 1)) <- FALSE
#'
#' @family loop.methods
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("loopLength", "PTSample", function(sample){
  return(2*rawToUnsignedInt(sample@wlooplen))
})

#' @rdname loopLength
#' @name loopLength<-
#' @aliases loopLength<-,PTSample-method
#' @export
setReplaceMethod("loopLength", c("PTSample", "ANY"), function(sample, value){
  value <- value[[1]]
  value <- as.integer(round(value/2))
  if (is.na(value) || value == "off" || value == F)
  {
    sample@wloopstart <- unsignedIntToRaw(0, 2)
    sample@wlooplen <- unsignedIntToRaw(1, 2)
    if (length(sample@left) == 0) sample@wlooplen <- unsignedIntToRaw(0, 2)
  } else
  {
    if (value == 0 && length(sample) == 0)
      sample@wlooplen <- unsignedIntToRaw(0, 2)
    else
    {
      if (value < 1 || value > (0xffff)) stop("Loop length out of range [2 - (2*0xffff)]!")
      if ((value + rawToUnsignedInt(sample@wloopstart))*2 > length(sample@left)) stop("Loop start plus length is greater than sample length")
      sample@wlooplen <- unsignedIntToRaw(value, 2)
    }
  }
  return (sample)
})

setMethod("show", "PTSample", function(object){
  print(object)
})

#' @rdname print
#' @aliases print,PTSample-method
#' @export
setMethod("print", "PTSample", function(x, ...){
  cat("\nPTSample Object:\n")
  cat(paste("\tSample name:" , rawToCharNull(x@name), "\n", sep = "\t\t\t"))
  cat(paste("\tSample length (samples):", length(x@left), "\n", sep = "\t"))
  cat(paste("\tSample length (seconds):", length(x@left)/noteToSampleRate(), "\n", sep = "\t"))
  cat(paste("\tSample volume (0-64):", as.integer(x@volume), "\n", sep = "\t\t"))
  cat(paste("\tLoop start position:", 2*rawToUnsignedInt(x@wloopstart), "\n", sep = "\t\t"))
  cat(paste("\tLoop length:", 2*rawToUnsignedInt(x@wlooplen), "\n", sep = "\t\t\t"))
  cat(paste("\tFinetune:", fineTune(x), "\n", sep = "\t\t\t"))
})

setGeneric("playSample", function(x, silence = 0, wait = T,
                                  note = "C-3", loop = 1, ...){
  standardGeneric("playSample")
})

#' Play audio samples
#'
#' Method to play \code{\link{PTSample}}s or all such samples from
#' \code{\link{PTModule}} objects as audio.
#'
#' This method plays \code{\link{PTSample}}s and such samples from
#' \code{\link{PTModule}} objects, using the \code{\link[audio]{play}} method
#' from the audio package. Default \code{\link{fineTune}} and \code{\link{volume}}
#' as specified for the \code{\link{PTSample}} will be applied when playing
#' the sample.
#' @rdname playSample
#' @name playSample
#' @aliases playSample,PTSample-method
#' @param x Either a \code{\link{PTSample}} or a \code{\link{PTModule}} object.
#' In the latter case, all samples in the module will be played in order.
#' @param silence Especially for short samples, the \code{\link[audio]{play}} routine
#' can be a bit buggy: playing audible noise, ticks or parts from other samples at the end of the sample.
#' By adding silence after the sample, this problem is evaded. Use this argument
#' to specify the duration of this silence in seconds. When, \code{x} is a
#' \code{\link{PTModule}} object, the silence will also be inserted in
#' between samples.
#' @param wait A \code{logical} value. When set to \code{TRUE} the playing
#' routine will wait with executing any code until the playing is finished.
#' When set to \code{FALSE}, subsequent R code will be executed while playing.
#' @param note A \code{character} string specifying the note to be used for
#' calculating the playback sample rate (using \code{\link{noteToSampleRate}}).
#' It should start with the note (ranging from `A' up to `G') optionally followed
#' by a hash sign (`#') if a note is sharp (or a dash (`-') if it's not) and finally
#' the octave number (ranging from 1 up to 3). A valid notation would for instance
#' be `F#3'.
#' The \code{\link{fineTune}} as specified for the sample will also be used as
#' an argument for calculating the playback rate. A custom \code{finetune}
#' can also be passed as an argument to \code{\link{noteToSampleRate}}.
#' @param loop A positive \code{numeric} indicating the duration of a looped
#' sample in seconds. A looped sample will be played at least once, even if
#' the specified duration is less than the sum of \code{\link{loopStart}}
#' position and the \code{\link{loopLength}}.
#' See \code{\link{loopStart}} and \code{\link{loopLength}} for details on how
#' to set (or disable) a loop.
#' @param ... Further arguments passed on to \code{\link{noteToSampleRate}}.
#' Can be used to change the video mode, or finetune argument for the call to that method.
#' @return Returns nothing but plays the sample(s) as audio.
#' @examples
#' \dontrun{
#' data("mod.intro")
#'
#' ## play all samples in mod.intro:
#' playSample(mod.intro, 0.2, loop = 0.5)
#'
#' ## play a chromatic scale using sample number 3:
#' for (note in c("A-2", "A#2", "B-2", "C-3", "C#3",
#'                "D-3", "D#3", "E-3", "F-3", "F#3",
#'                "G-3", "G#3"))
#' {
#'   playSample(PTSample(mod.intro, 3), note = note, silence = 0.05, loop = 0.4)
#' }
#'
#' ## play the sample at a rate based on a specific
#' ## video mode and finetune:
#' playSample(PTSample(mod.intro, 3), video = "NTSC", finetune = -5)
#' }
#'
#' @author Pepijn de Vries
#' @family sample.operations
#' @family sample.rate.operations
#' @family play.audio.routines
#' @export
setMethod("playSample", "PTSample", function(x, silence, wait, note, loop, ...){
  finetune <- fineTune(x)
  vl <- volume(x)/0x40
  silence <- abs(as.numeric(silence[[1]]))
  wait <- as.logical(wait[[1]])
  note <- as.character(note[[1]])
  loop <- abs(as.numeric(loop[[1]]))
  if (loop == 0) stop ("'loop' should be greater than 0.")
  wf <- NULL
  if ("finetune" %in% names(list(...)))
    sr <- noteToSampleRate(note, ...)
  else
    sr <- noteToSampleRate(note, finetune, ...)
  if (loopState(x))
  {
    n_samp <- round(loop*sr)
    if (loopStart(x) + loopLength(x) > n_samp) n_samp <- loopStart(x) + loopLength(x)
    wf <- loopSample(x, n_samples = n_samp)
  }
  x <- as(x, "Wave")
  if (!is.null(wf)) x@left <- wf
  rm(wf)
  x@samp.rate <- sr
  if(silence > 0) x <- tuneR::bind(x,
                                   silence(silence, samp.rate = sr,
                                           bit = 8, pcm = T, xunit = "time"))
  if (wait)
  {
    audio::wait(audio::play(vl*(x@left - 128)/128,
                            rate = sr))
  } else
  {
    audio::play(vl*(x@left - 128)/128,
                rate = sr)
  }
  invisible()
})

setGeneric("read.sample", function(filename, what = c("wav", "mp3", "8svx", "raw")) standardGeneric("read.sample"))

#' Read an audio file and coerce to a PTSample object
#'
#' Reads audio files from "wav" and "mp3" files, using \code{\link[tuneR]{tuneR}}
#' methods. Commodore Amiga native formats "8svx" and "raw" can also be read.
#'
#' This method provides a wrapper for the \code{\link[tuneR]{readWave}} and
#' \code{\link[tuneR]{readMP3}} methods from \code{\link[tuneR]{tuneR}}. It also
#' provides the means to import audio from file formats native to the Commodore
#' Amiga. Simple \href{https://en.wikipedia.org/wiki/8SVX}{8svx} files (also known
#' as "iff" files) can be read, it has not yet been tested with more complex 8svx files.
#' It was also common practice to store audio samples as raw data on the
#' Commodore Amiga, where each byte simply represented a signed integer value
#' of the waveform.
#'
#' All audio will be coerced to 8 bit mono with a maximum length of
#' \code{2*0xffff} = {131070} bytes (= samples) as per ProTracker standards.
#' @rdname read.sample
#' @name read.sample
#' @aliases read.sample,character-method
#' @param filename A \code{character} string representing the filename to be read.
#' @param what A \code{character} string indicating what type of file is to be
#' read. Can be one of the following: "\code{wav}" (default), "\code{mp3}",
#' "\code{8svx}" or "\code{raw}".
#' @return Returns a \code{PTSample} object based on the file read.
#' @examples
#' \dontrun{
#' data("mod.intro")
#'
#' ## create an audio file which we can then read:
#' write.sample(PTSample(mod.intro, 2), "snaredrum.iff", "8svx")
#'
#' ## read the created sample:
#' snare <- read.sample("snaredrum.iff", "8svx")
#' print(snare)
#' }
#'
#' @note As per ProTracker standards, a sample should have an even length
#' (in bytes). If a sample file has an odd length, a \code{raw} \code{0x00} value
#' is added to the end.
#' @family sample.operations
#' @author Pepijn de Vries
#' @family io.operations
#' @export
setMethod("read.sample", c("character", "ANY"), function(filename, what = c("wav", "mp3", "8svx", "raw")){
  samp_name <- substr(basename(filename), 1, 22)
  if (match.arg(what) == "wav")
  {
    result <- tuneR::readWave(filename, from = 1, to = 2*0xffff)
    result <- PTSample(result)
    name(result) <- samp_name
    return(result)
  }
  if (match.arg(what) == "mp3")
  {
    result <- tuneR::readMP3(filename)
    result <- PTSample(result)
    name(result) <- samp_name
    return(result)
  }
  readRaw <- function (con)
  {
    result <- NULL
    repeat
    {
      l1 <- length(result)
      result <- c(result, readBin(con, "raw", 1024))
      l2 <- length(result)
      if ((l2 - l1) < 1024 || length(result) > 2*0xffff) break
    }
    if (length(result) > 2*0xffff)
    {
      warning("Sample is too long. It is clipped!")
      result <- result[1:(2*0xffff)]
    }
    # length should be even:
    if ((length(result) %% 2) == 1)
    {
      result <-  c(result, raw(0))
    }
    return(result)
  }
  VHDR_samp_rate       <- noteToSampleRate()
  if (match.arg(what) == "8svx")
  {
    # read Amiga IFF 8SVX (or raw) audio sample from file connection
    # and return as raw data
    # function based on file format specs provided here:
    # http://www.fileformat.info/format/iff/corion.htm
    # http://wiki.amigaos.net/wiki/8SVX_IFF_8-Bit_Sampled_Voice
    con <- file(filename, "rb")

    VHDR_compress        <- 0
    BODY_data            <- NULL

    FORM_id_char         <- readBin(con, "raw", 4)
    if (!all(FORM_id_char == charToRaw("FORM")))
    {
      warning ("Not an IFF file! Attempting to read as raw file")
      # reset connection to start of file:
      seek(con, 0)
      # keep reading chunks of 1 Kb of raw data untill end of file is reached
      BODY_data <- readRaw(con)
    } else
    {
      FORM_size_remainder  <- rawToUnsignedInt(readBin(con, "raw", 4))
      FORM_file_type       <- readBin(con, "raw", 4)
      if (!all(FORM_file_type == charToRaw("8SVX"))) stop ("File appears to be an iff but not an 8SVX file!")
      repeat
      {
        CHUNK_id_char      <- readBin(con, "raw", 4)
        CHUNK_size         <- rawToUnsignedInt(readBin(con, "raw", 4))
        if (all(CHUNK_id_char == charToRaw("BODY")))
        {
          BODY_data        <- readBin(con, "raw", CHUNK_size)
          break
        } else if (all(CHUNK_id_char == charToRaw("VHDR")))
        {
          VHDR_samp_high <- rawToUnsignedInt(readBin(con, "raw", 4))
          VHDR_samp_low  <- rawToUnsignedInt(readBin(con, "raw", 4))
          VHDR_samp_cyc  <- rawToUnsignedInt(readBin(con, "raw", 4))
          VHDR_samp_rate <- rawToUnsignedInt(readBin(con, "raw", 2))
          VHDR_num_oct   <- rawToUnsignedInt(readBin(con, "raw", 1))
          VHDR_compress  <- rawToUnsignedInt(readBin(con, "raw", 1))

          ############################################################
          # only delta-Fibonacci compression is currently supported
          ############################################################

          if (VHDR_compress > 1) stop("The compression type is currently not supported")
          if (VHDR_compress == 1) warning("Sample is compressed with a delta-Fibonacci algorithm. Decompression is implemented but not fully tested in this version of ProTrackR.")
          VHDR_volume    <- rawToUnsignedInt(readBin(con, "raw", 4))
        } else if (all(CHUNK_id_char == charToRaw("NAME")))
        {
          samp_name <- rawToCharNull(readBin(con, "raw", CHUNK_size))
        } else
        {
          # non relevant chunk, put in trashcan
          trashcan       <- readBin(con, "raw", CHUNK_size)
          rm(trashcan)
        }
      }
    }
    if (is.null(BODY_data)) stop("File BODY data is missing!")
    if (VHDR_compress == 1) BODY_data <- .unpackFibonacciDelta(BODY_data)
    close(con)
    result <- new("PTSample",
                  left = as.integer(rawToSignedInt(BODY_data) + 128),
                  samp.rate = VHDR_samp_rate,
                  wlooplen = as.raw(0:1))
    name(result) <- samp_name
    return (result)
  }
  if (match.arg(what) == "raw")
  {
    con <- file(filename, "rb")
    BODY_data <- readRaw(con)
    close(con)
    result <- new("PTSample",
                  left = as.integer(rawToSignedInt(BODY_data) + 128),
                  samp.rate = VHDR_samp_rate,
                  wlooplen = as.raw(0:1))
    name(result) <- samp_name
    return (result)
  }
})

setGeneric("write.sample", function(sample, filename, what = c("wav", "8svx", "raw")) standardGeneric("write.sample"))

#' Write a PTSample object to an audio file
#'
#' Write a \code{PTSample} as a "wav", "8svx" or "raw" audio file.
#'
#' This method provides a wrapper for the \code{\link[tuneR]{writeWave}} method
#' from \code{\link[tuneR]{tuneR}}. It also provides the means to export audio
#' to file formats native to the Commodore Amiga. \code{PTSample}s can be
#' exported as simple (uncompressed) \href{https://en.wikipedia.org/wiki/8SVX}{8svx}
#' files also known as "iff" files). In addition they can be exported as raw data,
#' where each byte simply represented a signed integer value of the waveform.
#'
#' @rdname write.sample
#' @name write.sample
#' @aliases write.sample,PTSample,character-method
#' @param sample A \code{PTSample} object that needs to be exported to an audio
#' file.
#' @param filename A \code{character} string representing the filename to which
#' the audio needs to be saved.
#' @param what A \code{character} string indicating what type of file is to be
#' exported. Can be one of the following: "\code{wav}" (default),
#' "\code{8svx}" or "\code{raw}".
#' @return Saves the audio to a file, but returns nothing.
#' @examples
#' \dontrun{
#' data("mod.intro")
#'
#' ## Export the second sample of mod.intro as a wav file:
#' write.sample(PTSample(mod.intro, 2), "snaredrum.wav", "wav")
#'
#' ## Export the second sample of mod.intro as an 8svx file:
#' write.sample(PTSample(mod.intro, 2), "snaredrum.iff", "8svx")
#'
#' ## Export the second sample of mod.intro as a raw file:
#' write.sample(PTSample(mod.intro, 2), "snaredrum.raw", "raw")
#' }
#'
#' @family sample.operations
#' @author Pepijn de Vries
#' @family io.operations
#' @export
setMethod("write.sample", c("PTSample", "character", "ANY"), function(sample, filename, what = c("wav", "8svx", "raw")){
  if (match.arg(what) == "wav")
  {
    tuneR::writeWave(sample, filename)
  }
  if (match.arg(what) == "8svx")
  {
    con <- file(filename, "wb")
    writeBin(charToRaw("FORM"), con)
    #hier moet de bestandsgrootte minus 8 byte (104 is grootte header):
    writeBin(unsignedIntToRaw(length(sample@left) + 104 - 8, 4), con)
    writeBin(charToRaw("8SVX"), con)
    writeBin(charToRaw("VHDR"), con)
    # VHDR size
    writeBin(unsignedIntToRaw(20, 4), con)
    # VHDR samp high
    writeBin(unsignedIntToRaw(length(sample@left), 4), con)
    # VHDR samp low
    writeBin(unsignedIntToRaw(0, 4), con)
    # VHDR samp cyc
    writeBin(unsignedIntToRaw(32, 4), con)
    # VHDR samp rate
    writeBin(unsignedIntToRaw(sample@samp.rate, 2), con)
    # VHDR samp oct
    writeBin(unsignedIntToRaw(1, 1), con)
    # VHDR samp compress
    writeBin(unsignedIntToRaw(0, 1), con)
    # VHDR samp volume
    writeBin(unsignedIntToRaw(0x0000ffff, 4), con)
    writeBin(charToRaw("ANNO"), con)
    writeBin(unsignedIntToRaw(20, 4), con)
    writeBin(c(charToRaw("ProTrackR"), raw(11)), con)
    writeBin(charToRaw("NAME"), con)
    writeBin(unsignedIntToRaw(20, 4), con)
    writeBin(sample@name[1:20], con)
    writeBin(charToRaw("BODY"), con)
    writeBin(unsignedIntToRaw(length(sample@left), 4), con)
    writeBin(signedIntToRaw(sample@left - 128), con)
    close(con)
  }
  if (match.arg(what) == "raw")
  {
    con <- file(filename, "wb")
    writeBin(signedIntToRaw(sample@left - 128), con)
    close(con)
  }
  invisible()
})

setGeneric("name", function(x) standardGeneric("name"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))

#' Obtain or replace the name of a PTModule or PTSample
#'
#' The name of both a \code{\link{PTModule}} and
#' \code{\link{PTSample}} are stored as \code{raw} data.
#' This method returns the name as a \code{character} string, or it can
#' be used to assign a new name to a \code{\link{PTModule}} or
#' \code{\link{PTSample}}.
#'
#' The name of a \code{\link{PTModule}} and
#' \code{\link{PTSample}} is stored as a \code{vector} of
#' \code{raw} data with a length of 20 or 22 respectively. This method
#' provides the means for getting the name as a \code{character} string
#' or to safely redefine the name of a \code{\link{PTModule}} or
#' \code{\link{PTSample}} object. To do so,
#' the provided name (\code{value}) is converted to a \code{raw} \code{vector}
#' of length 20 or 22 respectively. Long names may therefore get clipped.
#'
#' @docType methods
#' @rdname name
#' @name name
#' @aliases name,PTSample-method
#' @param x A \code{\link{PTModule}} or a \code{\link{PTSample}}
#' object for which to obtain or replace the name.
#' @param value A \code{character} string which should be used to replace the
#' name of \code{\link{PTModule}} or \code{\link{PTSample}} \code{x}.
#' @return For \code{name}, the name of the \code{\link{PTModule}} or
#' \code{\link{PTSample}} object as a \code{character} string is returned.
#'
#' For \code{name<-}, object \code{x} with an updated name is returned.
#' @examples
#' data("mod.intro")
#'
#' ## get the name of mod.intro:
#' name(mod.intro)
#'
#' ## I don't like the name, let's change it:
#' name(mod.intro) <- "I like this name better"
#'
#' ## Note that the provided name was too long and is truncated:
#' name(mod.intro)
#'
#' ## print all sample names in the module:
#' unlist(lapply(as.list(1:31), function(x)
#'   name(PTSample(mod.intro, x))))
#'
#' @family character.operations
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("name", "PTSample", function(x){
  return(rawToCharNull(x@name))
})

#' @rdname name
#' @name name<-
#' @aliases name<-,PTSample,character-method
#' @export
setReplaceMethod("name", c("PTSample", "character"), function(x, value){
  if (length(value) > 1) warning("Provided name has more than 1 element. Only first element used.")
  value <- as.character(value)[[1]]
  value <- charToRaw(value)
  if (length(value) > 22)
  {
    warning("Name is too long and will be truncated.")
    value <- value[1:22]
  }
  if (length(value) < 22) value <- c(value, raw(22 - length(value)))
  x@name <- value
  return (x)
})

setGeneric("sampleLength", function(sample) standardGeneric("sampleLength"))

#' Get the length of a PTSample
#'
#' Gets the length (in samples = bytes) of an audio fragment stored as a
#' \code{\link{PTSample}}.
#'
#' \code{\link{PTSample}}s are 8 bit mono audio fragments. This method
#' returns the length of this fragment expressed as number of samples (which
#' also equals the number of bytes).
#' @rdname sampleLength
#' @name sampleLength
#' @aliases sampleLength,PTSample-method
#' @param sample A \code{PTSample} object for which the length needs to be returned.
#' @return Returns a \code{numeric} value representing the number of samples
#' (bytes) the \code{PTSample} object \code{sample} is composed of.
#' @examples
#' data("mod.intro")
#'
#' ## Show the length of the second sample in mod.intro
#' sampleLength(PTSample(mod.intro, 2))
#'
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("sampleLength", "PTSample", function(sample){
  return(length(sample@left))
})

setGeneric("waveform", function(sample, start.pos = 1, stop.pos = sampleLength(sample), loop = TRUE) standardGeneric("waveform"))
setGeneric("waveform<-", function(sample, value) standardGeneric("waveform<-"))

#' Extract or replace a PTSample waveform
#'
#' Extract or replace the waveform of a \code{\link{PTSample}} object. The
#' waveform is represented by a \code{vector} of numeric values ranging from
#' 0 up to 255.
#'
#' Sample waveforms are stored as 8 bit signed short integer values ranging
#' from -128 up to +127 in original ProTracker files. However, as the
#' \code{\link{PTSample}} class extends the \code{\link[tuneR]{Wave}} class,
#' the waveforms are represented by integer values ranging from 0 up to 255
#' in the \link{ProTrackR} package. As per ProTracker specifications,
#' samples are of 8 bit mono quality and can only have an even length with
#' a maximum of \code{2*0xffff} = \code{131070}. This method can be used to
#' extract a waveform or replace it.
#' @rdname waveform
#' @name waveform
#' @aliases waveform,PTSample-method
#' @param sample A \code{\link{PTSample}} object from which the waveform needs to
#' be extracted or replaced.
#' @param start.pos A \code{numeric} starting index, giving the starting
#' position for the waveform to be returned. Default value is \code{1}. This
#' index should be greater than zero.
#' @param stop.pos A \code{numeric} stopping index, giving the stopping
#' position for the waveform to be returned. Default value is
#' \code{sampleLength(sample)} This index should be greater than
#' \code{start.pos}.
#' @param loop A \code{logical} value indicating whether the waveform
#' should be modulated between the specified loop positions
#' (see \code{\link{loopStart}} and \code{\link{loopLength}}),
#' or the waveform should stop at the end of the sample (padded with \code{NA}
#' values beyond the sample length). Will do the first
#' when set to \code{TRUE} and the latter when set to \code{FALSE}.
#' @param value A \code{vector} of numeric values ranging from 0 up to 255,
#' representing the waveform that should be used to replace that of object
#' \code{sample}. The length should be even and not exceed \code{2*0xffff} =
#' \code{131070}. \code{\link{loopStart}} and \code{\link{loopLength}} will
#' be adjusted automatically when they are out of range for the new waveform.
#'
#' Use \code{NA} to generate an empty/blank \code{\link{PTSample}} object.
#' @return For \code{waveform}, the waveform of \code{sample} is returned
#' as a \code{vector} of \code{numeric} values ranging from 0 up to 255.
#' If '\code{loop}' is set to \code{FALSE}
#' and the starting position is beyond the sample length, \code{NA} values
#' are returned. If '\code{loop}' is set to \code{TRUE} and the starting
#' position is beyond the sample loop (if present, see
#' \code{\link{loopState}}), the waveform is modulated between the loop
#' positions.
#'
#' For \code{waveform<-}, a copy of object \code{sample} is returned in which
#' the waveform has been replaced with \code{value}.
#' @examples
#' data("mod.intro")
#'
#' ## Loop sample #1 of mod.intro beyond it's
#' ## length of 1040 samples:
#' wav1 <- waveform(PTSample(mod.intro, 1),
#'                  1, 5000)
#'
#' ## get the waveform from sample #2
#' ## of mod.intro:
#' wav2 <- waveform(PTSample(mod.intro, 2))
#'
#' ## create an echo effect using
#' ## the extracted waveform:
#' wav2 <- c(wav2, rep(128, 1000)) +
#'         c(rep(128, 1000), wav2)*0.25 - 25
#'
#' ## assign this echoed sample to
#' ## sample #2 in mod.intro:
#' waveform(PTSample(mod.intro, 2)) <- wav2
#'
#' ## Blank out sample #1 in mod.intro:
#' waveform(PTSample(mod.intro, 1)) <- NA
#'
#' @family integer.operations
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("waveform", "PTSample", function(sample, start.pos, stop.pos, loop){
  start.pos  <- as.integer(abs(start.pos[[1]]))
  stop.pos   <- as.integer(abs(stop.pos[[1]]))
  if (start.pos < 1) stop("Starting position should be greater than or equal to 1...")
  if (start.pos > stop.pos) stop("Starting position should be greater than stopping position...")
  loop       <- as.logical(loop[[1]])
  samp_range <- start.pos:stop.pos
  if (loop && loopState(sample))
  {
    ls <- loopStart(sample)
    samp_range[samp_range > (ls + 1)] <-
      ((samp_range[samp_range > (ls + 1)] - (ls + 1)) %% loopLength(sample)) + ls + 1
  }
  return (as.integer(sample@left)[samp_range])
})

#' @rdname waveform
#' @name waveform<-
#' @aliases waveform<-,PTSample-method
#' @export
setReplaceMethod("waveform", c("PTSample", "ANY"), function(sample, value){
  value <- as.numeric(value)
  if (loopLength(sample) == 0 && length(value) > 0) loopLength(sample) <- 2
  if (any(is.na(value)) && length(value) > 1) stop ("NAs are not allowed in the data, if length > 1!")
  if (!any(is.na(value)) && (length(value)%%2) == 1)
  {
    warning("Length of data is odd. A value of 128 is added to the end.")
    value <- c(value, 128)
  }
  if (length(value) > 2*0xffff)
  {
    warning("Data exceeds maximum length (131070). Data will be truncated!")
    value <- value[1:(2*0xffff)]
  }
  if (!any(is.na(value)) && (any(value < 0) || any(value > 255)))
  {
    warning("Some values are out of range [0-255], data will be normalised to the required range")
    min_v <- min(value)
    max_v <- max(value)
    value <- as.integer(round(255*(value - min_v)/(max_v - min_v)))
  }
  if (loopStart(sample) > length(value))
  {
    warning("Sample loop start is outside the new range. It is set to 0.")
    sample@wloopstart <- raw(2)
  }
  if ((loopStart(sample) + loopLength(sample)) > length(value))
  {
    warning("Sample loop end is outside the new range. It's set to its maximum.")
    loopend <- as.integer((length(value) - loopStart(sample))/2)
    if (loopend == 0) loopend <- 1
    sample@wlooplen <- unsignedIntToRaw(loopend, 2)
  }
  if (any(is.na(value)) || length(value) == 0)
  {
    sample@left <- integer(0)
    if (loopStart(sample) == 0 && loopLength(sample) == 2) loopLength(sample) <- 0
  } else
  {
    sample@left <- value
  }
  return(sample)
})

setGeneric("loopSample", function(sample, times, n_samples) standardGeneric("loopSample"))

#' Looped waveform of a sample
#'
#' Generate a looped \code{\link{waveform}} of a \code{\link{PTSample}} object.
#'
#' For playing routines, it can be useful to generate repeats of a sample loop.
#' This method returns the waveform of a \code{\link{PTSample}} where the
#' loop is repeated `\code{times}' times or has a length of `\code{n_samples}'.
#' @rdname loopSample
#' @name loopSample
#' @aliases loopSample,PTSample-method
#' @param sample A \code{\link{PTSample}} object that needs to be looped.
#' @param times A positive \code{integer} value indicating the number of
#' times a sample loop should be repeated. This argument is ignore if
#' \code{n_samples} is specified.
#' @param n_samples A positive \code{integer} value indicating the desired length
#' of the looped waveform in number of samples. This argument overrules the
#' \code{times} argument.
#' @return Returns a \code{\link{waveform}} represented by a \code{numeric}
#' \code{vector} of values ranging from 0 up to 255. Has a length of
#' \code{n_samples} when that argument is specified.
#' @examples
#' data("mod.intro")
#'
#' ## Loop sample number 4 10 times:
#' wform <- loopSample(PTSample(mod.intro, 4), times = 10)
#' plot(wform, type = "l")
#'
#' ## Loop sample number 4, such that its
#' ## final length is 5000 samples:
#' wform <- loopSample(PTSample(mod.intro, 4), n_samples = 5000)
#' plot(wform, type = "l")
#'
#' @family loop.methods
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("loopSample", c("PTSample", "ANY", "ANY"), function(sample, times, n_samples){
  if (missing(times) && missing(n_samples)) stop ("Either 'times' or 'n_samples' should be specified.")
  if (!loopState(sample)) stop("No loop set to sample...")
  if (!missing(times))
  {
    times <- as.integer(abs(times[[1]]))
    n_samples <- loopStart(sample) + times*loopLength(sample)
  }
  if (!missing(n_samples)) n_samples <- as.integer(abs(n_samples[[1]]))
  return (waveform(sample, 1, n_samples))
})

setGeneric("loopState", function(sample) standardGeneric("loopState"))

#' Get PTSample loop state
#'
#' Determines whether a loop is specified for a \code{\link{PTSample}} object.
#'
#' The loop state is not explicitly stored in a \code{\link{PTSample}} object.
#' It can be derived from the \code{\link{loopStart}} position and
#' \code{\link{loopLength}}. This method is provided as a convenient method
#' to get the state. Use either \code{\link{loopStart}} or \code{\link{loopLength}}
#' to change the state.
#' @rdname loopState
#' @name loopState
#' @aliases loopState,PTSample-method
#' @param sample A \code{\link{PTSample}} object for which the loop state needs
#' to be determined.
#' @return Returns a \code{logical} value indicating whether a loop is (\code{TRUE})
#' or isn't (\code{FALSE}) specified for the \code{sample}.
#' @examples
#' data("mod.intro")
#'
#' ## Get the loop status of sample number 1
#' ## (it has a loop):
#' loopState(PTSample(mod.intro, 1))
#'
#' ## Get the loop status of sample number 2
#' ## (it has no loop):
#' loopState(PTSample(mod.intro, 2))
#' @family loop.methods
#' @family sample.operations
#' @author Pepijn de Vries
#' @export
setMethod("loopState", c("PTSample"), function(sample){
  return (!(loopLength(sample) <= 2 && loopStart(sample) == 0))
})
