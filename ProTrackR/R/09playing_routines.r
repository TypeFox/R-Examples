## http://pastebin.com/pg95YduC
## https://bel.fi/alankila/modguide/interpolate.txt
## XXX max output rate differs per channel?:
## http://eab.abime.net/showthread.php?t=70783

## The amiga hardware reference manual:
## https://archive.org/stream/Amiga_Hardware_Reference_Manual_1985_Commodore_a/Amiga_Hardware_Reference_Manual_1985_Commodore_a_djvu.txt
## The minimum period value you should use is 124
## ticks per sample NTSC (123 PAL) and the maximum is 65535.
## These limits apply to both PAL and NTSC machines.
## It doesn't say if these are hard limits, or whether these
## are just the specs, that you possible go beyond.
## in de WINUAE emulator lijkt de period value lager te kunnen dan de 124
## each channel different limit: http://eab.abime.net/showthread.php?t=70783
## WINUAE seem to use noteToPeriod("B-3") as the limit for all channels

## Olav Sorensen:
## Porta effects have period limits of 113 - 856 hard coded in ProTracker
## Other effect have no limits coded in ProTracker, the limit of 113 is also a 'hard' physical limit on a real Amiga
setGeneric("playMod", function(mod, wait = T, ...) standardGeneric("playMod"))

#' Play PTModule objects
#'
#' Converts \code{\link{PTModule}} objects into audio
#' \code{\link[tuneR]{Wave}}s, and plays them.
#'
#' Unfortunately, it was not feasible to create a routine that can directly
#' interpret \code{\link{PTModule}} objects and play them simultaneously.
#' Instead, the audio first needs to be rendered after which it can be played.
#' This method therefore first calls \code{\link{modToWave}} and then
#' \code{\link{playWave}}. Rendering may take some time and requires some
#' balance between speed, quality and accuracy. See the documentation of the
#' \code{\link{modToWave}} method for the control you have on these aspects.
#' @rdname playMod
#' @name  playMod
#' @aliases playMod,PTModule-method
#' @param mod A \code{\link{PTModule}} object to be played.
#' @param wait A \code{logical} value. When set to \code{TRUE} the playing
#' routine will wait with executing any code until the playing is
#' finished. When set to \code{FALSE}, subsequent R code will be
#' executed while playing.
#' @param ... Arguments that are passed on to \code{\link{modToWave}}.
#' @return A \code{\link[tuneR]{Wave}} object, generated from the
#' \code{mod} object, is returned.
#' @examples
#' \dontrun{
#' data("mod.intro")
#'
#' ## play the module and capture the audio Wave
#' wav <- playMod(mod.intro)
#' }
#' @author Pepijn de Vries
#' @family play.audio.routines
#' @family module.operations
#' @export
setMethod("playMod", "PTModule", function(mod, wait, ...){
  wait <- as.logical(wait[[1]])
  wav <- modToWave(mod, ...)
  playWave(wav, wait)
  return(wav)
})

setGeneric("playWave", function(wave, wait = T) standardGeneric("playWave"))

#' Play Wave objects
#'
#' Use the command line \code{\link[audio]{play}} function from the
#' \code{audio} package to play \code{\link[tuneR]{Wave}} objects.
#'
#' As the \code{\link{tuneR}} package play-function relies on external
#' players, this method is provided as a convenient approach to play
#' samples in the R console, using the \code{audio} package. Wave
#' objects are played at the rate as specified in the object.
#' @rdname playWave
#' @name playWave
#' @aliases playWave,Wave-method
#' @param wave An object of class \code{\link[tuneR]{Wave}}.
#' @param wait A logical value. When set to \code{TRUE} the playing
#' routine will wait with executing any code until the playing is
#' finished. When set to \code{FALSE}, subsequent R code will be
#' executed while playing.
#' @return Returns an \code{\link[audio]{$.audioInstance}}.
#' @examples
#' \dontrun{
#' data(mod.intro)
#'
#' ## PTSample objects can also be
#' ## played with this function as they
#' ## are a child of the Wave object:
#' playWave(PTSample(mod.intro, 2))
#' }
#' @family play.audio.routines
#' @author Pepijn de Vries
#' @export
setMethod("playWave", "Wave", function(wave, wait){
  wait <- as.logical(wait[[1]])
  result <- NULL
  play_expr <- function(x){
    if (length(x@right) == 0) y <- x@left else y <- rbind(x@left, x@right)
    #The 1.99 doesn't feel right. However, if we use 2, 'clicks'
    #will appear in the audio...
    audio::play(1.99*(y/((2^x@bit) - 1) +
                        ifelse(x@bit == 8, -0.5, 0)),
                rate = x@samp.rate)
  }
  if (wait) audio::wait(result <- play_expr(wave)) else result <- play_expr(wave)
  return (result)
})

setGeneric("modToWave",
           function(mod,
                    video             = c("PAL", "NTSC"),
                    target.rate       = 44100,
                    target.bit        = 16,
                    stereo.separation = 1,
                    low.pass.filter   = TRUE,
                    tracks            = 1:4,
                    ...){
             standardGeneric("modToWave")
           })

#' Convert a PTModule object into an audio Wave object
#'
#' Converts a \code{\link{PTModule}} object into a \code{\link[tuneR]{Wave}}
#' object, which can be played, further analysed, modified and saved.
#'
#' Before the \code{\link{PTModule}} object can be converted into a
#' \code{\link[tuneR]{Wave}} object, the rows of the
#' \code{\link{PTPattern}} objects in the module need to be put
#' in the right order. This method does that by calling
#' \code{\link{playingtable}}.
#'
#' Once the rows of the pattern tables are in the right order, all selected
#' \code{\link{PTTrack}} objects of the module are looped by this function
#' and the routines described below are applied to each track.
#'
#' On the Commodore Amiga the chip responsible for audio output (Paula),
#' the audio playback of samples can be controlled by the user in two ways:
#' the playback rate of the sample can be changed by specifying `period'
#' values (see e.g. \code{\link{periodToSampleRate}}) and specifying a
#' volume which is linearly scaled between 0 (silent) and 64 (maximum).
#'
#' So, for each track, the correct period and volume values are determined
#' based on the note, effect command and sample information in the module.
#'
#' Then, the \code{\link{PTSample}} objects are resampled, using the
#' period values and volume values as determined in the previous step.
#'
#' Next audio filters are applied to mimic original Commodore Amiga
#' sound. Finaly, the wave data for each separate track is mixed to
#' one (mono) or two (stereo) of the output channels.
#'
#' Converting ProTracker modules into wave objects can be time consuming.
#' The time required to convert an object obviously depends on your
#' machine's capacities and the length of the module but also the
#' complexity of the module. To speed up the conversion you could
#' reduce the target sample rate or turn off the low pass filter.
#' On modern machines, the time required for conversion should generally
#' be less than the playback time of the module.
#'
#' You can save the resulting \code{\link[tuneR]{Wave}} object by calling
#' \code{\link[tuneR]{writeWave}}.
#' @note As audio can be mixed with this package at frequencies much greater than the
#' Commodore Amiga's audio output rate, some aliasing of the sound could occur.
#' This results in high frequency audio, that would not be produced on an Amiga.
#' The current version of this package does not filter out these artefacts.
#' This should not be a problem if you're not concerned with producing an
#' accurate Amiga timbre.
#' @rdname modToWave
#' @name modToWave
#' @aliases modToWave,PTModule-method
#' @param mod An object of class \code{\link{PTModule}}
#' @param video The video mode of a Commodore Amiga affects timing routines and
#' the playback sample rate. This mode can be specified with this argument and
#' is represented by a \code{character} string that can have either the value
#' `\href{https://en.wikipedia.org/wiki/PAL}{PAL}'
#' or `\href{https://en.wikipedia.org/wiki/NTSC}{NTSC}'. PAL is used by default.
#' @param target.rate A positive \code{integer} sample rate for the target \code{Wave}.
#' Should be at least 2000.
#' Default value is 44100 Hz, which is conform CD quality. 22050 Hz will also
#' produce a decent sound quality and saves you some working memory.
#' @param target.bit Number of bits for the target \code{Wave}. Should be a \code{numeric}
#' value of either 8, 16, 24 or 32. Default is 16, which is conform CD quality
#' (the quality doesn't really improve at higher bit values, as the original
#' samples are of 8 bit quality).
#' @param stereo.separation A \code{numeric} value between 0 and 1.
#' When set to 1 (default), stereo channels (Amiga channels 1 and 4 on left,
#' and channels 2 and 3 on right) are completely separated. When set to
#' less than 1, stereo channels are mixed, where the number gives the
#' fraction of separation of the channels. When set to 0, both channels
#' are completely mixed and a mono \code{Wave} is returned.
#' @param low.pass.filter A \code{logical} value indicating whether low pass
#' filters should be applied when generating wave data. The Commodore Amiga
#' had hardware audio filters. One (low pass 6 db/Oct tuned at
#' 4.9 kHz) that filters all audio and one (low pass 12 db/Oct tuned at
#' approximately 3.3 kHz) that can be turned on and off at will with effect
#' command E00/E01 (see also \link{ProTrackR} documentation,
#' section on effect commands). These filters are only applied when the
#' \code{low.pass.filter} argument is set to \code{TRUE} and the
#' \code{target.rate} is set to values > 4.9 kHz. If you don't want to simulate
#' this typical Amiga sound, turn the filters off to save processing time.
#' @param tracks Either \code{logical} or \code{numeric} values indicating
#' which of the 4 \code{\link{PTTrack}}s are to be converted. By default
#' all 4 tracks are selected.
#' @param ... Additional arguments that are passed to \code{\link{playingtable}}.
#' @return A \code{\link[tuneR]{Wave}} object, generated from the
#' \code{mod} object is returned.
#' @examples
#' \dontrun{
#' data(mod.intro)
#' wav <- modToWave(mod.intro)
#' }
#' @author Pepijn de Vries
#' @family module.operations
#' @export
setMethod("modToWave", "PTModule", function(mod, video, target.rate, target.bit, stereo.separation, low.pass.filter, tracks, ...){
  video <- match.arg(video)
  verbose <- F
  tracks <- sort(unique((1:maximumTrackCount)[tracks]))
  moreArgs <- list(...)
  if ("verbose" %in% names(moreArgs)) verbose <- moreArgs[["verbose"]]
  target.rate <- as.integer(target.rate[[1]])
  if (target.rate <= 2000) stop("target.rate should be an integer of at least 2000...")
  target.bit <- as.integer(target.bit[[1]])
  if (!(target.bit %in% c(8, 16, 24, 32))) stop("target.bit should be any of the following values: 8, 16, 24, or 32.")
  low.pass.filter <- as.logical(low.pass.filter[[1]])
  stereo.separation <- as.numeric(stereo.separation[[1]])
  if (stereo.separation < 0 || stereo.separation > 1) stop("stereo.separation should be a value between 0 and 1.")

  pt       <- playingtable(mod, video = video, ...)
  if (!verbose) cat("generating channel wave data:\n")
  result   <- apply(t(as.matrix(tracks)), 2,
                    function(x) .generate.channel.data(mod, pt, x, target.rate,
                                                       video, low.pass.filter,
                                                       verbose))
  if (!verbose) cat("mixing channels...")
  ## filtering can result in values out of range (<0 or >255)
  ## let's normalise the data when nessecary:
  rmin <- min(result)
  if (rmin < 0) result <- 255*(result - rmin)/(255 - rmin)
  rmax <- max(result)
  if (rmax > 255) result <- 255*result/rmax
  result <- cbind(rowMeans(result[, tracks %in% (2:3), drop = F]),
                  rowMeans(result[, tracks %in% c(1,4), drop = F]))
  result <- (((2^(target.bit - 8)))*(256/255) - (1/255))*result
  result <- apply(result, 2, as.integer)
  if (target.bit > 8)
  {
    result <- result - as.integer(ceiling((2^target.bit)/2))
  }
  result <- Wave(left      = as.integer(result[,1]),
                 right     = as.integer(result[,2]),
                 bit       = target.bit,
                 samp.rate = target.rate)
  result@left[is.na(result@left)] <- ifelse(target.bit == 8, 127, 0)
  result@right[is.na(result@right)] <- ifelse(target.bit == 8, 127, 0)
  if (stereo.separation <= 0)
    result <- tuneR::mono(result)
  else if(stereo.separation < 1)
  {
    result <- Wave(left      = as.integer(result@left*(0.5 + 0.5*stereo.separation) +
                                            result@right*(0.5 - 0.5*stereo.separation)),
                   right     = as.integer(result@right*(0.5 + 0.5*stereo.separation) +
                                            result@left*(0.5 - 0.5*stereo.separation)),
                   bit       = target.bit,
                   samp.rate = target.rate)
  }
  if (!verbose) cat("\t\t\tdone\n")
  return(result)
})

setGeneric("playingtable",
           function(mod,
                    starting.position = 1,
                    max.duration      = 2*60,
                    speed             = 6,
                    tempo             = 0x7D,
                    video             = c("PAL", "NTSC"),
                    play.once         = T,
                    verbose           = F){
             standardGeneric("playingtable")
           })

#' Generate a table for playing a PTModule object
#'
#' This method generates a table (\code{data.frame}) in which information
#' from the pattern tables are put in the right order and in a comprehensive
#' format.
#'
#' This method generates a table (\code{data.frame}) in which information
#' from the pattern tables (\code{\link{PTPattern}}) are put in the right
#' order, taking into account pattern breaks, position jumps and pattern
#' loops (see also \link{ProTrackR} documentation,
#' section on effect commands). The information is put in
#' a comprehensive format in a \code{data.frame}, with the following columns:
#' \describe{
#' \item{row}{Row number index of the original
#' \code{\link{PTPattern}} object.}
#' \item{filter}{A \code{logical} value indicating whether the
#' Amiga hardware audio filter was either turned on or off using
#' effect command E00/E01 (see also \link{ProTrackR} documentation,
#' section on effect commands).}
#' \item{speed}{Number of `ticks' per row as set with the Fxy
#' effect commands in the module.}
#' \item{tempo}{The tempo as specified by the Fxy commands in the module.}
#' \item{delay}{The delay that should be applied to the row as
#' specified with the EEx effect command in the module.}
#' \item{effect.track1..4}{The effect code (\code{raw}) as specified in each of the
#' 4 tracks in the module.}
#' \item{effect.mag.track1..4}{The effect magnitude (\code{raw}) as specified
#' for each of the 4 tracks in the module.}
#' \item{sample.nr.track1..4}{The sample index number (\code{numeric}) as specified for
#' each of the 4 tracks in the module.}
#' \item{note.track1..4}{The note (\code{factor}) as specified for each of
#' the four tracks in the module.}
#' \item{position}{The positions index number (\code{numeric}) from the
#' \code{\link{patternOrder}} table in the module.}
#' \item{duration}{Playback duration of the corresponding row in seconds.}
#' \item{cum_duration}{Cumulative playback duration of the corresponding row in seconds.}
#' }
#' @rdname playingtable
#' @name playingtable
#' @aliases playingtable,PTModule-method
#' @param mod An object of class \code{\link{PTModule}}.
#' @param starting.position A \code{numeric} starting position index.
#' Determines where in the \code{\link{patternOrder}} table of the module to
#' start generating the playingtable.
#' @param max.duration A \code{numeric} value indicating the maximum
#' length in seconds of the pattern information returned. By default set
#' to 120 seconds (2 minutes). As some modules can be very long,
#' or contain infinite loops or position jumps, the maximum duration
#' is required to break out of the routine for generating the table.
#' @param speed Default speed to use when it is not specified in the
#' pattern data. See \link{ProTrackR} documentation for more info on
#' `speed' and `tempo'.
#' @param tempo Default tempo to use when it is not specified in the
#' pattern data.  See \link{ProTrackR} documentation for more info on
#' `speed' and `tempo'.
#' @param video The video mode of a Commodore Amiga affects timing routines.
#' This mode can be specified with this argument and
#' is represented by a \code{character} string that can have either the value
#' `\href{https://en.wikipedia.org/wiki/PAL}{PAL}'
#' or `\href{https://en.wikipedia.org/wiki/NTSC}{NTSC}'. PAL is used by default.
#' @param play.once A \code{logical} value. When set to \code{TRUE},
#' the routine will stop adding data to the table when the starting
#' position (\code{starting.position}) is reach once again. Warning:
#' may not work correctly when the last pattern contains a pattern
#' break. Will be overruled when the \code{maximum.duration} is reached
#' before the end of the song.
#' @param verbose A \code{logical} value. Suppresses a progress report
#' from being printed to the \code{\link[base]{sink}} when set to \code{TRUE}.
#' The default value is \code{FALSE}.
#' @return Returns a \code{data.frame} with pattern rows put in the right
#' order. Information contained in the returned table is described in the
#' 'Details' section
#' @examples
#' \dontrun{
#' data(mod.intro)
#' pt <- playingtable(mod.intro)
#' }
#' @author Pepijn de Vries
#' @family module.operations
#' @export
setMethod("playingtable", "PTModule", function(mod,
                                               starting.position,
                                               max.duration,
                                               speed,
                                               tempo,
                                               video,
                                               play.once,
                                               verbose){
  starting.position <- as.numeric(starting.position[[1]])
  if (!(starting.position %in% 1:patternOrderLength(mod)))
    stop("starting.position is out of the range of the patternOrder table of the module.")
  max.duration      <- as.numeric(max.duration[[1]])
  if (max.duration < 1)
    stop("max.duration should be at least 1 second.")
  speed             <- as.integer(speed[[1]])
  if (!(speed %in% 0x01:0x1F))
    stop("Default speed should range between 1 and 31.")
  tempo             <- as.integer(tempo[[1]])
  if (!(tempo %in% 0x20:0xFF))
    stop("Default speed should range between 32 and 255.")
  video             <- match.arg(video)
  play.once         <- as.logical(play.once[[1]])
  verbose           <- as.logical(verbose[[1]])

  remember.start.pos <- starting.position
  filter_fin <- F
  threshold  <- 0.1/44100
  # first create tables per pattern:
  if (!verbose) cat("Processing pattern tables...")

  pat_play_tables <- lapply(mod@patterns, function(pattern){
    result  <- data.frame(row = 1:maximumPatternTableRowCount)

    fx   <- rawShift(rawShift(pattern@data[,4*(1:maximumTrackCount) - 1], 4), -4)
    fx.m <- pattern@data[,4*(1:maximumTrackCount)]
    snr  <- matrix(as.integer(rawShift(rawShift(pattern@data[,4*(1:maximumTrackCount) - 3], -4), 4)) +
                     as.integer(rawShift(pattern@data[,4*(1:maximumTrackCount) - 1], -4)),
                   ncol = maximumTrackCount)
    note <- matrix(periodToChar(loNybble(as.vector(pattern@data[,4*(1:maximumTrackCount) - 3]))*256 +
                                  as.numeric(pattern@data[,4*(1:maximumTrackCount) - 2])),
                   ncol = maximumTrackCount)

    colnames(fx)   <- paste("effect.track",     1:maximumTrackCount, sep ="")
    colnames(fx.m) <- paste("effect.mag.track", 1:maximumTrackCount, sep ="")
    colnames(snr)  <- paste("sample.nr.track",  1:maximumTrackCount, sep ="")
    colnames(note) <- paste("note.track",       1:maximumTrackCount, sep ="")

    highest_track        <- function(x){
      x <- x[!is.na(x)]
      if (length(x) == 0) return (NA) else return(x[length(x)])
    }

    .specificEffectMagnitudes <- function(select.effect,
                                          mag.range,
                                          return.how = function(x) x,
                                          select.additional)
      # internal function to get specific effect magnitudes
    {
      selection  <- fx == select.effect &
        fx.m >= mag.range[1] &
        fx.m <= mag.range[2]
      if (!missing(select.additional)) selection <- selection & select.additional
      mags       <- matrix(NA, maximumPatternTableRowCount, maximumTrackCount)
      mags[selection] <- as.numeric(fx.m[selection])
      mags       <- apply(mags, 1, return.how)
      mags       <- mags
      mags[!apply(selection, 1, any)] <- NA
      return(mags)
    }

    ##############################################
    ## E0x Filter
    ## start
    ##############################################
    result$filter <- .specificEffectMagnitudes(as.raw(0x0E),
                                               as.raw(c(0x00, 0x0F)),
                                               highest_track)
    result$filter[!is.na(result$filter)] <-
      loNybble(as.raw(result$filter[!is.na(result$filter)]))%%2 == 0
    result$filter <- as.logical(result$filter)
    ##############################################
    ## E0x Filter
    ## end
    ##############################################

    ##############################################
    ## Fxy speed/tempo
    ## start
    ##############################################
    result$speed  <- .specificEffectMagnitudes(as.raw(0x0F),
                                               as.raw(c(0x00, 0x1F)),
                                               highest_track)
    result$tempo  <- .specificEffectMagnitudes(as.raw(0x0F),
                                               as.raw(c(0x20, 0xFF)),
                                               highest_track)
    ##############################################
    ## Fxy speed/tempo
    ## end
    ##############################################


    ##############################################
    ## E6x loop
    ## start
    ##############################################
    effect_sel    <- fx == as.raw(0x0E) & fx.m >= as.raw(0x60) & fx.m <= as.raw(0x6F)
    loop          <- matrix(as.integer(loNybble(as.vector(fx.m))),
                            ncol = maximumTrackCount)
    loop[!effect_sel] <- NA
    colnames(loop) <- paste("loop.track",       1:maximumTrackCount, sep ="")
    ##############################################
    ## E6x loop
    ## start
    ##############################################

    ##############################################
    ## EEx delay
    ## start
    ##############################################
    result$delay  <- .specificEffectMagnitudes(as.raw(0x0E),
                                               as.raw(c(0xE0, 0xEF)),
                                               highest_track)
    result$delay[is.na(result$delay)] <- 0
    result$delay  <- loNybble(as.raw(result$delay))
    ##############################################
    ## EEx delay
    ## end
    ##############################################

    ##############################################
    ## Bxy position jump
    ## start
    ##############################################
    # Bxy xy 00-99 (as hexadeximals) jump to position xy + 1
    # If multiple tracks have a specified a break at the same row,
    # the highest track has priority, remainder is ignored.
    # If the specified position is out of range, jump to position 1.
    result$position.jump <- .specificEffectMagnitudes(as.raw(0x0B),
                                                      as.raw(c(0x00, 0xFF)),
                                                      highest_track)
    ##############################################
    ## Bxy position jump
    ## end
    ##############################################

    ##############################################
    ## Dxy pattern break
    ## start
    ##############################################
    # Dxy xy 00-63 (as decimals) jump to row xy + 1 in next pattern.
    # if xy is out of range, will jump to row 1 of next pattern.
    # When x is in the hex range ('a' - 'f') the break will jump to row 1.
    # when y is in the hex range ('a' - 'f') its value will be added
    # to the decimal value of x0. e.g. 1f will jump to row 10 + 15 + 1 = 26 of the
    # next pattern. If multiple tracks have a specified a break at the same row,
    # the highest track has priority, remainder is ignored.
    # Dxy effect left of Bxy effects should be ignored.
    right_of_Bxy <- t(apply(fx, 1, function(x){
      result <- rep(T, length(x))
      bright <- suppressWarnings(max(which(x == as.raw(0x0b))))
      if (!is.infinite(bright)) result[1:bright] <- F
      return(result)
    }))
    result$pattern.break   <- .specificEffectMagnitudes(as.raw(0x0D),
                                                        as.raw(c(0x00, 0xFF)),
                                                        highest_track,
                                                        right_of_Bxy)
    selection <- !is.na(result$pattern.break)
    result$pattern.breaklo[selection] <- loNybble(as.raw(result$pattern.break[selection]))
    result$pattern.breakhi[selection] <- hiNybble(as.raw(result$pattern.break[selection]))
    result$pattern.break[selection] <-
      result$pattern.breakhi[selection]*10 +
      result$pattern.breaklo[selection]

    result$pattern.break[which(result$pattern.break > 63)] <- 0
    result$pattern.breaklo <- NULL
    result$pattern.breakhi <- NULL

    ## if a pattern break is at the same row as a row delay command (EEx),
    ## the break will actually start a row later.
    effect_sel <- with(result, !is.na(pattern.break) &
                         apply(fx == as.raw(0x0e) &
                                 fx.m >= as.raw(0xe0) & fx.m <= as.raw(0xef), 1, any))
    result$pattern.break[effect_sel] <- result$pattern.break[effect_sel] + 1
    ##############################################
    ## Dxy pattern break
    ## end
    ##############################################

    result <- cbind(result, fx, fx.m, snr, note, loop)

    return(result)
  })
  if (!verbose)
  {
    cat("\t\tdone\n")
    cat("Processing pattern order...")
  }

  # Put the tables of each pattern in the right order:
  result          <- data.frame()
  row.skip        <- 0
  cum_duration    <- 0
  vblank_duration <- ifelse(match.arg(video) == "PAL", 1/50, 1/60)
  # flag indicating for each channel whether there is a pattern loop active:
  loop.on         <- rep(F, maximumTrackCount)
  # counter for pattern loops in each channel. Will count down to zero when
  # a loop is active, and will continue playing once it reaches zero.
  # when multiple loop ends are positioned in the same row, all
  # counters need to reach zero
  loop.counter    <- rep(0, maximumTrackCount)
  # starting position for each loop in each channel. Zero by default or
  # otherwise specified with effect command E60.
  loop.start.pos  <- rep(1, maximumTrackCount)
  while (TRUE)
  {
    pat_tab          <- pat_play_tables[[1 + patternOrder(mod)[starting.position]]]
    pat_tab$position <- starting.position
    pat_tab$row      <- 1:maximumPatternTableRowCount

    ## look for loop-starts (E60) from top to bottom using a while loop.
    loops                <- pat_tab[,grepl("loop.track", names(pat_tab), fixed = T)]
    loop.pos             <- which(apply(loops, 1, function(x) any(!is.na(x))))
    pos                  <- 1 + row.skip
    loop.index           <- ifelse (length(loop.pos) > 0,
                                    min(which(loop.pos > row.skip)),
                                    1)
    pat_tab$duration     <- NA
    pat_tab$cum_duration <- NA
    pat_tab_ext <- data.frame()
    if (length(loop.pos) > 0)
    {
      while (TRUE)
      {
        pat_tab_ext <- rbind(pat_tab_ext, pat_tab[pos:loop.pos[loop.index],])
        row.skip <- 0
        if (any(!is.na(pat_tab_ext$position.jump)) ||
            any(!is.na(pat_tab_ext$pattern.break))) break
        if (any(loops[loop.pos[loop.index],] > 0, na.rm = T))
        {
          ## if loop is currently off, set the counter...
          loop.counter[which(!loop.on & loops[loop.pos[loop.index],] > 0)] <-
            unlist(loops[loop.pos[loop.index], which(!loop.on & loops[loop.pos[loop.index],] > 0)]) + 1

          ## ... and turn the loop on
          loop.on[which(!loop.on & loops[loop.pos[loop.index],] > 0)] <- T

          ## When the loop is turned on, subtract 1 from the counter
          loop.counter[which(loop.on & loops[loop.pos[loop.index],] > 0)]  <-
            loop.counter[which(loop.on & loops[loop.pos[loop.index],] > 0)] - 1
          ## If the counter reaches 0, turn off the loop
          loop.on[loop.on & loop.counter == 0] <- F

          pos <- suppressWarnings(max(loop.start.pos[which(loop.on & loops[loop.pos[loop.index],] > 0)]))
          if (is.infinite(pos) || all(!loop.on))
          {
            pos <- loop.pos[loop.index] + 1
            loop.index <- loop.index + 1
          } else
          {
            loop.index <- min(which(loop.pos > pos))
          }
        } else if (any(loops[loop.pos[loop.index],] == 0, na.rm = T))
        {
          loop.start.pos[which(loops[loop.pos[loop.index],] == 0)] <- loop.pos[loop.index]
          pos <- loop.pos[loop.index] + 1
          loop.index <- loop.index + 1
        }
        if (loop.index > length(loop.pos))
        {
          if (loop.pos[loop.index - 1] < maximumPatternTableRowCount)
            pat_tab_ext <- rbind(pat_tab_ext,
                                 pat_tab[(loop.pos[loop.index - 1] + 1):maximumPatternTableRowCount,])
          break
        }
        #if E60 -> set loopstart
        #if E6x -> check if loop is already active. if so, counter - 1;
        #          if not, turn loop on and set the counter.
        pat_tab_ext <- .fill.pattern.table(pat_tab_ext, cum_duration, speed, tempo, vblank_duration, filter_fin)
        cum_duration <- pat_tab_ext$cum_duration[length(pat_tab_ext$cum_duration)]
        if ((cum_duration - max.duration) > threshold) break
      }
      pat_tab <- pat_tab_ext
      rm(pat_tab_ext)
    } else
    {
      pat_tab <- pat_tab[pos:maximumPatternTableRowCount,]
      row.skip <- 0
    }
    pat_tab$loop.track1 <- NULL
    pat_tab$loop.track2 <- NULL
    pat_tab$loop.track3 <- NULL
    pat_tab$loop.track4 <- NULL

    # Dxy or Bxy
    # find pattern breaks or position jumps and take action
    p.jump           <- suppressWarnings(min(which(!is.na(pat_tab$position.jump))))
    p.jump[is.infinite(p.jump)] <- 999
    p.brk            <- suppressWarnings(min(which(!is.na(pat_tab$pattern.break))))
    p.brk[is.infinite(p.brk)]   <- 999
    if (p.jump <= p.brk && p.jump != 999)
    {
      starting.position <- pat_tab$position.jump[!is.na(pat_tab$position.jump)][1]
      pat_tab           <- pat_tab[1:p.jump,]
    }
    if (p.brk <= p.jump && p.brk != 999)
    {
      row.skip <- pat_tab$pattern.break[!is.na(pat_tab$pattern.break)][1]
      ## when pattern break was used in combination with EEx, the break position
      ## can beyond the end of the next pattern, hence the following wrap:
      if (row.skip >= maximumPatternTableRowCount)
      {
        row.skip <- row.skip - maximumPatternTableRowCount
        starting.position <- starting.position + 1
        if (starting.position > patternOrderLength(mod)) starting.position <- 1
      }
      pat_tab           <- pat_tab[1:p.brk,]
    }

    pat_tab      <- .fill.pattern.table(pat_tab, cum_duration, speed, tempo, vblank_duration, filter_fin)
    cum_duration <- pat_tab$cum_duration[length(pat_tab$cum_duration)]
    result       <- rbind(result, pat_tab)
    speed0 <- any(result$speed == 0)
    if (speed0)
    {
      result <- result[1:min(which(result$speed == 0)),]
      final_speed <- result$speed[result$speed != 0]
      final_speed <- final_speed[length(final_speed)]
      result$speed[result$speed == 0] <- final_speed
      break
    }
    cum_duration <- result$cum_duration[nrow(result)]
    if ((cum_duration - max.duration) > threshold)
    {
      result <- result[(result$cum_duration - max.duration) <= threshold,]
      break
    }
    if (speed0) break
    starting.position <- starting.position + 1
    if (starting.position > patternOrderLength(mod)) starting.position <- 1
    speed <- result$speed[nrow(result)]
    tempo <- result$tempo[nrow(result)]
    filter_fin <- result$filter[nrow(result)]
    ## If the new starting position equals the old one, stop playing
    if (starting.position == remember.start.pos && play.once) break
  }
  # pattern.breaks and position jumps have been processed and is
  # no longer relevant. These columns are removed.
  result$pattern.break <- NULL
  result$position.jump <- NULL
  noteAsFact <- function(x) factor(as.character(x),
                                   levels = c("---", as.vector(outer(names(ProTrackR::period_table)[-1:-2],
                                                                     1:3, function(x, y) paste(x, y, sep = "")))))
  result$note.track1 <- noteAsFact(result$note.track1)
  result$note.track2 <- noteAsFact(result$note.track2)
  result$note.track3 <- noteAsFact(result$note.track3)
  result$note.track4 <- noteAsFact(result$note.track4)
  if (!verbose) cat("\t\tdone\n")
  return(result)
})

.fill.pattern.table <- function(pat_tab, cum_duration, speed, tempo, vblank_duration, filter_fin)
{
  .repeat.before = function(x)
  {
    ind = which(!is.na(x))
    if(is.na(x[1]))
      ind = c(1,ind)
    rep(x[ind], times = diff(c(ind, length(x) + 1)))
  }
  pat_tab$filter <- .repeat.before(pat_tab$filter)
  pat_tab$filter[is.na(pat_tab$filter)] <- filter_fin

  pat_tab$speed <- .repeat.before(pat_tab$speed)
  # if speed is not specified in the module, use the speed argument
  pat_tab$speed[is.na(pat_tab$speed)] <- speed

  pat_tab$tempo <- .repeat.before(pat_tab$tempo)
  # if tempo is not specified in the module, use the tempo argument
  pat_tab$tempo[is.na(pat_tab$tempo)] <- tempo

  r.sel <- which(is.na(pat_tab$duration))
  if (length(r.sel) > 0)
  {
    pat_tab$duration[r.sel] <- with(pat_tab[r.sel,], (delay + 1)*vblank_duration*speed*125/tempo)
    pat_tab$cum_duration[r.sel] <- unlist(lapply(as.list(1:(diff(range(r.sel)) + 1)), function(x){
      sum(pat_tab$duration[r.sel][1:x])
    }))
    pat_tab$cum_duration[r.sel] <- pat_tab$cum_duration[r.sel] + cum_duration
  }
  return (pat_tab)
}

.generate.channel.data <- function(mod, pt, track.nr, target.rate, video = c("PAL", "NTSC"), low.pass.filter = T, verbose = F)
{
  if (!verbose) cat(paste("Track ", track.nr, ":\n", sep =""))
  video <- match.arg(video)
  vblank_duration <- ifelse(video == "PAL", 1/50, 1/60)

  ramp.down <- function(x)
  {
    x <- 255 - round(510*(x%%64)/63)
    x[x > 0] <- 255 - x[x > 0]
    return(x)
  }
  # the ramp down for vibrato seems to differ from that of tremolo
  ramp.down.vib <- function(x)
  {
    x <- round(510*((x + 32)%%64)/63) - 255
    return(x)
  }

  square.wave <- function(x)
  {
    ifelse(x%%64 < 32, 255, -255)
  }

  result                 <- cbind(pt[,grepl(paste(".track", track.nr, sep =""),
                                            names(pt))], tempo = pt$tempo, filter = pt$filter)
  result                 <- result[as.vector(unlist(mapply(function(x, y) rep(y,x),
                                                           pt$speed*(1 + pt$delay), 1:nrow(pt)))),]
  names(result)          <- gsub("([.]track)[0-9]+$", "", names(result))
  result$tick            <- unlist(lapply(as.list(pt$speed*(1 + pt$delay)), function(x) 1:x))
  result$target.nsamples <- target.rate*125*vblank_duration/result$tempo
  decim <- 1/(result$target.nsamples %% 1)
  decim <- ifelse((1:nrow(result))%%(decim) < 1 & !is.infinite(decim), 1, 0)
  result$target.nsamples <- floor(result$target.nsamples) + decim
  rm(decim)

  ##############################################
  ## E9x retrigger each xth tick, even if no new note is played
  ## start
  ##############################################
  retrig  <- with(result, effect == as.raw(0x0e) & hiNybble(effect.mag) == 0x09 &
                    (tick - 1) %% loNybble(effect.mag) == 0)

  retrig  <- retrig | with(result, (note != "---" & !(effect %in% as.raw(c(0x03, 0x05)))) & tick == 1)
  ##############################################
  ## E9x retrigger each xth tick, even if no new note is played
  ## end
  ##############################################

  ##############################################
  ## EDx
  ## start
  ##############################################
  ## XXX it seems that this command should be ignored when x > speed. Normaly,
  ## this isn't a problem, but it is when used in combination with EEx. This
  ## still needs to be implemented.
  effect_sel <- with(result, effect == as.raw(0x0e) & hiNybble(effect.mag) == 0x0d)
  retrig[effect_sel] <- with(result[effect_sel,], tick == (1 + loNybble(effect.mag)))
  ##############################################
  ## EDx
  ## end
  ##############################################

  retrig  <- diff(c(1, which(retrig)))
  retrig  <- retrig[retrig != 0]
  retrig  <- unlist(lapply(as.list(retrig), function(x) 1:x))
  len.rem <- nrow(result) - length(retrig)
  if(len.rem > 0) retrig <- c(retrig, 1:len.rem)
  result$retrigger.sample <- retrig
  rm(retrig, len.rem)

  if (!verbose) cat("    Processing volume effects...")
  ## Set volume and volume effects:

  ## retrigger volume when a new note is played (that is not sample.nr == 0),
  ## except when porta to note is used,
  ## and retrigger when a sample number is given without a note.
  retrig_vol <- with(result, sample.nr > 0 & (retrigger.sample == 1 | tick == 1))

  effect_sel <- with(result, note != "---" | (effect == as.raw(0x0e) &
                                                effect.mag %in% as.raw(0x90:0x9F)))
  sample.nr  <- result$sample.nr
  sample.nr[-1][sample.nr[-1] == 0] <- NA
  sample.nr  <- .fill.parameter(sample.nr, is.na(sample.nr), function(x, y, z){rep(z, times = y - x + 1)})
  result$sample.nr.original <- result$sample.nr
  result$sample.nr[which(effect_sel & result$sample.nr == 0)] <-
    sample.nr[which(effect_sel & result$sample.nr == 0)]
  rm(sample.nr)

  result$volume[retrig_vol] <-
    unlist(lapply(as.list(result$sample.nr[retrig_vol]),
                  function(x) ifelse(x == 0, 0, volume(mod@samples[[x]]))))
  rm(retrig_vol, effect_sel)

  ##############################################
  ## ECx - cut volume after x ticks to zero
  ## start
  ##############################################
  effect_sel <- result$effect == as.raw(0x0E) & hiNybble(result$effect.mag) == 0x0C
  result$volume[effect_sel] <-
    unlist(lapply(as.list(result$sample.nr[effect_sel]),
                  function(x) ifelse(x == 0, 0, volume(mod@samples[[x]]))))
  effect_sel <- result$effect == as.raw(0x0E) & hiNybble(result$effect.mag) == 0x0C & result$tick > loNybble(result$effect.mag)
  result$volume[effect_sel] <- 0
  ##############################################
  ## ECx - cut volume after x ticks to zero
  ## end
  ##############################################

  ##############################################
  ## Cxy - set volume (start)
  ##############################################
  result$volume[result$effect == as.raw(0x0c)] <-
    as.integer(result$effect.mag[result$effect == as.raw(0x0c)])

  result$volume[1][is.na(result$volume[1])] <- 0
  ##############################################
  ## Cxy - set volume (end)
  ##############################################

  while (any(is.na(result$volume)))
  {
    ## For effect commands that do not affect volume, fill with last known volume value
    effect_sel <- with(result,
                       !(effect %in% as.raw(c(0x0a, 0x0c, 0x05, 0x06))) &
                         !(effect == as.raw(0x0e) & effect.mag >= as.raw(0xA0) & effect.mag <= as.raw(0xBF)))
    result$volume <-.fill.parameter(result$volume, effect_sel, function(x, y, z){rep(z, times = y - x + 1)})

    ##############################################
    ## Volume slide effects (start)
    ## Axy
    ## 5xy
    ## 6xy
    ##############################################
    effect_sel <- with(result, effect %in% as.raw(c(0x0A, 0x05, 0x06)))
    if (any(effect_sel))
    {
      effect.mag <- -loNybble(result$effect.mag)
      effect.mag[result$effect.mag > 0x0F] <- hiNybble(result$effect.mag[result$effect.mag > 0x0F])
      effect.mag[result$tick == 1] <- 0
      result$volume <- .fill.parameter(result$volume, effect_sel,
                                       function(x, y, z, mag){
                                         cums <- cumsum(mag[x:y])
                                         cums[z + cums < 0]    <- -z
                                         cums[z + cums > 0x40] <- 0x40 - z
                                         return(cums + z)
                                       }, effect.mag)
      rm(effect.mag)
    }
    ##############################################
    ## Volume slide effects (end)
    ## Axy
    ## 5xy
    ## 6xy
    ##############################################

    ##############################################
    ## Volume fine slide effects (start)
    ## EAx
    ## EBx
    ##############################################
    effect_sel <- with(result, effect == as.raw(0x0E) & hiNybble(effect.mag) %in% c(0x0A, 0x0B))
    if (any(effect_sel))
    {
      effect.mag <- with(result, ifelse(hiNybble(effect.mag) == 0x0A, 1, -1)*loNybble(result$effect.mag))
      effect.mag[result$tick != 1] <- 0
      result$volume <- .fill.parameter(result$volume, effect_sel,
                                       function(x, y, z, mag){
                                         cums <- cumsum(mag[x:y])
                                         cums[z + cums < 0]    <- -z
                                         cums[z + cums > 0x40] <- 0x40 - z
                                         return(cums + z)
                                       }, effect.mag)
      rm(effect.mag)
    }
    ##############################################
    ## Volume fine slide effects (end)
    ## EAx
    ## EBx
    ##############################################

  }
  ##############################################
  ## E7x - Set tremolo waveform (start)
  ##############################################
  effect_sel <- with(result, effect == as.raw(0x0e) & effect.mag %in% as.raw(0x70:0x7F))
  trem.wf    <- rep(NA, nrow(result))
  trem.wf[effect_sel] <- loNybble(result$effect.mag[effect_sel])
  trem.wf[1][is.na(trem.wf[1])] <- 0
  trem.wf <- .fill.parameter(trem.wf, !effect_sel, function(x, y, z){rep(z, times = y - x + 1)})
  ##############################################
  ## E7x - Set tremolo waveform (end)
  ##############################################

  ##############################################
  ## 7xy - Tremolo (start)
  ##############################################
  effect_sel    <- result$effect == as.raw(0x07)
  if (any(effect_sel))
  {
    table.pos.lo  <- loNybble(result$effect.mag[effect_sel]) == 0
    effect.mag.lo <- integer(nrow(result))
    effect.mag.lo[effect_sel] <- loNybble(result$effect.mag[effect_sel])
    effect.mag.lo[effect_sel][table.pos.lo] <- NA
    effect.mag.lo <-.fill.parameter(effect.mag.lo, is.na(effect.mag.lo), function(x, y, z){rep(z, times = y - x + 1)})
    table.pos.hi  <- hiNybble(result$effect.mag[effect_sel]) == 0
    effect.mag.hi <- integer(nrow(result))
    effect.mag.hi[effect_sel] <- hiNybble(result$effect.mag[effect_sel])
    effect.mag.hi[effect_sel][table.pos.hi] <- NA
    effect.mag.hi <-.fill.parameter(effect.mag.hi, is.na(effect.mag.hi), function(x, y, z){rep(z, times = y - x + 1)})
    effect.mag <- as.raw(effect.mag.hi*0x10 + effect.mag.lo)
    tremolo <- integer(nrow(result))
    tremolo[effect_sel][diff(c(-999, which(effect_sel))) == 1] <- NA
    tremolo[result$retrigger.sample == 1] <- 0
    tremolo.pos <- tremolo
    tremolo.pos <- .fill.parameter(tremolo.pos, is.na(tremolo.pos) & trem.wf < 4,
                                   function(x, y, z, tck){
                                     trm <- integer(y-x)
                                     sel <- tck[x:y] != 1
                                     trm[sel] <- seq(0, sum(sel) - 1)
                                     return(trm)
                                   }, result$tick)
    if (any(is.na(tremolo.pos)))
    {
      tremolo.pos.start <- tremolo.pos + 1
      tremolo.pos.start[tremolo.pos.start == 1] <- NA
      tremolo.pos.start[1][is.na(tremolo.pos.start[1])] <- 0
      tremolo.pos.start <- .fill.parameter(tremolo.pos.start, is.na(tremolo.pos.start), function(x, y, z){rep(z, times = y - x + 1)})
      tremolo.pos[is.na(tremolo.pos)] <- .fill.parameter(tremolo.pos[is.na(tremolo.pos)], rep(T, sum(is.na(tremolo.pos))),
                                                         function(x, y, z, tck){
                                                           trm <- integer(y-x)
                                                           sel <- tck[is.na(tremolo.pos)][x:y] != 1
                                                           start.val <- tremolo.pos.start[is.na(tremolo.pos)][x]
                                                           if (length(start.val) == 0 || is.na(start.val)) start.val <- 0
                                                           trm[sel] <- seq(start.val, start.val + sum(sel) - 1)
                                                           return(trm)
                                                         }, result$tick)
      rm(tremolo.pos.start)
    }
    tremolo <- .fill.parameter(tremolo, is.na(tremolo),
                               function(x, y, z, mag, tck, t.waveform, trem.pos){
                                 trm <- integer(y-x)
                                 spd <- hiNybble(as.raw(mag[x:y]))
                                 dep <- loNybble(as.raw(mag[x:y]))
                                 sel <- tck[x:y] != 1
                                 ## set waveform based on E7x command:
                                 twf <- proTrackerVibrato
                                 if ((t.waveform[x]%%4) == 1) twf <- ramp.down
                                 if ((t.waveform[x]%%4) %in% 2:3) twf <- square.wave
                                 trm[sel] <-
                                   as.integer(dep[sel]*twf(trem.pos[x:y][sel]*spd[sel])/64)
                                 return(trm)
                               }, effect.mag, result$tick, trem.wf, tremolo.pos)
    result$volume[effect_sel] <- result$volume[effect_sel] + tremolo[effect_sel]
    rm(tremolo, tremolo.pos, trem.wf)
  }
  result$volume[result$volume < 0x00] <- 0x00
  result$volume[result$volume > 0x40] <- 0x40
  ##############################################
  ## 7xy - Tremolo (end)
  ##############################################

  if (!verbose) cat("\tdone\n")
  if (!verbose) cat("    Processing period effects...")

  result$period <- NA
  effect_sel <- with(result, retrigger.sample == 1 &
                       !(note == "---" &
                           effect == as.raw(0x0e) &
                           effect.mag %in% as.raw(0x90:0x9F)))
  result$period[effect_sel] <-
    noteToPeriod(result$note[effect_sel],
                 unlist(lapply(as.list(result$sample.nr[effect_sel]),
                               function(x) ifelse(x == 0, 0, fineTune(mod@samples[[x]])))))
  result$period[effect_sel][result$period[effect_sel] < noteToPeriod("B-3")] <- noteToPeriod("B-3")

  ##############################################
  ## E5x set finetune (start)
  ##############################################
  ## only if sample is retriggered? (XXX not true, E5x, has some quirky behaviour, needs to be fixed):
  effect_sel <- with(result, retrigger.sample == 1 & effect == as.raw(0x0e) &
                       hiNybble(effect.mag) == 0x05)
  result$period[effect_sel] <-
    noteToPeriod(result$note[effect_sel],
                 nybbleToSignedInt(result$effect.mag[effect_sel], "low"))
  result$period[effect_sel][result$period[effect_sel] < noteToPeriod("B-3")] <- noteToPeriod("B-3")
  rm(effect_sel)
  ##############################################
  ## E5x set finetune (end)
  ##############################################

  ##############################################
  ## 3xy Porta to note initializing target
  ## 5xy
  ##############################################
  effect_sel    <- result$effect %in% as.raw(c(0x03, 0x05))
  snr           <- result$sample.nr
  snr[snr == 0] <- NA
  snr[1][is.na(snr[1])] <- 0
  snr           <-.fill.parameter(snr, is.na(snr), function(x, y, z){rep(z, times = y - x + 1)})
  target        <- paste(as.character(result$note), sprintf("%02i", snr))

  target2 <- target
  target2[effect_sel][substr(target2[effect_sel], 1, 3) == "---"] <- NA
  target2[1][is.na(target2[1])] <- paste("---", sprintf("%02i", snr[1]))

  target[substr(target, 1, 3) == "---" | !effect_sel] <- NA

  target2 <-.fill.parameter(target2, is.na(target2), function(x, y, z){rep(z, times = y - x + 1)})
  rm(snr)

  potential.target <- rep(T, nrow(result))
  potential.target[result$note == "---" | result$tick != 1] <- NA
  ##############################################
  ## 3xy Porta to note end initializing targets
  ## 5xy
  ##############################################

  while (any(is.na(result$period)))
  {
    ## For effect commands that do not affect period, fill with last known volume value
    effect_sel <- with(result,
                       !(effect %in% as.raw(c(0x01, 0x02, 0x03, 0x05))) &
                         !(effect == as.raw(0x0e) & hiNybble(effect.mag) %in% c(0x01, 0x02)))
    result$period <-.fill.parameter(result$period, effect_sel, function(x, y, z){rep(z, times = y - x + 1)})

    ##############################################
    ## 3xy porta to note
    ## 5xy porta to note + volume slide
    ## start
    ##############################################
    effect_sel  <- result$effect %in% as.raw(c(0x03, 0x05))
    effect.mag <- as.integer(result$effect.mag)
    effect.mag[effect_sel & effect.mag == 0] <- NA

    ## 5xy Does not affect the magnitude of the porta:
    effect.mag[result$effect == as.raw(0x05)] <- NA

    effect.mag[effect_sel] <-.fill.parameter(effect.mag[effect_sel], is.na(effect.mag[effect_sel]), function(x, y, z){rep(z, times = y - x + 1)})
    effect.mag[!effect_sel | result$tick == 1] <- 0

    target[1][is.na(target[1])] <- paste("---", sprintf("%02i", result$sample.nr[1]))

    ## tnr2 = a modified version of 'target_not_reached', taking
    ## into account whether the target has been reset
    tnr2 <- rep(T, nrow(result))
    ##  once a target is reached porta to note should stop sliding
    ## if target is not reached but a sample is retriggered, the target should not be reset!
    if (any(effect_sel))
    {
      target_not_reached <- rep(F, nrow(result))
      target.finetune <- as.numeric(substr(target[!is.na(target)], 5, 6))
      target.finetune <- unlist(lapply(as.list(as.numeric(substr(target[!is.na(target)], 5, 6))), function(x) ifelse(x == 0, 0, fineTune(mod@samples[[x]]))))
      temp <- noteToPeriod(substr(target[!is.na(target)], 1, 3), finetune = target.finetune)
      temp[temp < noteToPeriod("B-3")] <- noteToPeriod("B-3")
      target_not_reached[!is.na(target)]  <- result$period[!is.na(target)] != temp
      rm(temp)
      target_not_reached[is.na(result$period)] <- NA
      rm(target.finetune)

      ## If it's not known whether the target is reached, don't do anything yet!:
      unknown_target <- which(result$retrigger.sample == 1 & c(F, is.na(target_not_reached[-nrow(result)])) & !effect_sel)
      target[unknown_target] <- "???"

      target_not_reached[is.na(target_not_reached)] <- T
      ## sample is retriggered and target is reached and effect is not 3xy or 5xy:
      ## reset the target:
      reset_target <- which(result$retrigger.sample == 1 & c(F, !target_not_reached[-nrow(result)]) & !effect_sel)
      target[reset_target] <- target2[reset_target]

      ## if no new target is specified and the target is reached don't slide
      repeat{
        tnr2 <- rep(NA, nrow(result))
        tnr2[!(is.na(potential.target))] <- T
        tnr2 <- .fill.parameter(tnr2, is.na(potential.target), function (x,y,z, tnr, fx){
          res <- rep(T, y - x + 1)
          res[fx[x:y]] <- as.logical(cumprod(tnr[x:y][fx[x:y]]))
          res
        }, target_not_reached, effect_sel)

        pt <- potential.target
        pt[is.na(pt)] <- F
        potential.target[which(pt & !effect_sel & !c(T,tnr2[-length(tnr2)]))] <- NA
        if (length(which(pt)) == length(which(potential.target))) break
      }
    }
    target  <-.fill.parameter(target,  is.na(target),  function(x, y, z){rep(z, times = y - x + 1)})
    target[target == "???" | !is.na(result$period)] <- NA
    result$period <- .fill.parameter(result$period, effect_sel & !is.na(target),
                                     function(x, y, z, mag, target, tnr){
                                       samp <- as.numeric(substr(target, 5, 6))
                                       ft <- unlist(lapply(as.list(samp[x:y]), function(x) ifelse(x == 0, 0, fineTune(mod@samples[[x]]))))
                                       target <- substr(target, 1, 3)
                                       tar.per <- noteToPeriod(target[x:y], ft)
                                       tar.per[tar.per < noteToPeriod("B-3")] <- noteToPeriod("B-3")
                                       sgn <- sign(tar.per - z)
                                       cums <- sgn*cumsum(mag[x:y])
                                       cums[which(sgn*(tar.per - (z + cums)) < 0)] <- (tar.per - z)[which(sgn*(tar.per - (z + cums)) < 0)]
                                       # XXX do I need to check whether the period is out of range? Probably not possible since it slides to valid target
                                       cums[!tnr[x:y]] <- 0
                                       return(cums + z)
                                     }, effect.mag, target, tnr2)
    rm(effect.mag)
    ##############################################
    ## 3xy porta to note
    ## 5xy porta to note + volume slide
    ## end
    ##############################################

    ##############################################
    ## Porta up
    ## 1xy
    ## start
    ##############################################
    effect_sel <- (result$effect == as.raw(0x01))
    if (any(effect_sel))
    {
      effect.mag <- -1*as.integer(result$effect.mag)
      effect.mag[result$tick == 1] <- 0
      result$period <- .fill.parameter(result$period, effect_sel,
                                       function(x, y, z, mag){
                                         if (is.na(z)) return(rep(NA, y - x + 1))
                                         while (T)
                                         {
                                           #browser()
                                           cums <- cumsum(mag[x:y])
                                           res <- bitwAnd(0xFFF, cums + z)

                                           sel <- which((res < noteToPeriod("B-3")) & mag[x:y] < 1)
                                           if (length(sel) > 0 && !any(is.na(sel)))
                                           {
                                             end_sel <- which(diff(c(sel[1] - 1, sel)) > 1)
                                             if (length(end_sel) > 0) sel <- sel[-(end_sel[[1]]:length(sel))]
                                             if (length(sel) == 1 && sel < length(x:y) && mag[x:y][sel + 1] < 0 && res[sel + 1] > res[sel])
                                               mag[x:y][sel + 1] <- mag[x:y][sel + 1] + diff(c(0, res[sel] - noteToPeriod("B-3")))
                                             mag[x:y][sel] <-
                                               mag[x:y][sel] - diff(c(0, res[sel] - noteToPeriod("B-3")))
                                           }
                                           cums <- cumsum(mag[x:y])
                                           res <- bitwAnd(0xFFF, cums + z)
                                           if (!any(res < noteToPeriod("B-3"))) break
                                         }
                                         return(res)
                                       }, effect.mag)
      rm(effect.mag)
    }
    ##############################################
    ## Porta up
    ## 1xy
    ## end
    ##############################################

    ##############################################
    ## Porta down
    ## 2xy
    ## start
    ##############################################
    effect_sel <- (result$effect == as.raw(0x02))
    if (any(effect_sel))
    {
      effect.mag <- as.integer(result$effect.mag)
      effect.mag[result$tick == 1] <- 0
      result$period <- .fill.parameter(result$period, effect_sel,
                                       function(x, y, z, mag){
                                         if (is.na(z)) return(rep(NA, y - x + 1))
                                         if (z > 856) z <- 856
                                         while (T)
                                         {
                                           cums <- cumsum(mag[x:y])
                                           res <- cums + z
                                           sel <- which((res > 856) & mag[x:y] >= 0)
                                           if (length(sel) > 0 && !any(is.na(sel)))
                                           {
                                             end_sel <- which(diff(c(sel[1] - 1, sel)) > 1)
                                             if (length(end_sel) > 0) sel <- sel[-(end_sel[[1]]:length(sel))]
                                             mag[x:y][sel] <-
                                               mag[x:y][sel] + diff(c(0, 856 - res[sel]))
                                             cums <- cumsum(mag[x:y])
                                             res <- cums + z
                                           }
                                           if (!any(res > 856 & mag[x:y] > 0)) break
                                         }
                                         return(res)
                                       }, effect.mag)
      rm(effect.mag)
    }
    ##############################################
    ## Porta down
    ## 2xy
    ## end
    ##############################################

    ##############################################
    ## Fine porta up/down
    ## E1x
    ## E2x
    ## start
    ##############################################
    effect_sel <- with(result, effect == as.raw(0x0E) & hiNybble(effect.mag) %in% c(0x01, 0x02))
    if (any(effect_sel))
    {
      effect.mag <- with(result, ifelse(hiNybble(effect.mag) == 0x01, -1, 1)*loNybble(result$effect.mag))
      effect.mag[result$tick != 1] <- 0
      result$period <- .fill.parameter(result$period, effect_sel,
                                       function(x, y, z, mag){
                                         if (is.na(z)) return(rep(NA, y - x + 1))
                                         while (T)
                                         {
                                           cums <- cumsum(mag[x:y])
                                           res <- cums + z
                                           sel <- which((res < noteToPeriod("B-3")) & mag[x:y] < 1)
                                           if (length(sel) > 0 && !any(is.na(sel)))
                                           {
                                             end_sel <- which(diff(c(sel[1] - 1, sel)) > 1)
                                             if (length(end_sel) > 0) sel <- sel[-(end_sel[[1]]:length(sel))]
                                             mag[x:y][sel] <-
                                               mag[x:y][sel] - diff(c(0, res[sel] - noteToPeriod("B-3")))
                                           }
                                           cums <- cumsum(mag[x:y])
                                           res <- cums + z
                                           sel <- which((res > 856) & mag[x:y] > -1)
                                           if (length(sel) > 0 && !any(is.na(sel)))
                                           {
                                             end_sel <- which(diff(c(sel[1] - 1, sel)) > 1)
                                             if (length(end_sel) > 0) sel <- sel[-(end_sel[[1]]:length(sel))]
                                             mag[x:y][sel] <-
                                               mag[x:y][sel] + diff(c(0, 856 - res[sel]))
                                           }
                                           if (!(any(res < noteToPeriod("B-3")) || any(res > 0xffff))) break
                                         }
                                         return(res)
                                       }, effect.mag)
      rm(effect.mag)
    }
    ##############################################
    ## Fine porta up/down
    ## E1x
    ## E2x
    ## end
    ##############################################

  }
  rm(target, target2)

  ##############################################
  ## Fix period values above range - start
  ##############################################

  #browser()
  out.of.range <- rep(NA, nrow(result))
  out.of.range[1] <- F
  out.of.range[result$retrigger.sample == 1] <- F
  out.of.range[result$period > 856] <- T
  out.of.range <- .fill.parameter(out.of.range, is.na(out.of.range), function(x, y, z){rep(z, times = y - x + 1)})
  ## XXX the value of 0xffff does not seem to be correct
  ## the value seem to differ for different test cases, but is always high.
  ## Need to figure out how to get the right value
  result$period[out.of.range &
                  !(result$effect == as.raw(0x01) & result$tick > 1) &
                  !(result$effect == as.raw(0x02) & result$tick > 1)] <- 0xffff
  rm(out.of.range)

  ##############################################
  ## Fix period values above range - end
  ##############################################

  ##############################################
  ## E4x - Set vibrato waveform (start)
  ##############################################
  effect_sel <- with(result, effect == as.raw(0x0e) & effect.mag %in% as.raw(0x40:0x4F))
  vibr.wf    <- rep(NA, nrow(result))
  vibr.wf[effect_sel] <- loNybble(result$effect.mag[effect_sel])
  vibr.wf[1][is.na(vibr.wf[1])] <- 0
  vibr.wf <- .fill.parameter(vibr.wf, !effect_sel, function(x, y, z){rep(z, times = y - x + 1)})
  ##############################################
  ## E4x - Set vibrato waveform (end)
  ##############################################

  ##############################################
  ## 4xy - Vibtrato
  ## 6xy - Vibrato + volume slide (the latter part is already dealt with above)
  ## start
  ##############################################

  effect_sel    <- result$effect %in% as.raw(c(0x04, 0x06))
  if (any(effect_sel))
  {
    table.pos.lo  <- loNybble(result$effect.mag[effect_sel]) == 0
    effect.mag.lo <- integer(nrow(result))
    effect.mag.lo[effect_sel] <- loNybble(result$effect.mag[effect_sel])
    effect.mag.lo[effect_sel][table.pos.lo] <- NA
    effect.mag.lo <-.fill.parameter(effect.mag.lo, is.na(effect.mag.lo), function(x, y, z){rep(z, times = y - x + 1)})
    table.pos.hi  <- hiNybble(result$effect.mag[effect_sel]) == 0
    effect.mag.hi <- integer(nrow(result))
    effect.mag.hi[effect_sel] <- hiNybble(result$effect.mag[effect_sel])
    effect.mag.hi[effect_sel][table.pos.hi] <- NA
    effect.mag.hi <-.fill.parameter(effect.mag.hi, is.na(effect.mag.hi), function(x, y, z){rep(z, times = y - x + 1)})
    effect.mag <- as.raw(effect.mag.hi*0x10 + effect.mag.lo)
    vibrato <- integer(nrow(result))
    vibrato[effect_sel][diff(c(-999, which(effect_sel))) == 1] <- NA
    vibrato[result$retrigger.sample == 1] <- 0
    vibrato.pos <- vibrato
    vibrato.pos <- .fill.parameter(vibrato.pos, is.na(vibrato.pos) & vibr.wf < 4,
                                   function(x, y, z, tck){
                                     vib <- integer(y-x)
                                     sel <- tck[x:y] != 1
                                     vib[sel] <- seq(0, sum(sel) - 1)
                                     return(vib)
                                   }, result$tick)
    if (any(is.na(vibrato.pos)))
    {
      vibrato.pos.start <- vibrato.pos + 1
      vibrato.pos.start[vibrato.pos.start == 1] <- NA
      vibrato.pos.start[1][is.na(vibrato.pos.start[1])] <- 0
      vibrato.pos.start <- .fill.parameter(vibrato.pos.start, is.na(vibrato.pos.start), function(x, y, z){rep(z, times = y - x + 1)})
      vibrato.pos[is.na(vibrato.pos)] <- .fill.parameter(vibrato.pos[is.na(vibrato.pos)], rep(T, sum(is.na(vibrato.pos))),
                                                         function(x, y, z, tck){
                                                           vib <- integer(y-x)
                                                           sel <- tck[is.na(vibrato.pos)][x:y] != 1
                                                           start.val <- vibrato.pos.start[is.na(vibrato.pos)][x]
                                                           if (length(start.val) == 0 || is.na(start.val)) start.val <- 0
                                                           vib[sel] <- seq(start.val, start.val + sum(sel) - 1)
                                                           return(vib)
                                                         }, result$tick)
      rm(vibrato.pos.start)
    }
    vibrato <- .fill.parameter(vibrato, effect_sel,
                               function(x, y, z, mag, tck, v.waveform, vib.pos){
                                 vib <- integer(y-x)
                                 spd <- hiNybble(as.raw(mag[x:y]))
                                 dep <- loNybble(as.raw(mag[x:y]))
                                 sel <- tck[x:y] != 1
                                 ## set waveform based on E4x command:
                                 vwf <- proTrackerVibrato
                                 if ((v.waveform[x]%%4) == 1) vwf <- ramp.down.vib
                                 if ((v.waveform[x]%%4) %in% 2:3) vwf <- square.wave
                                 vib[sel] <-
                                   as.integer(dep[sel]*vwf(vib.pos[x:y][sel]*spd[sel])/128)
                                 return(vib)
                               }, effect.mag, result$tick, vibr.wf, vibrato.pos)
    result$period[effect_sel] <- result$period[effect_sel] + vibrato[effect_sel]
    rm(vibrato)
  }

  ##############################################
  ## 4xy - Vibtrato
  ## 6xy - Vibrato + volume slide (the latter part is already dealt with above)
  ## end
  ##############################################

  ##############################################
  ## 0xy Arpeggio (start)
  ##############################################

  ## XXX check if period values are out of range!
  effect_sel  <- result$effect == raw(1) & result$effect.mag > raw(1)
  if (any(effect_sel))
  {
    effect.mag  <- integer(nrow(result))
    period_tab2 <- expand.grid(finetune=-8:7, note=names(ProTrackR::period_table[,-1:-2]), octave= 1:3)
    period_tab2$period <- apply(period_tab2, 1, function(x){
      result <- noteToPeriod(paste(x["note"], x["octave"], sep =""),
                             finetune = x["finetune"])
      result[result < noteToPeriod("B-3")] <- noteToPeriod("B-3")
      result
    })
    target1     <- unlist(lapply(as.list(result$period), function(x){
      res <- which(abs(period_tab2$period - x) == min(abs(period_tab2$period - x)) &
                     abs(period_tab2$finetune) == min(abs(period_tab2$finetune)))
      if (length(res) == 0) return (NA) else return(res)
    }))
    target2     <- period_tab2$period[target1 + 16*loNybble(result$effect.mag)] -
      period_tab2$period[target1 + 0]
    target1     <- period_tab2$period[target1 + 16*hiNybble(result$effect.mag)] -
      period_tab2$period[target1 + 0]
    target1[is.na(target1)] <- 0
    target2[is.na(target2)] <- 0
    effect.mag[(result$tick - 1) %% 3 == 1] <- target1[(result$tick - 1) %% 3 == 1]
    effect.mag[(result$tick - 1) %% 3 == 2] <- target2[(result$tick - 1) %% 3 == 2]
    result$period[effect_sel] <- result$period[effect_sel] + effect.mag[effect_sel]
    rm(period_tab2, target1, target2)
  }
  ##############################################
  ## 0xy Arpeggio (end)
  ##############################################

  ##############################################
  ## Fix period values below range - start
  ##############################################

  result$period[result$period < noteToPeriod("B-3")] <- noteToPeriod("B-3")

  ##############################################
  ## Fix period values below range - end
  ##############################################

  if (!verbose) cat("\tdone\n")
  if (!verbose) cat("    Resampling...")

  # XXX when sample number switches from one to another while it is not retriggered
  # The sample will stop playing at the end of its loop (or the end of the sample
  # if there's no loop) and will continue playing at the start of the loop of the
  # next sample.
  # If no new sample is played set sample number to last played sample

  result$sample.rate <- periodToSampleRate(result$period, video = video)

  result$sample.pos  <- NA
  result$sample.pos[result$retrigger.sample == 1] <- 1

  ##############################################
  ## 9xy sample offset (start)
  ##############################################
  ## protracker buggy with 9xy command:
  ## ftp://ftp.modland.com/pub/documents/format_documentation/Tracker%20differences%20for%20Coders.txt
  ## also see test case ptoffset.mod...

  effect_sel <- result$effect == as.raw(0x09)
  if (any(effect_sel))
  {
    base.magnitude <- rep(NA, nrow(result))
    base.magnitude[result$sample.nr.original > 0 & effect_sel] <- as.integer(result$effect.mag[result$sample.nr.original > 0 & effect_sel])
    base.magnitude[1][is.na(base.magnitude[1])] <- 0

    snr <- (result$sample.nr.original != 0)
    snr[result$retrigger.sample != 1] <- NA
    snr[result$retrigger.sample == 1 & effect_sel] <- T
    snr[1][is.na(snr[1])] <- F
    snr <- .fill.parameter(snr, is.na(snr), function(x, y, z){rep(z, times = y - x + 1)})

    mod_bool <- rep(NA, nrow(result))
    mod_bool[with(result, note == "---" & effect_sel)] <- T
    mod_bool[with(result, (retrigger.sample == 1 | effect_sel) & sample.nr.original > 0)] <- F
    mod_bool <- .fill.parameter(mod_bool, is.na(mod_bool), function(x, y, z){rep(z, times = y - x + 1)})
    effect_sel2 <- with(result, !mod_bool & diff(c(-999, snr)) != 0 & sample.nr.original == 0 & note != "---" & !effect_sel & tick == 1) | with(result, effect_sel & sample.nr.original == 0 & tick == 1)
    rm(mod_bool)
    modify.magnitude <- rep(NA, nrow(result))
    effect.mag       <- rep(NA, nrow(result))
    effect.mag[effect_sel] <- as.integer(result$effect.mag[effect_sel])
    effect.mag <- .fill.parameter(effect.mag, is.na(effect.mag), function(x, y, z){rep(z, times = y - x + 1)})
    modify.magnitude[effect_sel2] <- as.integer(effect.mag[effect_sel2])

    base.magnitude.fill <- base.magnitude
    base.magnitude.fill <- .fill.parameter(base.magnitude, is.na(base.magnitude), function(x, y, z){rep(z, times = y - x + 1)})

    modify.magnitude[which(modify.magnitude == 0)] <- base.magnitude.fill[which(modify.magnitude == 0)]
    modify.magnitude[is.na(modify.magnitude)] <- 0
    #now fill the base magnitudes with the cumulative sum of modify magnitudes:

    effect.fill <- rep(NA, nrow(result))
    effect.fill[with(result, note != "---" & sample.nr.original != 0)] <- effect_sel[with(result, note != "---" & sample.nr.original != 0)]
    effect.fill[effect_sel] <- T
    effect.fill <- .fill.parameter(effect.fill, is.na(effect.fill), function(x, y, z){rep(z, times = y - x + 1)})

    base.magnitude[(result$sample.nr.original > 0 & !effect_sel) | (result$note != "---" & !effect.fill)] <- 0
    base.magnitude <- .fill.parameter(base.magnitude,
                                      is.na(base.magnitude),
                                      function(x, y, z, mod.mag){
                                        rep(z, times = y - x + 1) + cumsum(mod.mag[x:y])
                                      }, modify.magnitude)

    result$sample.pos[result$retrigger.sample == 1] <-
      base.magnitude[result$retrigger.sample == 1]*0x100 + 1
    rm(base.magnitude, base.magnitude.fill, effect_sel2, snr, effect.fill)
  } else
  {
    result$sample.pos[result$retrigger.sample == 1] <- 1
  }

  ##############################################
  ## 9xy sample offset (end)
  ##############################################

  result$sample.pos <- .fill.parameter(result$sample.pos, result$retrigger.sample != 1,
                                       function(x, y, z, rt, tp){
                                         z + cumsum(rt[(x:y) - 1]*125*vblank_duration/tp[(x:y) - 1])
                                       }, result$sample.rate, result$tempo)

  result$sample.nr[result$retrigger.sample != 1 & result$sample.nr == 0] <- NA
  result$sample.nr <- .fill.parameter(result$sample.nr, is.na(result$sample.nr),
                                      function(x, y, z){rep(z, times = y - x + 1)})

  ## Sample switching starts here. not implemented correctly yet
  result$sample.switch <- NA
  result$sample.switch[result$retrigger.sample == 1] <- F
  result$sample.switch[diff(c(0, result$sample.nr)) != 0 &
                         result$retrigger.sample != 1 &
                         result$sample.nr != 0] <- T
  result$sample.switch[1][is.na(result$sample.switch[1])] <- F
  result$sample.switch <- .fill.parameter(result$sample.switch, is.na(result$sample.switch),
                                          function(x, y, z){rep(z, times = y - x + 1)})

  ## take one line before to know from which sample to switch from
  samp.switch <- subset(result, result$sample.switch | c(result$sample.switch[-1], F))
  if(nrow(samp.switch) > 0)
  {
    samp.switch$index <- which(result$sample.switch | c(result$sample.switch[-1], F))
    group.start <- diff(c(-999, samp.switch$index)) > 1
    samp.switch$group <- NA
    samp.switch$group[group.start] <- 1:sum(group.start)
    samp.switch$group <- .fill.parameter(samp.switch$group, is.na(samp.switch$group),
                                         function(x, y, z){rep(z, times = y - x + 1)})
    for (i in unique(samp.switch$group))
    {
      samp.group <- subset(samp.switch, samp.switch$group == i)
      snr <- result$sample.nr[samp.group$index[[1]] - 1]
      if (snr > 0)
      {
        sample.from <- mod@samples[[snr]]
        ls <- loopStart(sample.from)
        st <- loopState(sample.from)
        if (st)
        {
          samp.group$sample.pos[samp.group$sample.pos > (ls + 1)] <-
            ((samp.group$sample.pos[samp.group$sample.pos > (ls + 1)] - (ls + 1)) %% loopLength(sample.from)) + ls + 1
          samp.group$sample.nr[-1][diff(samp.group$sample.pos) >= 0] <-
            samp.group$sample.nr[1]
          samp.group$sample.switch <- c(F, diff(samp.group$sample.pos) < 0)
        } else
        {
          sl <- sampleLength(sample.from)
          samp.group$sample.nr[samp.group$sample.pos < sl] <- samp.group$sample.nr[1]
          samp.group$sample.switch <- samp.group$sample.pos >= sl
        }
      }
      samp.switch[samp.switch$group == i,] <- samp.group
    }
    idx <- samp.switch$index
    samp.switch$index <- NULL
    samp.switch$group <- NULL
    result[idx,] <- samp.switch
  }
  # Don't use new sample without retrigger until the previous sample reaches its
  # (loop) end.

  ## current implementation of sample switching ends here...

  resamp_sel <- result$retrigger.sample == 1 | c(F, result$period[-nrow(result)] != result$period[-1])
  result$resamp.start <- NA
  result$resamp.start[resamp_sel] <- 1:sum(resamp_sel)
  result$resamp.start <- .fill.parameter(result$resamp.start, !resamp_sel, function(x, y, z){rep(z, times = y - x + 1)})
  rm(resamp_sel)

  duration <- 125*vblank_duration/result$tempo
  resamp.table <- stats::aggregate(cbind(target.nsamples, duration)~resamp.start, result, sum)
  resamp.table <- merge(resamp.table, aggregate(cbind(sample.nr, sample.rate, sample.pos)~resamp.start, result, function(x) x[[1]]))
  rm(duration)

  resamp.agg <- aggregate(row.nr~sample.nr+target.nsamples+sample.rate+sample.pos, cbind(resamp.table, row.nr = 1:nrow(resamp.table)), function(x) x)

  resamp.order <- lapply(as.list(1:nrow(resamp.table)), function(x) which(unlist(lapply(resamp.agg$row.nr, function(y) x %in% y))))

  channel <- apply(resamp.agg, 1, function(x){
    # Is this check on sample.pos really required?:
    if (is.infinite(x[["sample.pos"]]) || x[["sample.nr"]] == 0 || is.infinite(x[["sample.rate"]])) return(rep(128, x[["target.nsamples"]]))
    smp       <- mod@samples[[x[["sample.nr"]]]]
    startpos  <- x[["sample.pos"]]
    stoppos   <- startpos + x[["target.nsamples"]]*x[["sample.rate"]]/target.rate
    start2    <- (startpos %% 1)*target.rate/x[["sample.rate"]]
    out.range <- floor(1 + start2):(floor(start2) + x[["target.nsamples"]])
    wf        <- waveform(smp, floor(startpos), ceiling(stoppos))
    wf[is.na(wf)] <- 128
    wf        <- resample(wf, x[["sample.rate"]], target.rate, method = "constant")
    wf        <- wf[out.range]
    #    if (any(is.na(wf))) print(paste(x[["sample.nr"]], x[["sample.rate"]], startpos, stoppos, range(out.range)))
    return(list(wf))
  })
  channel <- unlist(lapply(resamp.order, function(x) unlist(channel[[x]])))

  rm(resamp.table)

  chan.vol <- as.vector(unlist(apply(result, 1, function(x){
    rep(as.integer(x["volume"]), x["target.nsamples"])
  })))

  channel <- as.integer(128 + ((channel - 128)*chan.vol/0x40))
  if (!verbose) cat("\t\t\tdone\n")

  ##############################################
  ## E0x turn filter on/off
  ## start
  ##############################################
  ## filter frequencies:
  ## http://forum.audacityteam.org/viewtopic.php?f=42&t=70941
  ## http://forum.arduino.cc/index.php?topic=21484.0
  ## https://bel.fi/alankila/modguide/interpolate.txt
  ## x == even -> filter on
  ## x == odd  -> filter off
  ## can only be applied when target.rate > 4900 (which is the band which is filtered)
  ## For later consideration:
  ## should the filter be applied here? Or before setting volume, of after mixing?
  channel <- channel - 128

  if (target.rate > 4900 && low.pass.filter)
  {
    if (!verbose) cat("    Applying low pass filter...")
    sel_filter <- as.vector(unlist(apply(result, 1, function(x){
      rep(x["filter"], as.integer(x["target.nsamples"]))
    })))
    sel_filter <- grepl("TRUE", sel_filter, fixed = T)
    ## "LED" filter
    butterworth <- signal::butter(2, 3275/target.rate, "low")
    channel[sel_filter] <- signal::filter(butterworth, channel[sel_filter])
    ## standard filter:
    butterworth <- signal::butter(1, 4900/target.rate, "low")
    channel[!sel_filter] <- signal::filter(butterworth, channel[!sel_filter])
    if (!verbose) cat("\t\tdone\n")
  }
  ##############################################
  ## E0x turn filter on/off
  ## end
  ##############################################
  return(channel + 128)
}

.fill.parameter <- function(par, selection, expr, ...)
{
  sel2        <- selection & is.na(par)
  start.range <- which(sel2)[which(diff(c(-999, which(sel2))) != 1)]
  stop.range  <- c(which(sel2)[which(diff(c(which(sel2))) != 1)], which(sel2)[sum(sel2)])
  start.sel   <- start.range - 1
  start.sel[start.sel == 0] <- 1
  start.par   <- par[start.sel]
  par[sel2]   <- as.vector(unlist(mapply(expr,
                                         start.range,
                                         stop.range,
                                         start.par,
                                         MoreArgs = list(...))))
  return(par)
}
