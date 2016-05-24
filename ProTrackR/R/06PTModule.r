validity.PTModule <- function(object)
{
  if (length(object@pattern.order.length) != 1)              return (F)
  if (!(as.integer(object@pattern.order.length) > 0 &&
        as.integer(object@pattern.order.length) < 129))      return (F)
  if (length(object@name)          != 20)                    return (F)
  if (length(object@pattern.order) != 128)                   return (F)
  if (length(object@tracker.byte)  != 1)                     return (F)
  if (object@tracker.byte != as.raw(0x7F))                   return (F)
  # We're only being compatible with ProTracker, which holds 31 samples
  if (length(object@samples)       != 31)                    return (F)
  if (!all(unlist(lapply(object@samples, class)) == "PTSample"))
    return (F)
  if (!all(unlist(lapply(object@samples, validObject))))     return (F)
  if (!all(unlist(lapply(object@patterns, class)) == "PTPattern"))
    return (F)
  if (!all(unlist(lapply(object@patterns, validObject))))    return (F)
  if (length(object@tracker.flag)  != 4)                     return (F)
  if (!all(object@tracker.flag     == charToRaw("M.K.")) &&
      !all(object@tracker.flag     == charToRaw("M!K!")))    return (F)
  tf <- all(object@tracker.flag    == charToRaw("M!K!"))
  if (length(object@patterns) > 64   && !tf)                 return (F)
  if (length(object@patterns) > 100  &&  tf)                 return (F)
  if (length(object@patterns) < 1)                           return (F)
  if (!all(unlist(lapply(object@patterns, class)) == "PTPattern"))
    return (F)
  if ((max(as.integer(object@pattern.order)) + 1) > length(object@patterns)) return (F)
  return (T)
}

#' The PTModule class
#'
#' The PTModule class provides a container to store and modify and use ProTracker
#' module files.
#'
#' MOD is a computer file format used primarily to represent music. A MOD file
#' contains a set of instruments in the form of samples, a number of patterns
#' indicating how and when the samples are to be played, and a list of what
#' patterns to play in what order. The simplified structure of a module class is
#' visualised in the scheme below. Details are given in the slot descriptions
#' below.
#'
#' \if{html}{\figure{protrackerscheme.png}{ProTracker conceptual scheme}}
#' \if{latex}{\figure{protrackerscheme.pdf}{options: width=6in}}
#'
#' This class is designed to hold all relevant information of a ProTracker
#' module (MOD) for which ProTracker 2.3a documentation was used. The ProTrackR
#' package may be compatible with earlier or later versions, but this was not
#' tested. Use \code{\link{read.module}} and \code{\link{write.module}} to import
#' and export objects of class \code{PTModule}.
#'
#' @slot name A \code{vector} of length 20 of class "\code{raw}", representing
#' the name of the \code{PTModule}. The name
#' of a module can be extracted or replaced with the \code{\link{name}} method.
#' @slot pattern.order A \code{vector} of length 128 of class "\code{raw}". The
#' \code{raw} values represent the indices of \code{PTPattern} tables and indicate
#' in which order these patterns need to be played. Note that the \code{raw} values
#' are conform the indices used in ProTracker, starting at zero. In R, indices of
#' objects start at one. Users need to compensate for this discrepancy theirselves.
#'
#' The pattern order table can be extracted or replaced with the
#' \code{\link{patternOrder}} method.
#' @slot pattern.order.length A single value of class "\code{raw}". Indicates
#' the length of the visible (and playable) part of the pattern order table.
#'
#' Use the \code{\link{patternOrderLength}} method to extract or replace the length
#' of a pattern order table of a module.
#' @slot tracker.byte A single "\code{raw}" value. Gives an indication of which
#' Tracker was used to produce a module file. In ProTracker modules, this byte
#' is set to 0x7f, which is also used in \code{PTModule} objects. This value
#' should not be changed.
#' @slot tracker.flag A \code{vector} of length 4 of class "\code{raw}", indicates
#' the version of a module, which basically reflects how many patterns the module
#' can hold. For details, and extracting and replacing this flag see the
#' \code{\link{trackerFlag}} method.
#' @slot samples List of length 31 of class "\code{\link{PTSample}}".
#' @slot patterns List of class "\code{\link{PTPattern}}" (the pattern tables).
#' The list should have at least 1 element, and can have a maximum of 64 or
#' 100 elements (depending on the state of the \code{\link{trackerFlag}}).
#'
#' @name PTModule-class
#' @rdname PTModule-class
#' @aliases PTModule
#' @references
#' \url{https://en.wikipedia.org/wiki/MOD_(file_format)}
#'
#' \url{http://wiki.multimedia.cx/index.php?title=Protracker_Module}
#'
#' \url{http://coppershade.org/articles/More!/Topics/Protracker_File_Format/}
#' @examples
#' ## create an empty PTModule class object:
#' mod.empty <- new("PTModule")
#'
#' ## get an example PTModule class object
#' ## provided with the ProTrackR package:
#' data("mod.intro")
#'
#' @family module.operations
#' @author Pepijn de Vries
#' @exportClass PTModule
setClass("PTModule",
         representation(name                 = "raw",
                        pattern.order        = "raw",
                        pattern.order.length = "raw",
                        tracker.byte         = "raw",
                        tracker.flag         = "raw",
                        samples              = "list",
                        patterns             = "list"),
         prototype(name = raw(20),
                   pattern.order        = rep(as.raw(0x00), 128),
                   pattern.order.length = as.raw(0x01),
                   tracker.byte         = as.raw(0x7F),
                   tracker.flag         = charToRaw("M.K."),
                   # We're only being compatible with ProTracker by supporting 31 samples
                   samples              = lapply(rep("PTSample", 31), new),
                   patterns             = list(new("PTPattern"))),
         validity = validity.PTModule)

#' Plot a PTModule object
#'
#' Plots the waveforms of the (non-empty) \code{\link{PTSample}}s in a
#' \code{\link{PTModule}} object.
#'
#' A plotting routine based on the \code{\link[lattice]{xyplot}} from the
#' lattice-package. Plots each (non-empty) waveform in a separate panel. Use arguments
#' of the \code{\link[lattice]{xyplot}} function to customise the plot.
#' @docType methods
#' @rdname plot
#' @name plot
#' @aliases plot,PTModule,missing-method
#' @param x A \code{\link{PTModule}} object for which the
#' waveforms of the \code{\link{PTSample}}s need to be plotted.
#' @param y \code{missing}. Argument from the generic plotting method, don't use.
#' @param plot.loop.positions A \code{logical} value indicating whether
#' loop positions need to be visualised. For looped samples, the starting
#' and ending positions are marked by a vertical green and red line, respectively.
#' @param ... Arguments that are passed on to \code{\link[lattice]{xyplot}}.
#' @return Returns an object of class \code{trellis}. See documentation of
#' \code{\link[lattice]{xyplot}} for more details.
#'
#' @examples
#' ## get the example PTModule provided with the ProTrackR package
#' data("mod.intro")
#'
#' ## The most basic way to plot the module samples:
#' plot(mod.intro)
#'
#' ## By using xyplot arguments, we can make it look nicer:
#' plot(mod.intro, type = "l", layout = c(1,4),
#'      scales = list(x = list(relation = "free")))
#' @author Pepijn de Vries
#' @export
setMethod("plot", c("PTModule", "missing"), function(x, y, plot.loop.positions = T, ...){
  amplitude   <- NULL
  samp_name   <- NULL
  `time (s)`  <- NULL
  plot.loop.positions <- as.logical(plot.loop.positions)
  for (i_sample in 1:31)
  {
    if (length(x@samples[[i_sample]]@left) > 0)
    {
      len         <- length(x@samples[[i_sample]]@left)
      amplitude   <- c(amplitude, x@samples[[i_sample]]@left)
      samp_name   <- c(samp_name,
                       rep(paste(sprintf("%02i", i_sample),
                                 rawToCharNull(x@samples[[i_sample]]@name)), len))
      `time (s)`  <- c(`time (s)`,
                       seq(0, (len - 1)/noteToSampleRate(),
                           length.out = len))
    }
  }

  if (plot.loop.positions)
  {
    has_data <- unlist(lapply(x@samples, function(x) (length(x@left) > 0)))
    start    <- unlist(lapply(x@samples, function(x) rawToUnsignedInt(x@wloopstart)))
    len      <- unlist(lapply(x@samples, function(x) rawToUnsignedInt(x@wlooplen)))
    panel.loop <- function(...)
    {
      if (!(start[has_data][lattice::packet.number()] == 0 &&
            len[has_data][lattice::packet.number()] == 1))
      {
        lattice::panel.abline(v =                2*start[has_data][lattice::packet.number()] /
                                noteToSampleRate(), col = "green")
        lattice::panel.abline(v = 2*((len[has_data] - 0.5) + start[has_data])[lattice::packet.number()] /
                                noteToSampleRate(), col = "red")
      }
    }
    lattice::xyplot(amplitude~`time (s)`|samp_name,
                    panel = function(...){lattice::panel.xyplot(...); panel.loop(...)}, ...)
  } else
  {
    lattice::xyplot(amplitude~`time (s)`|samp_name, ...)
  }
})

setMethod("show", "PTModule", function(object){
  print(object)
})

#' @rdname print
#' @aliases print,PTModule-method
#' @export
setMethod("print", "PTModule", function(x, ...){
  cat("\nPTModule Object:\n")
  cat(paste("\tModule name:" , rawToCharNull(x@name), "\n", sep = "\t\t\t"))
  cat(paste("\tNumber of samples:",
            sum(unlist(lapply(x@samples, function(x) length(x@left) > 0))), "\n", sep = "\t\t"))
  cat(paste("\tNumber of patterns:", 1+max(as.integer(x@pattern.order)), "\n", sep = "\t\t"))
  cat(paste("\tPattern order table length:", as.integer(x@pattern.order.length), "\n", sep = "\t"))
})

#' @rdname playSample
#' @aliases playSample,PTModule-method
#' @export
setMethod("playSample", "PTModule", function(x, silence, wait, note, loop, ...){
  silence  <- abs(as.numeric(silence[[1]]))
  wait     <- as.logical(wait[[1]])
  note     <- as.character(note[[1]])
  loop     <- abs(as.numeric(loop[[1]]))
  if (loop == 0) stop ("'loop' should be greater than 0.")
  wf_all   <- NULL
  for (i_sample in 1:31)
  {
    samp   <- PTSample(x, i_sample)
    if (length(samp) > 0)
    {
      ft     <- fineTune(samp)
      vl     <- volume(samp)
      wf     <- samp@left
      if ("finetune" %in% names(list(...)))
        sr <- noteToSampleRate(note, ...)
      else
        sr <- noteToSampleRate(note, ft, ...)
      if (loopState(samp))
      {
        n_samp <- round(loop*sr)
        if (loopStart(samp) + loopLength(samp) > n_samp) n_samp <- loopStart(samp) + loopLength(samp)
        wf <- loopSample(samp, n_samples = n_samp)
      }
      ## seewave not available for OS X replace resampling with custom
      ## resampling algorithm
      wf     <- resample(vl*(wf - 128)/(128*0x40), sr, 44100, method = "constant")
      wf_all <- c(wf_all, wf, rep(0, silence*44100))
    }
  }
  if (wait)
  {
    audio::wait(audio::play(wf_all, rate = 44100))
  } else
  {
    audio::play((wf_all - 128)/128, rate = 44100)
  }
  invisible()
})

setGeneric("read.module", function(file, ignore.validity = F) standardGeneric("read.module"))

#' Read a ProTracker module file
#'
#' Reads a ProTracker module file and coerces it to a \code{\link{PTModule}}
#' object.
#'
#' The routine to read ProTracker modules is based on the referenced version
#' of ProTracker 2.3A. This means that the routine may not be able to
#' read files produced with later ProTracker versions, or earlier versions with
#' back-compatibility issues. So far I've successfully tested this method
#' on all modules I've composed with ProTracker version 2.3A (which I believe
#' was one of the more popular versions of ProTracker back in the days).
#'
#' It should also be able to read most of the .mod files in
#' \href{http://modarchive.org/}{The Mod Archive}.
#'
#' @docType methods
#' @name  read.module
#' @rdname read.module
#' @aliases read.module,character-method
#' @param file either a filename or a file connection, that
#' allows reading binary data (see e.g., \code{\link[base]{file}} or \code{\link[base]{url}}).
#' @param ignore.validity A \code{logical} value indicating whether the
#' validity of the \code{PTModule} should be ignored. When set to
#' \code{FALSE} (default), the validity of the read object is checked; an
#' error is thrown when the object is not valid. When this argument is set to
#' \code{TRUE}, the validity of the object will not be checked and a potentially
#' invalid object is returned. As the validity check of \code{PTModule}
#' objects is very strict, it can be useful to ignore this check. This way
#' you can try to read a broken module file, try to fix it such that it becomes valid and
#' save (with \code{\link{write.module}}) it again.
#' @return Returns a \code{PTModule} object read from the provided ProTracker file
#'
#' @examples
#' \dontrun{
#'
#' ## first create an module file from example data:
#' data("mod.intro")
#' write.module(mod.intro, "intro.mod")
#'
#' ## read the module:
#' mod  <-  read.module("intro.mod")
#'
#' ## or create a connection yourself:
#' con  <- file("intro.mod", "rb")
#'
#' ## note that you can also read from URL connections!
#' mod2 <- read.module(con)
#'
#' ## don't forget to close the file:
#' close(con)
#' }
#' @references \url{http://wiki.multimedia.cx/index.php?title=Protracker_Module}
#'
#' \url{http://coppershade.org/articles/More!/Topics/Protracker_File_Format/}
#' @family io.operations
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("read.module", c("character"), function(file, ignore.validity){
  ignore.validity <- as.logical(ignore.validity)[[1]]
  con <- file(file, "rb")
  mod <- read.module(con, ignore.validity)
  close(con)
  return(mod)
})

#' @rdname read.module
#' @aliases read.module,ANY-method
#' @export
setMethod("read.module", c("ANY"), function(file, ignore.validity)
{
  ignore.validity <- as.logical(ignore.validity)[[1]]
  # function to read module from file
  # connection following the specs listed here:
  # http://wiki.multimedia.cx/index.php?title=Protracker_Module
  # check if we got a connection we can use (otherwise throw error):
  con <- file
  if (!("connection" %in% class(con))) stop ("argument con is not a file connection!")
  con_info <- summary(con)
  if (!(con_info$text == "binary" && con_info$`can read` == "yes")) stop("Unsuitable connection provided. read.module() requires a binary connection from which can be read.")

  mod <- new("PTModule")

  # read module name from file
  mod@name <- readBin(con, "raw", 20)

  samp.lengths <- list()

  # loop the 31 samples and read all required info,
  # except for the wave data itself, from the file
  for (i_sample in 1:length(mod@samples))
  {
    # read the sample name from the file
    mod@samples[[i_sample]]@name         <- readBin(con, "raw", 22)

    # sample length in words, multiply by two for length in bytes
    samp.lengths[[i_sample]]             <- readBin(con, "raw", 2)

    # read the fine tune value for the sample:
    mod@samples[[i_sample]]@finetune     <- readBin(con, "raw", 1)

    # read the default volume value for the sample,
    # ranging from 0x00 (min) to 0x40 (max):
    mod@samples[[i_sample]]@volume       <- readBin(con, "raw", 1)

    # the sample loop start position in words (0 when loop is off):
    mod@samples[[i_sample]]@wloopstart   <- readBin(con, "raw", 2)

    # the sample loop end position in words (1 when loop is off):
    mod@samples[[i_sample]]@wlooplen     <- readBin(con, "raw", 2)
  }
  # clean up memory:
  rm(i_sample)

  # read the length of the table specifying the order of patterns to be played:
  mod@pattern.order.length <- readBin(con, "raw", 1)

  # read a byte that may give us some info on
  # the tracker used to create the module:
  mod@tracker.byte         <- readBin(con, "raw", 1)

  # read the table specifying the order of patterns to be played:
  mod@pattern.order        <- readBin(con, "raw", 128)

  # read a tag that can help us to determine which
  # tracker was used to create the module:
  mod@tracker.flag         <- readBin(con, "raw", 4)

  # Go for ProTracker compatibility. Forget about other trackers:
  pattern_count     <- max(as.numeric(mod@pattern.order)) + 1

  # read pattern data from file. Loop through each pattern
  for (i_pattern in 1:pattern_count)
  {
    mod@patterns[[i_pattern]] <- new("PTPattern",
                                     data = matrix(readBin(con, "raw", maximumPatternTableRowCount*maximumTrackCount*4),
                                                   ncol = maximumTrackCount*4,
                                                   nrow = maximumPatternTableRowCount,
                                                   byrow = T))
  }

  # read sample wave data
  for (i_sample in 1:length(samp.lengths))
  {
    mod@samples[[i_sample]]@left <- as.integer(128 + rawToSignedInt(readBin(con, "raw", 2*rawToUnsignedInt(samp.lengths[[i_sample]]))))
  }

  # test if the end of file is reached (should be the case):
  test_byte <- readBin(con, "raw", 1)
  if (length(test_byte) > 0) warning("\nFinished reading module data from file,\nbut end of file is not reached.")

  if (!validity.PTModule(mod) & ignore.validity) stop("File is not a valid ProTracker Module.")

  return (mod)
})

setGeneric("write.module", def = function(mod, file){
  standardGeneric("write.module")
})

#' Export an PTModule object as a ProTracker module file
#'
#' Export an \code{\link{PTModule}} object as a ProTracker module file,
#' conform ProTracker 2.3A specifications.
#'
#' The routine to write ProTracker modules is based on the referenced version
#' of ProTracker 2.3A. This means that the routine may not be able to
#' write files that ar compatible with later or earlier ProTracker versions.
#'
#' @docType methods
#' @name write.module
#' @rdname write.module
#' @aliases write.module,PTModule,ANY-method
#' @param mod A valid PTModule object to be saved as a ProTracker *.mod file
#' @param file either a filename to write to, or a file connection, that
#' allows to write binary data (see \code{\link[base]{file}}).
#' @return Writes to a module file but returns nothing.
#'
#' @examples
#' ## get the PTModule object provided with the ProTrackR package
#' data("mod.intro")
#'
#' ## save the object as a valid ProTracker module file:
#' write.module(mod.intro, "intro.mod")
#'
#' ## or create the connection yourself:
#' con <- file("intro2.mod", "wb")
#' write.module(mod.intro, con)
#'
#' ## don't forget to close the connection after you're done:
#' close(con)
#' @references \url{http://wiki.multimedia.cx/index.php?title=Protracker_Module}
#'
#' \url{http://coppershade.org/articles/More!/Topics/Protracker_File_Format/}
#' @family io.operations
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("write.module", c("PTModule", "ANY"), function(mod, file)
{
  con <- file
  if (!("connection" %in% class(con))) stop ("argument con is not a file connection!")
  con_info <- summary(con)
  if (!(con_info$text == "binary" && con_info$`can write` == "yes")) stop("Unsuitable connection provided. write.module() requires a connection to which binary data can be written.")

  # read module name from file
  writeBin(mod@name, con)

  # loop the 31 samples and read all required info,
  # except for the wave data itself, from the file
  for (i_sample in 1:length(mod@samples))
  {
    # read the sample name from the file
    writeBin(mod@samples[[i_sample]]@name, con)

    # sample length in words, multiply by two for length in bytes
    writeBin(unsignedIntToRaw(length(mod@samples[[i_sample]]@left)/2, 2), con)

    # read the fine tune value for the sample:
    writeBin(mod@samples[[i_sample]]@finetune, con)

    # read the default volume value for the sample,
    # ranging from 0x00 (min) to 0x40 (max):
    writeBin(mod@samples[[i_sample]]@volume, con)

    # the sample loop start position in words (0 when loop is off):
    writeBin(mod@samples[[i_sample]]@wloopstart, con)

    # the sample loop end position in words (1 when loop is off):
    writeBin(mod@samples[[i_sample]]@wlooplen, con)
  }
  # clean up memory:
  rm(i_sample)

  # read the length of the table specifying the order of patterns to be played:
  writeBin(mod@pattern.order.length, con)

  # read a byte that may give us some info on
  # the tracker used to create the module:
  writeBin(mod@tracker.byte, con)

  # read the table specifying the order of patterns to be played:
  writeBin(mod@pattern.order, con)

  # read a tag that can help us to determine which
  # tracker was used to create the module:
  writeBin(mod@tracker.flag, con)

  # Go for ProTracker compatibility. Forget about other trackers:
  pattern_count     <- max(as.numeric(mod@pattern.order)) + 1

  # read pattern data from file. Loop through each pattern
  for (i_pattern in 1:pattern_count)
  {
    writeBin(as.vector(t(mod@patterns[[i_pattern]]@data)), con)
  }

  # read sample wave data
  for (i_sample in 1:length(mod@samples))
  {
    writeBin(signedIntToRaw(mod@samples[[i_sample]]@left - 128), con)
  }
  invisible()
})

#' @export
#' @rdname write.module
#' @aliases write.module,PTModule,character-method
setMethod("write.module", c("PTModule", "character"), function(mod, file)
{
  if (!validity.PTModule(mod)) stop("provided object is not a valid PTModule object.")
  con <- file(file, "wb")
  write.module(mod, con)
  close(con)
  invisible()
})

#' @rdname name
#' @aliases name,PTModule-method
#' @export
setMethod("name", "PTModule", function(x){
  return(rawToCharNull(x@name))
})

#' @rdname name
#' @aliases name<-,PTModule,character-method
#' @export
setReplaceMethod("name", c("PTModule", "character"), function(x, value){
  if (length(value) > 1) warning("Provided name has more than 1 element. Only first element used.")
  value <- as.character(value)[[1]]
  value <- charToRaw(value)
  if (length(value) > 20)
  {
    warning("Name is too long and will be truncated.")
    value <- value[1:20]
  }
  if (length(value) < 20) value <- c(value, raw(20 - length(value)))
  x@name <- value
  return (x)
})

setGeneric("patternOrder", function(x, full = FALSE) standardGeneric("patternOrder"))

#' Get the pattern order table
#'
#' The pattern order table is a \code{vector} of \code{numeric} indices of
#' \code{\link{PTPattern}} tables, which determines in which order the patterns
#' need to be played. This method returns this \code{vector}.
#'
#' The actual length of the \code{vector} containing the pattern order is 128
#' as per ProTracker standards. Only part of this \code{vector} is `visible'
#' and will be used to determine in which order pattern tables are to be played.
#' This method can be used to return either the visible or full (all 128) part
#' of the table. It can also be used to assign a new patter order table.
#'
#' Note that \code{\link{PTPattern}} indices start at 0, as per ProTracker
#' standards, whereas R start indices at 1. Hence, add 1 to the indices obtained
#' with \code{patternOrder}, in order to extract the correct
#' \code{\link{PTPattern}} from a \code{\link{PTModule}}.
#'
#' The maximum index plus 1 in the full pattern order table should equal
#' the number of pattern tables (see \code{\link{patternLength}}) in the
#' \code{\link{PTModule}}. Is you assign a new pattern order, with a lower
#' maximum, \code{\link{PTPattern}} objects will get lost (see also examples)!
#'
#' @rdname patternOrder
#' @name patternOrder
#' @aliases patternOrder,PTModule-method
#' @param x A \code{\link{PTModule}} object for which the pattern order table
#' needs to be returned or modified.
#' @param full A \code{logical} value indicating whether the full (\code{TRUE},
#' default), or only the visible (\code{FALSE}) part of the pattern order table
#' should be returned. This argument will also affect how new pattern order
#' tables are assigned (see \code{value}).
#' @param value A \code{numeric} \code{vector} (maximum length: 128) holding
#' \code{\link{PTPattern}} indices minus 1 for the new pattern order table.
#'
#' When \code{full = TRUE}, the \code{vector} will be padded with zeros to a
#' length of 128, and the \code{\link{patternOrderLength}} will be set to the
#' length of \code{value}. When \code{full = FALSE}, \code{value} will only
#' repplace the part of the order table up to the length of \code{value}. The
#' remainder of the table is not changed. The \code{\link{patternOrderLength}}
#' is also not modified in this case.
#' @return For \code{patternOrder}, a \code{vector} of \code{numeric}
#' \code{\link{PTPattern}} indices is returned.
#'
#' For \code{patternOrder<-}, an updated version of object \code{x} is returned,
#' in which the pattern order table is modified based on \code{value}.
#' @note The maximum number of \code{\link{PTPattern}}s cannot exceed either 64 or
#' 100 (depending on the \code{\link{trackerFlag}}). This means that values in
#' the order table should also not exceed these values minus 1.
#' @examples
#' data("mod.intro")
#'
#' ## get the visible part of the patternOrder table:
#' patternOrder(mod.intro)
#'
#' ## get the full patternOrder table:
#' patternOrder(mod.intro, full = TRUE)
#'
#' ## add 1 to get extract the right PTPattern from
#' ## mod.intro:
#' first.pattern.played <-
#'   (PTPattern(mod.intro, patternOrder(mod.intro)[1] + 1))
#'
#' ## set a different playing order:
#' patternOrder(mod.intro) <- c(0:3, 0:3, 0:3)
#'
#' ## The assignment above uses a value that
#' ## longer than the patternOrderLength.
#' ## This means that a part ends up in the
#' ## 'invisible' part of the order table:
#' patternOrder(mod.intro)
#' patternOrder(mod.intro, full = TRUE)
#'
#' ## Let's do the same assignment, but update
#' ## the visible part of the table as well:
#' patternOrder(mod.intro, full = TRUE) <- c(0:3, 0:3, 0:3)
#'
#' ## note that the maximum of the order table plus 1
#' ## equals the patternLength of mod.intro (always the case
#' ## for a valid PTModule object):
#' max(patternOrder(mod.intro, full = TRUE) + 1) ==
#'   patternLength(mod.intro)
#'
#' ## Let's do something dangerous. If the replacement
#' ## indices do not hold a maximum value that equals
#' ## the patternLength minus 1, PTPatterns will get lost,
#' ## in order to maintain the validity of mod.intro:
#' patternOrder(mod.intro) <- rep(0, 12)
#'
#' @author Pepijn de Vries
#' @family pattern.operations
#' @family module.operations
#' @export
setMethod("patternOrder", c("PTModule"), function(x, full){
  full <- as.logical(full)
  if (full)
    return(as.integer(x@pattern.order)) else
      return(as.integer(x@pattern.order[1:as.integer(x@pattern.order.length)]))
})

setGeneric("patternOrder<-", function(x, full = FALSE, value) standardGeneric("patternOrder<-"))

#' @rdname patternOrder
#' @name patternOrder<-
#' @aliases patternOrder<-,PTModule,ANY,numeric-method
#' @export
setReplaceMethod("patternOrder", c("PTModule", "ANY", "numeric"), function(x, full, value){
  full <- as.logical(full)
  value <- abs(as.integer(as.vector(value)))
  if (any(is.na(value))) stop ("NAs are not allowed in replacement!")
  if (length(value) > 128)
  {
    warning("Replacement has more than 128 elements. Only the first 128 elements used")
  } else if(!full)
  {
    value <- c(value, as.numeric(x@pattern.order)[-(1:length(value))])
  }
  max.patterns <- ifelse(all(x@tracker.flag    == charToRaw("M!K!")), 100, 64)
  if (any(value >= max.patterns))
  {
    value[value >= max.patterns] <- max.patterns - 1
    warning(paste("Replacement contained one or more indices above the ",
                  "maximum number of patterns allowed (", max.patterns,
                  "). Ceiling is applied to the values.", sep = ""))
  }
  current <- as.integer(x@pattern.order[1:as.integer(x@pattern.order.length)])
  if (max(current) < max(value))
  {
    for (i in (max(current) + 2):(max(value)+ 1))
    {
      x@patterns[[i]] <- new("PTPattern")
    }
  } else if(max(current) > max(value))
  {
    warning(paste("Replacement pattern order contain less patterns ",
                  "than the original list. These patterns are now lost!",
                  sep = ""))
    for (i in (max(current)+ 1):(max(value) + 2))
    {
      x@patterns[[i]] <- NULL
    }
  }
  x@pattern.order <- as.raw(c(value, rep(0, 128 - length(value))))
  if (full) x@pattern.order.length <- as.raw(length(value))
  return (x)
})

setGeneric("patternLength", function(x) standardGeneric("patternLength"))

#' Get the number of PTPattern tables in a PTModule
#'
#' Get the number of \code{\link{PTPattern}} tables in a \code{\link{PTModule}}
#' object.
#'
#' The number of \code{\link{PTPattern}} tables in a \code{\link{PTModule}}
#' object should range from 1 up to either 64 or 100. The maximum depends on the
#' \code{\link{trackerFlag}} of the \code{\link{PTModule}} object.
#' @rdname patternLength
#' @name patternLength
#' @aliases patternLength,PTModule-method
#' @param x A \code{\link{PTModule}} object for which the number of
#' \code{\link{PTPattern}} tables need to be returned.
#' @return Returns a \code{numeric} value representing the number of
#' \code{\link{PTPattern}} tables in object \code{x}.
#' @examples
#' data("mod.intro")
#'
#' ## Get the number of pattern tables in mod.intro:
#' patternLength(mod.intro)
#'
#' @family pattern.operations
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("patternLength", "PTModule", function(x){
  return(length(x@patterns))
})

setGeneric("patternOrderLength", function(x) standardGeneric("patternOrderLength"))
setGeneric("patternOrderLength<-", function(x, value) standardGeneric("patternOrderLength<-"))

#' Get the length of the pattern order table
#'
#' The pattern order table is a \code{vector} of \code{numeric} indices of
#' \code{PTPattern} tables, which determines in which order the patterns
#' need to be played. This method returns the visible length of this
#' \code{vector}.
#'
#' The actual length of the \code{vector} containing the pattern order is 128
#' as per ProTracker standards. Only part of this \code{vector} is `visible'
#' and will be used to determine in which order pattern tables are to be played.
#' The length returned by this method is the length of this visible part of the
#' pattern order table. The length of this visible part can also be set with this
#' method.
#'
#' @rdname patternOrderLength
#' @name patternOrderLength
#' @aliases patternOrderLength,PTModule-method
#' @param x A \code{\link{PTModule}} object for which the length of the
#' visible part of the pattern order table is to be returned.
#' @param value A \code{numeric} value which is to be used to set the visible
#' length of the pattern order table.
#' @return For \code{patternOrderLength} the visible length of the pattern
#' order table of \code{\link{PTModule}} \code{x} is returned as a \code{numeric}
#' value, ranging from 1 up to 128.
#'
#' For \code{patternOrderLength<-} an updated version of object \code{x} is
#' returned, in which the visible length of the pattern order table is set
#' to \code{value}. Note that this does not change the pattern order table
#' itself, only which part is `visible'.
#' @examples
#' data("mod.intro")
#'
#' ## get the length of the pattern order table:
#' patternOrderLength(mod.intro)
#'
#' ## set the length of the pattern order table to 1:
#' patternOrderLength(mod.intro) <- 1
#'
#' ## note that the pattern order table remained intact:
#' patternOrder(mod.intro, full = TRUE)
#'
#' @family pattern.operations
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("patternOrderLength", "PTModule", function(x){
  return(as.numeric(x@pattern.order.length))
})

#' @rdname patternOrderLength
#' @name patternOrderLength<-
#' @aliases patternOrderLength<-,PTModule,numeric-method
#' @export
setReplaceMethod("patternOrderLength", c("PTModule", "numeric"), function(x, value){
  value <- as.integer(value)[[1]]
  if (value < 1) stop ("Pattern order length should be at least 1.")
  if (value > 128)
  {
    value <- 128
    warning("Provided length for pattern.order.length is larger than 128. It is set to the limit of 128.")
  }
  x@pattern.order.length <- as.raw(value)
  return(x)
})

setGeneric("trackerFlag", function(x) standardGeneric("trackerFlag"))
setGeneric("trackerFlag<-", function(x, value = c("M.K.", "M!K!")) standardGeneric("trackerFlag<-"))

#' Tracker flag indicating version compatibility
#'
#' Method to obtain a tracker flag, which indicates the version compatibility
#' of a ProTracker module (\code{\link{PTModule}} object).
#'
#' ProTrackR supports two tracker flags: "\code{M.K.}" and "\code{M!K!}". M.K.
#' are presumably the initials of programmers Mahony and Kaktus, unfortunately
#' documentation on this matter is ambiguous. In any case, modules with the
#' flag "\code{M.K.}" can hold up to 64 patterns, whereas modules with the flag
#' "\code{M!K!}" can hold up to 100 patterns. Use this method to obtain or
#' replace the tracker flag of a \code{\link{PTModule}}.
#' @rdname trackerFlag
#' @name trackerFlag
#' @aliases trackerFlag,PTModule-method
#' @param x A \code{\link{PTModule}} object for which the flag needs to
#' returned or replaced.
#' @param value A \code{character} string representing the tracker flag with which
#' that of object \code{x} needs to be replaced with. Should either be "\code{M.K.}"
#' or "\code{M!K!}". Note that if a current flag "\code{M!K!}" is
#' replaced by "\code{M.K.}", \code{\link{PTPattern}}s may get lost as the
#' latter supports less patterns.
#' @return For \code{trackerFlag}, the tracker flag of object \code{x} is returned.
#'
#' For \code{trackerFlag<-}, a copy of object \code{x} with an updated tracker
#' flag is returned.
#' @examples
#' data("mod.intro")
#'
#' ## the current trackerFlag of mod.intro is "M.K.",
#' ## meaning that it can hold a maximum of 64 patterns:
#' trackerFlag(mod.intro)
#'
#' patternOrder(mod.intro, full = TRUE) <- 0:63
#'
#' ## If we upgrade the trackerFlag of mod.intro to "M!K!"
#' ## it can hold a maximum of 100 patterns!:
#' trackerFlag(mod.intro) <- "M!K!"
#'
#' patternOrder(mod.intro, full = TRUE) <- 0:99
#'
#' ## Now let's do something dangerous:
#' ## current flag is "M!K!", by setting it
#' ## back to "M.K.", patterns 65:100 are lost...
#' trackerFlag(mod.intro) <- "M.K."
#'
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("trackerFlag", "PTModule", function(x){
  return(rawToCharNull(x@tracker.flag))
})

#' @rdname trackerFlag
#' @name trackerFlag<-
#' @aliases trackerFlag<-,PTModule-method
#' @export
setReplaceMethod("trackerFlag", c("PTModule"), function(x, value){
  if (match.arg(value) == "M.K.")
  {
    x@tracker.flag <- charToRaw("M.K.")
    if (length(x@patterns) > 64)
    {
      warning("PTModule holds more patterns than M.K. flag allows. Higher patterns are removed!")
      for (i_pat in length(x@patterns):65)
      {
        x@patterns[[i_pat]] <- NULL
      }
      x@pattern.order[as.integer(x@pattern.order) > 63] <- as.raw(0x00)
    }
    return (x)
  }
  if (match.arg(value) == "M!K!")
  {
    x@tracker.flag <- charToRaw("M!K!")
    return (x)
  }
})

setGeneric("deletePattern", function(x, index){
  standardGeneric("deletePattern")
})

#' Remove a PTPattern table from a PTModule object
#'
#' This method removes a \code{\link{PTPattern}} from a
#' \code{\link{PTModule}} object and updates the
#' \code{\link{patternOrder}} table accordingly.
#'
#' This method safely removes a \code{\link{PTPattern}} from a
#' \code{\link{PTModule}} object, guarding the validity of the
#' \code{\link{PTModule}} object. It therefore also updates
#' the \code{\link{patternOrder}} table, by renumbering the indices
#' listed there. The index of the removed object is replaced with a zero
#' in the \code{\link{patternOrder}} table.
#'
#' @note As per ProTracker specification, the pattern indices
#' stored in the \code{\link{PTModule}} and obtained with
#' \code{\link{patternOrder}} start at 0. Whereas R starts indexing at 1.
#' Beware of this discrepancy.
#' @rdname deletePattern
#' @name deletePattern
#' @aliases deletePattern,PTModule,numeric-method
#' @param x A \code{\link{PTModule}} from which a
#' \code{\link{PTPattern}} needs to be removed.
#' @param index A \code{numeric} index of the \code{\link{PTPattern}}
#' table that needs to be removed. The index should be between 1 and
#' \code{\link{patternLength}}. It's not possible to delete multiple
#' patterns simultaneously with this method. A \code{\link{PTModule}}
#' should always hold at least 1 pattern table, therefore, the last
#' \code{\link{PTPattern}} table cannot be deleted.
#' @return Returns a \code{\link{PTModule}} from which the selected
#' \code{\link{PTPattern}} is deleted.
#' @examples
#' data("mod.intro")
#' print(mod.intro)
#'
#' ## delete pattern #2 from mod.intro:
#'
#' mod.intro <- deletePattern(mod.intro, 2)
#' print(mod.intro)
#'
#' @family pattern.operations
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("deletePattern", c("PTModule", "numeric"), function(x, index){
  index <- as.integer(index)
  if (length(index) != 1 ) stop("Only a single index should be used.")
  if (index < 1) stop("Invalid index.")
  if (index > length(x@patterns)) stop("Index out of range.")
  if (length(x@patterns) < 2) stop("You can't delete the last pattern table!")
  p.order <- patternOrder(x, T)
  p.order[p.order == (index - 1)] <- 0
  p.order[p.order > (index - 1)] <- p.order[p.order > (index - 1)] - 1
  x@patterns[index] <- NULL
  x@pattern.order <- as.raw(p.order)
  return(x)
})

setGeneric("appendPattern", function(x, pattern){
  standardGeneric("appendPattern")
})

#' Append a PTPattern to a PTModule
#'
#' Appends a specified \code{\link{PTPattern}} to a
#' \code{\link{PTModule}}.
#'
#' Depending on the \code{\link{trackerFlag}}, a ProTracker module can hold
#' either 64 or 100 pattern tables. As long as the number of pattern tables
#' is below this maximum, new pattern tables can be added to the module with
#' this function.
#'
#' The \code{\link{patternOrder}} table should hold the maximum index of the
#' available pattern tables in a module, otherwise, the module is not valid.
#' As the maximum index increases, by appending a pattern table, the
#' \code{\link{patternOrder}} table should be updated. The
#' \code{\link{appendPattern}} method does this automatically, by replacing the first
#' non-unique index in the order table, outside the current order table's length,
#' with the new maximum index. If this is not possible, the highest element
#' in the order table is set to hold the maximum index.
#'
#' @note As per ProTracker specification, the pattern indices
#' stored in the \code{\link{PTModule}} and obtained with
#' \code{\link{patternOrder}} start at 0. Whereas R starts indexing at 1.
#' Beware of this discrepancy.
#' @rdname appendPattern
#' @name appendPattern
#' @aliases appendPattern,PTModule,PTPattern-method
#' @param x A \code{\link{PTModule}} object to which a
#' \code{\link{PTPattern}} is to be appended.
#' @param pattern A \code{\link{PTPattern}} object which is
#' to be appended to the \code{\link{PTModule}} \code{x}.
#' @return Returns a \code{\link{PTModule}}, to which the
#' \code{\link{PTPattern}} is appended.
#' @examples
#' data("mod.intro")
#'
#' ## append an empty pattern to mod.intro
#'
#' mod.intro <- appendPattern(mod.intro, new("PTPattern"))
#'
#' ## append a copy of pattern # 1 (this is pattern #0 in the
#' ## patternOrder table) to mod.intro
#'
#' mod.intro <- appendPattern(mod.intro, PTPattern(mod.intro, 1))
#'
#' @family pattern.operations
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("appendPattern", c("PTModule", "PTPattern"), function(x, pattern){
  max.patterns <- ifelse(all(x@tracker.flag    == charToRaw("M!K!")), 100, 64)
  if (length(x@patterns) > max.patterns) stop ("Can't insert pattern. Module already holds maximum number of patterns.")
  p.order     <- patternOrder(x, full = T)
  p.order.len <- patternOrderLength(x)
  x@patterns[[length(x@patterns) + 1]] <- pattern
  if (p.order.len == 128){
    warning("No place in the pattern order table to include the new pattern.
            Last item in the order table is replaced!")
    p.order[128] <- length(x@patterns) - 1
  } else
  {
    index <- which(duplicated(p.order) & ((1:128) > p.order.len))
    if (length(index) == 0) index <- 128
    index <- min(index)
    p.order[index] <- length(x@patterns) - 1
  }
  x@pattern.order <- as.raw(p.order)
  return(x)
})

setGeneric("moduleSize", function(x) standardGeneric("moduleSize"))

#' Get module file size
#'
#' Get the file size in bytes of a \code{\link{PTModule}} object, when it is to be saved
#' as an original module file with \code{\link{write.module}}.
#'
#' The ProTracker module has a 1084 byte sized header containing all (meta)
#' information on the patterns, their order and the audio samples. Each pattern
#' holds exactly 1 Kb of information and the length of the audio samples corresponds
#' with the size in bytes, as they are of 8 bit quality in mono. This function
#' calculates the file size of the \code{\link{PTModule}} object when it is to
#' be saved with \code{\link{write.module}}.
#' @rdname moduleSize
#' @name moduleSize
#' @aliases moduleSize,PTModule-method
#' @param x A \code{\link{PTModule}} object for which the file size is
#' to be calculated.
#' @return Returns potential uncompressed module file size in bytes represented
#' by a number of class \code{object_size}.
#' @examples
#' ## Calculate the file size for the example module 'mod.intro':
#'
#' data("mod.intro")
#' moduleSize(mod.intro)
#'
#' ## Note that this is not the same as the size the object
#' ## requires in R working memory:
#'
#' object.size(mod.intro)
#'
#' ## In working memory it takes more memory to store the module, than in a
#' ## file. This is because the S4 structure of the object consumes some
#' ## memory. In addition, samples are of 8 bit quality, corresponding with
#' ## a byte per sample. In the PTSample object it is stored as a
#' ## vector of integer values. In R, integer values are 32 bit, which
#' ## costs 4 times as much memory as the original 8 bit.
#'
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("moduleSize", "PTModule", function(x){
  header_size   <- sum(c(20, 1, 1, 128, 4))
  sample_header <- sum(c(22, 2, 1, 1, 2, 2))
  pattern_size  <- 64*4*4
  pat.len  <- length(x@patterns)
  samp.len <- sum(unlist(lapply(x@samples, function(x) length(x@left))))
  result <- header_size + sample_header*31 +
    pat.len*pattern_size + samp.len
  class(result) <- "object_size"
  return(result)
})

setGeneric("clearSong", function(mod) standardGeneric("clearSong"))

#' Clear all pattern info from module
#'
#' Remove all patterns (\code{\link{PTPattern}}) and \code{\link{patternOrder}}
#' info from a \code{\link{PTModule}} object.
#'
#' Conform the original ProTracker, this method removes all patterns
#' (\code{\link{PTPattern}}) and \code{\link{patternOrder}}
#' info from a module. You keep the audio \code{\link{PTSample}}s.
#'
#' @rdname clearSong
#' @name clearSong
#' @aliases clearSong,PTModule-method
#' @param mod A \code{\link{PTModule}} object from which all pattern (order)
#' info needs to be removed.
#' @return Returns a copy of of object \code{mod} in which all pattern (order)
#' info is removed.
#' @examples
#' data(mod.intro)
#'
#' ## 'clear.mod' is a copy of 'mod.intro' without the
#' ## pattern (order) info. It still has the audio samples.
#' clear.mod <- clearSong(mod.intro)
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("clearSong", "PTModule", function(mod){
  mod@patterns <- list(new("PTPattern"))
  suppressWarnings(patternOrder(mod, T) <- 0)
  return(mod)
})

setGeneric("clearSamples", function(mod) standardGeneric("clearSamples"))

#' Clear all samples from module
#'
#' Remove all \code{\link{PTSample}}s from a \code{\link{PTModule}} object.
#'
#' Conform the original ProTracker, this method removes all patterns
#' \code{\link{PTSample}}s from a module. You keep all patterns
#' (\code{\link{PTPattern}}) and \code{\link{patternOrder}} info.
#'
#' @rdname clearSamples
#' @name clearSamples
#' @aliases clearSamples,PTModule-method
#' @param mod A \code{\link{PTModule}} object from which all samples needs
#' to be removed.
#' @return Returns a copy of of object \code{mod} in which all samples are removed.
#' @examples
#' data(mod.intro)
#'
#' ## 'clear.mod' is a copy of 'mod.intro' without the
#' ## samples. It still holds all pattern tables and
#' ## pattern order info.
#' clear.mod <- clearSamples(mod.intro)
#' @family module.operations
#' @author Pepijn de Vries
#' @export
setMethod("clearSamples", "PTModule", function(mod){
  mod@samples <- lapply(as.list(1:31), function(x) new("PTSample"))
  return(mod)
})
