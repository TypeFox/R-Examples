#' ProTracker Period Table
#'
#' Table of ProTracker period values and corresponding, octave, tone and fine tune
#'
#' Table of ProTracker period values used in calculating the
#' playback sampling rate of samples for specific tones.
#' These are the values that are actually used by ProTracker,
#' they cannot be calculated directly due to
#' undocumented rounding inconsistencies. This lookup table is therefore
#' a requirement.
#'
#' @docType data
#' @name period_table
#' @format a \code{data.frame} with fourteen columns:
#' \itemize{
#'   \item{The column named `octave': \code{integer} value [1,3]}
#'   \item{The column named `finetune': \code{integer} value [-8, 7] used to tune a sample}
#'   \item{The columns named `C-' to `B-': represent the twelve (semi)tones.
#'   The values in these columns are the period values for the corresponding
#'   tone, octave and finetune.}
#' }
#' @family period.operations
#' @examples
#' data("period_table")
NULL

#' Paula clock table
#'
#' Table that provides audio output frequencies for the Commodore Amiga
#' original chipset.
#'
#' Paula was one of the custom chips on the original Commodore Amiga. This chip
#' was dedicated (amongst other tasks) to controlling audio playback. The
#' chip's output rate depended on the video mode used:
#' either `\href{https://en.wikipedia.org/wiki/PAL}{PAL}'
#' or `\href{https://en.wikipedia.org/wiki/NTSC}{NTSC}'. This table provides the
#' output rate for both video modes that can be used in calculating sample rates.
#'
#' @docType data
#' @name paula_clock
#' @format a \code{data.frame} with two columns:
#' \itemize{
#'   \item{`frequency' A \code{numeric} value representing Paula's output rate in Hz.}
#'   \item{`video' A \code{character} string representing the two video modes.}
#' }
#' @references \url{https://en.wikipedia.org/wiki/Original_Chip_Set#Paula}
#' @examples
#' data("paula_clock")
NULL

#' ProTracker Funk Table
#'
#' Small list of numbers used by an obscure audio effect in ProTracker
#'
#' This dataset is included for completeness sake. It is not yet used by any
#' class, method or function in the \code{\link{ProTrackR}} package. It may
#' very well be obsolete for recent ProTracker versions.
#' @docType data
#' @name funk_table
#' @format A \code{numeric} \code{vector} of length 16 holding values to be
#' used in ProTracker funk repeat effects.
#' @references \url{http://fossies.org/linux/uade/amigasrc/players/tracker/eagleplayers/mod32_protracker/PTK_versions.txt}
#' @examples data("funk_table")
NULL

#' Example of a PTModule object
#'
#' A \code{\link{PTModule}} object included in the package as example.
#'
#' This PTModule object is based on an original ProTracker module file
#' I've composed in the late nineteen nineties. It is used as example for many
#' of the \code{\link{ProTrackR}} methods and you can use it to test your own
#' code. It can also be exported back to the original ProTracker module file
#' by using \code{\link{write.module}}.
#' @docType data
#' @name mod.intro
#' @format A \code{\link{PTModule}} object containing 4
#' \code{\link{PTSample}} objects (and 27 empty \code{PTSample}
#' objects, adding up to the 31 samples a \code{PTModule} should hold) and 4
#' \code{\link{PTPattern}} objects.
#' @examples
#' data("mod.intro")
#' print(mod.intro)
#' plot(mod.intro)
#'
#' \dontrun{
#' playSample(mod.intro)
#'
#' ## Save as an original module file,
#' ## which can be played with ProTracker (or several modern audio players):
#' write.module(mod.intro, "intro.mod")
#' }
#' @author Pepijn de Vries
NULL
