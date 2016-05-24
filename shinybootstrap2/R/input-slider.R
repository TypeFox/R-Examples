#' Slider Input Widget
#'
#' Constructs a slider widget to select a numeric value from a range.
#'
#' @param inputId Specifies the \code{input} slot that will be used to access
#'   the value.
#' @param label A descriptive label to be displayed with the widget, or
#'   \code{NULL}.
#' @param min The minimum value (inclusive) that can be selected.
#' @param max The maximum value (inclusive) that can be selected.
#' @param value The initial value of the slider. A numeric vector of length
#'   one will create a regular slider; a numeric vector of length two will
#'   create a double-ended range slider. A warning will be issued if the
#'   value doesn't fit between \code{min} and \code{max}.
#' @param step Specifies the interval between each selectable value on the
#'   slider (\code{NULL} means no restriction).
#' @param round \code{TRUE} to round all values to the nearest integer;
#'   \code{FALSE} if no rounding is desired; or an integer to round to that
#'   number of digits (for example, 1 will round to the nearest 10, and -2 will
#'   round to the nearest .01). Any rounding will be applied after snapping to
#'   the nearest step.
#' @param format Customize format values in slider labels. See
#'   \url{https://code.google.com/p/jquery-numberformatter/} for syntax
#'   details.
#' @param locale The locale to be used when applying \code{format}. See details.
#' @param ticks \code{FALSE} to hide tick marks, \code{TRUE} to show them
#'   according to some simple heuristics.
#' @param animate \code{TRUE} to show simple animation controls with default
#'   settings; \code{FALSE} not to; or a custom settings list, such as those
#'   created using \code{animationOptions}.
#' @inheritParams selectizeInput
#' @family input elements
#' @seealso \code{\link{updateSliderInput}}
#'
#' @details
#'
#' Valid values for \code{locale} are: \tabular{ll}{ Arab Emirates \tab "ae" \cr
#' Australia \tab "au" \cr Austria \tab "at" \cr Brazil \tab "br" \cr Canada
#' \tab "ca" \cr China \tab "cn" \cr Czech \tab "cz" \cr Denmark \tab "dk" \cr
#' Egypt \tab "eg" \cr Finland \tab "fi" \cr France  \tab "fr" \cr Germany \tab
#' "de" \cr Greece \tab "gr" \cr Great Britain \tab "gb" \cr Hong Kong \tab "hk"
#' \cr India \tab "in" \cr Israel \tab "il" \cr Japan \tab "jp" \cr Russia \tab
#' "ru" \cr South Korea \tab "kr" \cr Spain \tab "es" \cr Sweden \tab "se" \cr
#' Switzerland \tab "ch" \cr Taiwan \tab "tw" \cr Thailand \tab "th" \cr United
#' States \tab "us" \cr Vietnam \tab "vn" \cr }
#'
#' @export
sliderInput <- function(inputId, label, min, max, value, step = NULL,
                        round=FALSE, format='#,##0.#####', locale='us',
                        ticks=TRUE, animate=FALSE, width=NULL) {

  if (identical(animate, TRUE))
    animate <- animationOptions()

  if (!is.null(animate) && !identical(animate, FALSE)) {
    if (is.null(animate$playButton))
      animate$playButton <- tags$i(class='icon-play')
    if (is.null(animate$pauseButton))
      animate$pauseButton <- tags$i(class='icon-pause')
  }

  # build slider
  sliderTag <- slider(inputId, min=min, max=max, value=value, step=step,
    round=round, locale=locale, format=format, ticks=ticks, animate=animate,
    width=width)

  if (is.null(label)) {
    sliderTag
  } else {
    tags$div(
      controlLabel(inputId, label),
      sliderTag
    )
  }
}

hasDecimals <- function(value) {
  truncatedValue <- round(value)
  return (!identical(value, truncatedValue))
}

#' @rdname sliderInput
#'
#' @param interval The interval, in milliseconds, between each animation step.
#' @param loop \code{TRUE} to automatically restart the animation when it
#'   reaches the end.
#' @param playButton Specifies the appearance of the play button. Valid values
#'   are a one-element character vector (for a simple text label), an HTML tag
#'   or list of tags (using \code{\link{tag}} and friends), or raw HTML (using
#'   \code{\link{HTML}}).
#' @param pauseButton Similar to \code{playButton}, but for the pause button.
#'
#' @export
animationOptions <- function(interval=1000,
                             loop=FALSE,
                             playButton=NULL,
                             pauseButton=NULL) {
  list(interval=interval,
       loop=loop,
       playButton=playButton,
       pauseButton=pauseButton)
}

# Create a new slider control (list of slider input element and the script
# tag used to configure it). This is a lower level control that should
# be wrapped in an "input" construct (e.g. sliderInput in bootstrap.R)
#
# this is a wrapper for: https://github.com/egorkhmelev/jslider
# (www/shared/slider contains js, css, and img dependencies)
slider <- function(inputId, min, max, value, step = NULL, ...,
                   round=FALSE, format='#,##0.#####', locale='us',
                   ticks=TRUE, animate=FALSE, width=NULL) {
  # validate inputId
  inputId <- as.character(inputId)
  if (!is.character(inputId))
    stop("inputId not specified")

  # validate numeric inputs
  if (!is.numeric(value) || !is.numeric(min) || !is.numeric(max))
    stop("min, max, and value must all be numeric values")
  else if (min(value) < min)
    stop(paste("slider initial value", value,
               "is less than the specified minimum"))
  else if (max(value) > max)
    stop(paste("slider initial value", value,
               "is greater than the specified maximum"))
  else if (min > max)
    stop(paste("slider maximum is greater than minimum"))
  else if (!is.null(step)) {
    if (!is.numeric(step))
      stop("step is not a numeric value")
    if (step > (max - min))
      stop("step is greater than range")
  }

  # step
  range <- max - min
  if (is.null(step)) {
    # short range or decimals means continuous decimal
    if (range < 2 || hasDecimals(min) || hasDecimals(max))
      step <- range / 250 # ~ one step per pixel
    else
      step = 1
  }

  # Default state is to not have ticks
  if (identical(ticks, TRUE)) {
    # Automatic ticks
    tickCount <- (range / step) + 1
    if (tickCount <= 26)
      ticks <- paste(rep('|', floor(tickCount)), collapse=';')
    else {
      ticks <- NULL
#       # This is a smarter auto-tick algorithm, but to be truly useful
#       # we need jslider to be able to space ticks irregularly
#       tickSize <- 10^(floor(log10(range/0.39)))
#       if ((range / tickSize) == floor(range / tickSize)) {
#         ticks <- paste(rep('|', (range / tickSize) + 1), collapse=';')
#       }
#       else {
#         ticks <- NULL
#       }
    }
  }
  else if (is.numeric(ticks) && length(ticks) == 1) {
    # Use n ticks
    ticks <- paste(rep('|', ticks), collapse=';')
  }
  else if (length(ticks) > 1 && (is.numeric(ticks) || is.character(ticks))) {
    # Explicit ticks
    ticks <- paste(ticks, collapse=';')
  }
  else {
    ticks <- NULL
  }

  # build slider
  dep <- htmlDependency("jslider", "1",
    c(file = system.file("www/jslider", package = "shinybootstrap2")),
    script = c("js/jquery.slider.min.js", "jslider-shiny.js"),
    stylesheet = "css/jquery.slider.min.css"
  )
  sliderFragment <- list(
    attachDependencies(
      tags$input(
        id=inputId, type="slider",
        name=inputId, value=paste(value, collapse=';'), class="jslider",
        'data-from'=min, 'data-to'=max, 'data-step'=step,
        'data-skin'='plastic', 'data-round'=round, 'data-locale'=locale,
        'data-format'=format, 'data-scale'=ticks,
        'data-smooth'=FALSE,
        'data-width'=validateCssUnit(width)
      ),
      dep
    )
  )

  if (identical(animate, TRUE))
    animate <- animationOptions()

  if (!is.null(animate) && !identical(animate, FALSE)) {
    if (is.null(animate$playButton))
      animate$playButton <- 'Play'
    if (is.null(animate$pauseButton))
      animate$pauseButton <- 'Pause'

    sliderFragment[[length(sliderFragment)+1]] <-
      tags$div(class='slider-animate-container',
               tags$a(href='#',
                      class='slider-animate-button',
                      'data-target-id'=inputId,
                      'data-interval'=animate$interval,
                      'data-loop'=animate$loop,
                      tags$span(class='play', animate$playButton),
                      tags$span(class='pause', animate$pauseButton)))
  }

  return(tagList(sliderFragment))
}
