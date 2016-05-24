#' Detect Intervals in a Ride.
#'
#' Section a ride file according to power output.
#'
#' Often a ride will contain intervals/efforts that are not in any way marked in
#' the device data (e.g. as "laps"). Using changepoint analysis, it is possible
#' to retrospectively identify these efforts. This is contingent on supplying
#' the number of changepoints to the underlying algorithm, simplified here as a
#' \code{"sections"} argument.
#'
#' For example, if there are two efforts amidst a
#' ride, this means we are looking to identify 5 \emph{sections} (i.e.
#' neutral-effort-neutral-effort-neutral). See examples.
#'
#' Depends on the package \code{"changepoint"}.
#'
#' @param data a \strong{formatted} dataset produced by \code{read*()}.
#' @param sections how many sections should be identified? Includes stoppages.
#' @param plot logical; if \code{TRUE}, graphically displays the resultant sections.
#' @param ... graphical parameters to be passed to \code{par()}. Ignored if
#'   \code{plot = FALSE}.
#'
#' @return if \code{plot = TRUE} nothing is returned. If \code{plot = FALSE}
#'   (default) a vector of section "levels" is returned.
#'
#' @examples
#' data(intervaldata)
#'
#' ## "intervaldata" is a ride that includes two efforts (2 & 5 minutes) and a cafe
#' ## stop. The efforts are marked in the lap column, which we can use as a
#' ## criterion.
#'
#' with(intervaldata, tapply(X = delta.t, INDEX = lap, sum)) / 60  # Minutes.
#'
#' ## The above shows the efforts were laps two and four. What was the power?
#' with(intervaldata, tapply(X = power.W, INDEX = lap, mean))[c(2, 4)]
#'
#' ## And for the sake of example, some other summary metrics...
#' l <- split(intervaldata, intervaldata$lap)
#' names(l) <- paste("Lap", names(l))  # Pretty names.
#' vapply(l, FUN.VALUE = numeric(3), FUN = function(x)
#'   c(t.min = ride_time(x$timer.s) / 60, NP = NP(x), TSS = TSS(x)))
#'
#' ## Could we have gotten the same information without the lap column?
#' ## Two efforts and a cafe stop == 7 sections.
#' interval_detect(intervaldata, sections = 7, plot = TRUE)
#'
#' ## An overzealous start to the first effort is being treated as a seperate section,
#' ## so let's allow for an extra section...
#' interval_detect(intervaldata, sections = 8, plot = TRUE)
#'
#' ## Looks okay, so save the output and combine the second and third sections.
#' intervaldata$intv <- interval_detect(intervaldata, sections = 8, plot = FALSE)
#' intervaldata$intv[intervaldata$intv == 3] <- 2
#'
#' ## Are the timings as expected?
#' with(intervaldata, tapply(X = delta.t, INDEX = intv, sum)) / 60  # Minutes.
#'
#' ## Close enough!
#'
#' i <- split(intervaldata, intervaldata$intv)
#' names(i) <- paste("Interval", seq_along(i))  # Pretty names.
#' toplot <- vapply(i, FUN.VALUE = numeric(3), FUN = function(x)
#'   c(t.min = ride_time(x$timer.s) / 60, NP = NP(x), TSS = TSS(x)))
#'
#' print(toplot)
#'
#' par(mfrow = c(3, 1))
#' mapply(function(r, ylab) barplot(
#'   toplot[r, c(1:3, 5:7)], names.arg = seq_along(toplot[r, c(1:3, 5:7)]),
#'   xlab = "Section", ylab = ylab),
#'   r = 1:3, ylab = c("Ride time (minutes)", "NP", "TSS"))
#'
#' @export
interval_detect <- function(data, sections, plot = FALSE, ...)
  UseMethod("interval_detect", data)
#' @export
interval_detect.default <- function(data, sections, plot = FALSE, ...)
  format_error()
#' @export
interval_detect.cycleRdata <- function(data, sections, plot = FALSE, ...) {
  stopifnot(sections >= 2)
  data  <- data[, c("timer.s", "power.W")]  # Isolate relevant columns.
  data  <- expand_stops(data)
  if (!is.null(attr(data, "new")))
    data[attr(data, "new"), "power.W"] <- 0

  minseglen <- 10 / attr(data, "deltat")  # 10 seconds.

  cpts <- suppressWarnings(changepoint::cpt.mean(
    data[, "power.W"], method = "BinSeg", minseglen = minseglen,
    Q = (sections - 1),  # Number of changepoints.
    class = FALSE        # Dont return "cpt" object.
  ))
  # Assemble "section" values.
  data$section <- 0
  data$section[c(1, cpts[-length(cpts)])] <- 1
  data$section <- cumsum(data$section)

  if (plot) {  # Doesn't return anything.
    opar <- par(no.readonly = TRUE); on.exit(par(opar))
    par(...)
    # Underlying power data.
    with(data, plot(
      x = timer.s, y = power.W, type = "l",
      col = "#333333", lwd = 0.5,
      xlab = "Time (sec)", ylab = "Power (watts)"
    ))
    # Overlay section lines (mean powers).
    data$avepower <- ave(data$power.W, data$section)
    s <- unique(data$section)
    for (i in s)
      with(data[data$section == i, ],
           lines(x = timer.s, y = avepower, lwd = 3, col = base_pal(i, s)))

    invisible()
  } else
    data[attr(data, "wo_expand"), "section"]
}
