#' Read cycling device data.
#'
#' Read data from a cycling head unit into the R environment; optionally
#' formatting it for use with other functions in this package. Critical power
#' and session RPE metrics can also be associated with the data and used by
#' other functions (e.g. \link{summary.cycleRdata}).
#'
#' Note that most functions within this package depend on imported data being
#' formatted; i.e. \code{read*("file_path", format = TRUE)}. Hence, unless the
#' raw data is of particular interest and/or the user wants to process it
#' manually, the format argument should be TRUE (default). When working with a
#' formatted dataset, do not change existing column names. The formatted data
#' structure is described in detail in \link{ridedata}.
#'
#' Garmin .fit file data is parsed with the java command line tool provided in
#' the \href{http://www.thisisant.com/resources/fit}{FIT SDK}. The latest source
#' code and licensing information can be found at the previous link.
#'
#' SRM device files (.srm) are also parsed at the command line, provided
#' \href{http://www.zuto.de/project/srmio/}{Rainer Clasen's srmio library} is
#' installed and available. The associated GitHub repo' can be found
#' \href{https://github.com/rclasen/srmio}{here}.
#'
#' @param file character; path to the file.
#' @param format logical; should data be formatted?
#' @param CP,sRPE optional; critical power and session RPE values to be
#'   associated with the data. Ignored if \code{format = FALSE}.
#'
#' @examples \dontrun{
#' fl  <- system.file("extdata/example_files.tar.gz",
#'                    package = "cycleRtools")
#' fls <- untar(fl, list = TRUE)
#' untar(fl)  # Extract to working directory.
#'
#' dat <- lapply(fls, read_ride, format = TRUE, CP = 300, sRPE = 5)
#'
#' file.remove(fls)
#' }
#'
#' @return a data frame object.
#'
#' @describeIn read_ride A wrapper for read_* functions that chooses the
#'   appropriate function based on file extension.
#'
#' @export
#'
#' @useDynLib cycleRtools
#' @importFrom Rcpp sourceCpp
#' @import xml2
#' @import stats
#' @import graphics
#' @import grDevices
#' @import utils
read_ride <- function(file = file.choose(), format = TRUE, CP = NULL, sRPE = NULL) {
  data <- switch(tools::file_ext(file),
                 "fit"  = read_fit(file, format, CP, sRPE),
                 "srm"  = read_srm(file, format, CP, sRPE),
                 "pwx"  = read_pwx(file, format, CP, sRPE),
                 "tcx"  = read_tcx(file, format, CP, sRPE),
                 {message("Unrecognised file extension, returning NULL."); NULL})
  data
}
#' @describeIn read_ride Read a Garmin (Ltd) device .fit file. This invokes
#'   \code{\link[base]{system2}} to execute the FitCSVTool.jar command line tool
#'   (see \href{http://www.thisisant.com/resources/fit}{FIT SDK}). Hence, this
#'   function requires that Java (JRE/JDK) binaries be on the system path.
#'
#' @export
read_fit <- function(file = file.choose(), format = TRUE, CP = NULL, sRPE = NULL) {
  if (system2("java", stdout = FALSE, stderr = FALSE) == 127)
    stop("java binaries not installed/on system path", call. = FALSE)
  FitCSVTool <- system.file("java/FitCSVTool.jar", package = "cycleRtools")
  if (Sys.info()["sysname"] == "Windows")
    FitCSVTool <- Sys.which(FitCSVTool)
  tmpf <- tempfile()
  message("Reading .fit file...")
  # Records --------------------------------------------------------------------
  system2("java", args = c(
    "-jar", FitCSVTool, "-b", file, tmpf,
    "--data", "record", "--defn", "none"
  ), stdout = FALSE)
  record_data <- read.csv(paste0(tmpf, "_data.csv"))
  # Laps -----------------------------------------------------------------------
  system2("java", args = c(
    "-jar", FitCSVTool, "-b", file, tmpf,
    "--data", "lap", "--defn", "none"
  ), stdout = FALSE)
  lap_data <- read.csv(paste0(tmpf, "_data.csv"))

  record_data$lap <- 0
  i <- match(lap_data$lap.start_time, record_data$record.timestamp.s.)
  record_data$lap[i] <- 1
  record_data$lap[1] <- 1
  record_data$lap    <- factor(cumsum(record_data$lap), ordered = TRUE)
  # Format----------------------------------------------------------------------
  if (format) {
    record_data <- format_fit(record_data)
    record_data <- assign_attrs(record_data, CP, sRPE)
  }
  record_data
}
#' @describeIn read_ride Read a Training Peaks .pwx file. Requires the "xml2" package
#'   to be installed.
#' @export
read_pwx <- function(file = file.choose(), format = TRUE, CP = NULL, sRPE = NULL) {
  x    <- read_xml(file)
  ns   <- xml_ns(x)
  cols <- suppressWarnings(
    xml_name(xml_children(xml_find_one(x, "//d1:sample", ns = ns)))
  )
  message("Reading .pwx file...") # --------------------------------------------

  data <- lapply(cols, function(c) {
    as.numeric(xml_text(xml_find_all(x, paste0("//d1:sample//d1:", c), ns)))
  })
  names(data) <- cols

  # Deal with missing fields ---------------------------------------------------
  len <- vapply(data, length, numeric(1))
  if (length(unique(len)) > 1) {
    message("Resolving missing data points...")
    issues  <- names(len[len < max(len, na.rm = TRUE)])
    .issues <- cols[vapply(issues, function(x) which(grepl(x, cols)), numeric(1))]
    nds     <- as_list(xml_find_all(x, "//d1:sample", ns))
    logi    <- vapply(.issues, function(i) {
      tr    <- lapply(nds, function(n)
        try(xml_find_one(n, i, ns), silent = TRUE))
      out   <- !vapply(tr, inherits, what = "try-error", logical(1))
      return(out)
    }, logical(length(nds)))
    colnames(logi)   <- issues
    to_correct       <- match(issues, names(data))
    data[to_correct] <- lapply(to_correct, function(i)
      ifelse(logi[, names(data[i])], data[[i]], NA))
  }
  # Timestamp ------------------------------------------------------------------
  ts <- xml_text(xml_find_one(x, "//d1:time", ns = ns))
  ts <- sub("T", " ", ts)
  ts <- strptime(ts, "%Y-%m-%d %H:%M:%S")
  data$timestamp <- ts + data$timeoffset

  data <- as.data.frame(data)
  if (format) {
    data <- format_pwx(data)
    data <- assign_attrs(data, CP, sRPE)
  }
  data
}
#' @describeIn read_ride Read an SRM (.srm) file. This requires
#'   \href{http://www.zuto.de/project/srmio/}{Rainer Clasen's srmio library} to
#'   be installed and on the system path.
#' @export
read_srm <- function(file = file.choose(), format = TRUE, CP = NULL, sRPE = NULL) {
  if (system2("srmcmd", stdout = FALSE, stderr = FALSE) == 127)
    stop(paste(
      "srmio library not installed/on the system path.",
      "visit http://www.zuto.de/project/srmio/", sep = "\n"),
      call. = FALSE)
  message("Reading .srm file...")
  tmpf     <- tempfile()
  system2("srmcmd", args = c("--read", file), stdout = tmpf)
  data <- read.table(tmpf, header = TRUE, sep = "\t")
  if (format) {
    data <- format_srm(data)
    data <- assign_attrs(data, CP, sRPE)
  }
  data
}
#' @describeIn read_ride Read a Garmin .tcx file. Requires the "xml2" package to be
#'   installed.
#' @export
read_tcx <- function(file = file.choose(), format = TRUE, CP = NULL, sRPE = NULL) {
  x    <- read_xml(file)
  ns   <- xml_ns(x)
  cols <- suppressWarnings(
    xml_name(xml_children(xml_find_one(x, "//d1:Trackpoint", ns)))
  )
  cols <- paste0("//d1:", cols)

  # Get extra data columns as appropriate --------------------------------------
  if (any(grepl("Position", cols))) {
    cols <- cols[!grepl("Position", cols)]
    tmp  <- suppressWarnings(
      xml_name(xml_children(xml_find_one(x, "//d1:Position", ns)))
    )
    cols <- c(cols, paste0("//d1:", tmp))
  }
  if (any(grepl("Extensions", cols))) {
    cols <- cols[!grepl("Extensions", cols)]
    tmp  <- suppressWarnings(
      xml_name(xml_children(xml_find_one(x, "//ns3:TPX", ns)))
    )
    cols <- c(cols, paste0("//ns3:", tmp))
  }

  # Only look in trackpoints.
  trcols <- paste0("//d1:Trackpoint", cols)

  message("Reading .tcx file...") # --------------------------------------------

  data <- lapply(trcols, function(c) {
    out <- xml_text(xml_find_all(x, c, ns))
    if (all(!is.na(suppressWarnings(as.numeric(out)))))
      out <- as.numeric(out)  # Don't coerce character columns.
    out
  })
  names(data) <- vapply(strsplit(cols, ":"), function(x) x[length(x)], character(1))

  # Deal with missing fields ---------------------------------------------------
  len <- vapply(data, length, numeric(1))
  if (length(unique(len)) > 1) {
    message("Resolving missing data points...")
    issues  <- names(len[len < max(len, na.rm = TRUE)])
    .issues <- cols[vapply(issues, function(x) which(grepl(x, cols)), numeric(1))]
    nds     <- as_list(xml_find_all(x, "//d1:Trackpoint", ns))
    logi    <- vapply(.issues, function(i) {
      tr    <- lapply(nds, function(n)
        try(xml_find_one(n, paste0(".", i), ns), silent = TRUE))
      out   <- !vapply(tr, inherits, what = "try-error", logical(1))
      return(out)
    }, logical(length(nds)))
    colnames(logi)   <- issues
    to_correct       <- match(issues, names(data))
    data[to_correct] <- lapply(to_correct, function(i)
      ifelse(logi[, names(data[i])], data[[i]], NA))
  }

  # Reformat timestamps --------------------------------------------------------
  data$Time   <- as.POSIXct(
    gsub("[[:upper:]]", " ", data$Time),
    origin = Sys.time() - as.numeric(Sys.time()),
    format = "%Y-%m-%d %H:%M:%S ")

  data <- as.data.frame(data)
  if (format) {
    data  <- format_tcx(data)
    data  <- assign_attrs(data, CP, sRPE)
  }
  data
}


# read_GC <- function(file = file.choose(), format = TRUE, CP = NULL, sRPE = NULL) {
#   data <- jsonlite::fromJSON(readLines(file))
# }

assign_attrs <- function(data, CP, sRPE) {
  if (!is.null(CP) && is.numeric(CP)) {
    if (all(!is.na(data$power.smooth.W)))
      data$Wexp.kJ <- with(data, Wbal_(timer.s, power.smooth.W, CP)) / 1000 # kJ.
    attr(data, "CP") <- CP
  }
  if (!is.null(sRPE) && is.numeric(sRPE))
    attr(data, "sRPE") <- sRPE

  class(data) <- c("cycleRdata", "data.frame")

  data
}
