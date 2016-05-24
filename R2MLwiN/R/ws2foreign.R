#' Converts data from one format to another within MLwiN
#'
#' This function converts data from one format to another within MLwiN, via
#' MLwiN macro language. The \code{foreign} package allows R to read data files
#' for most of these formats.
#'
#' MLwiN supports conversion between MLwiN (*.wsz, *.ws), Minitab (*.mtw), SAS
#' (*.xpt), SPSS (*.sav), and Stata (*.dta) files.
#'
#' @param wsfile A file name specifying the data file (with a specific
#' extension) to be converted.
#' @param foreignfile A file name specifying the data file (with a specific
#' extension) after conversion.
#' @param MLwiNPath A path to the MLwiN folder. By default, \code{MLwiNPath = NULL}
#' and path set by \code{options('MLwiN_path')}, the default for which can be
#' changed via \code{options(MLwiN_path = 'path/to/MLwiN vX.XX/')}).
#' @param x64 A logical value indicating whether the 64 bit version of MLwiN is
#' used, the default for this is determined by the version of R used. If
#' \code{FALSE}, the 32 bit version is called.
#'
#' @details
#' The converted data file (with a specific extension) will be saved to the file
#' specified by \code{foreignfile}.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso \code{\link[foreign]{read.dta}}
#'
#' @examples
#'
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' wsfile = normalizePath(file.path(getOption("MLwiN_path"), "samples/tutorial.ws"), winslash = "/")
#' ## the tutorial.dta file will be saved in the temporary folder
#' inputfile = normalizePath(file.path(tempdir(), "tutorial.dta"), winslash = "/", mustWork = FALSE)
#' ws2foreign(wsfile, foreignfile = inputfile)
#' library(foreign)
#' indata = read.dta(inputfile)
#' str(indata)
#' unlink(inputfile)
#' }
#'
#' @export
ws2foreign <- function(wsfile, foreignfile, MLwiNPath = NULL, x64 = NULL) {
  ## Convert MLwiN worksheet file to other data file which is used in Minitab, SAS, SPSS, or Stata
  temptfile <- normalizePath(tempfile("coversion_", fileext = ".txt"), winslash = "/", mustWork = FALSE)
  cat(file = temptfile)
  write("ECHO     0", temptfile, append = TRUE)
  write(paste("LOAD   '", wsfile, "'", sep = ""), temptfile, append = TRUE)
  write(paste("STOR   '", foreignfile, "'", sep = ""), temptfile, append = TRUE)
  write("EXIT", temptfile, append = TRUE)

  if (is.null(x64)) {
    if (.Machine$sizeof.pointer == 8) {
      x64 <- TRUE
    } else {
      x64 <- FALSE
    }
  }

  if (is.null(MLwiNPath)) {
    MLwiNPath <- getOption("MLwiN_path")
  }

  pathinfo <- file.info(MLwiNPath)
  if (is.na(pathinfo$isdir)) {
    stop(paste0(MLwiNPath, " does not exist"))
  }

  if (!isTRUE(pathinfo$isdir)) {
    if (file.access(MLwiNPath, mode = 1) == 0) {
      cmd <- MLwiNPath
    } else {
      stop(paste0(MLwiNPath, " is not executable"))
    }
  }

  if (isTRUE(pathinfo$isdir)) {
    if (x64) {
      cmd <- paste0(MLwiNPath, "/x64/mlnscript.exe")
      if (file.access(cmd, mode = 1) != 0) {
        cmd <- paste0(MLwiNPath, "/i386/mlnscript.exe")
        if (file.access(cmd, mode = 1) != 0) {
          cmd <- paste0(MLwiNPath, "/mlwin.exe")
          if (file.access(cmd, mode = 1) != 0) {
            stop("Cannot find valid MLwiN executable")
          }
        }
      }
    } else {
      cmd <- paste0(MLwiNPath, "/i386/mlnscript.exe")
      if (file.access(cmd, mode = 1) != 0) {
        cmd <- paste0(MLwiNPath, "/mlwin.exe")
        if (file.access(cmd, mode = 1) != 0) {
          stop("Cannot find valid MLwiN executable")
        }
      }
    }
  }

  args <- paste0("/nogui /run ", "\"", temptfile, "\"")
  system2(cmd, args = args, stdout = stdout, stderr = stderr)
  file.remove(temptfile)
  cat("\n")
}
