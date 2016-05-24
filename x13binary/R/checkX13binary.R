#' Check if X-13ARIMA-SEATS Runs Properly
#' 
#' Performs a test run of X-13ARIMA-SEATS. Fails if no output is produced.
#' 
#' @param fail.unsupported  logical, whether beeing on an unsupported platform
#'   leads to an error.
#' @param verbose  logical, should a message be returned on success?
#' @examples
#' checkX13binary()
#' 
#' @export
checkX13binary <- function(fail.unsupported = FALSE, verbose = TRUE){

  if (supportedPlatform()){
    if (.Platform$OS.type == "windows"){    
      x13.bin <- system.file("bin", "x13ashtml.exe", package="x13binary")
    } else {
      x13.bin <- system.file("bin", "x13ashtml", package="x13binary")
    }
    if (x13.bin == ""){
      stop("X-13 binary file not found")
    }

    tdir <- tempdir()
    file.copy(system.file("testdata", "Testairline.spc", package="x13binary"), tdir)
    if (.Platform$OS.type == "windows") {

      # change wd on win as X-13 writes `fort.6` to it
      owd <- getwd()
      on.exit(setwd(owd))
      setwd(tdir)

      # shell() gives a more verbose output on windows
      sout <- shell(paste0("\"", normalizePath(x13.bin), "\"", " Testairline"), 
                    intern = TRUE)

      if (isTRUE(attr(sout,"status") != 0)){
        if (verbose){
          packageStartupMessage("Rerunning with full console output:")
          shell(paste0("\"", normalizePath(x13.bin), "\"", " Testairline"))
        }
        stop("When running\n\n  ", x13.bin, 
             "\n\nCommand Prompt returned the following message:\n\n", 
             paste(strwrap(sout, indent = 2, exdent = 2), collapse = "\n"),
             "\n\n")
      }
    } else {
      sout <- system(paste(x13.bin, file.path(tdir, "Testairline")), intern = TRUE)
      if (isTRUE(attr(sout,"status") != 0)){
        stop("When running\n\n  ", x13.bin, 
             "\n\nthe system returned the following message:\n\n", 
             sout,
             "\n\n")
      }
      # drop error if output contains the word ERROR
      # (This does not necessarily lead to a non zero exit status)
      if (inherits(sout, "character")){
        if (any(grepl("ERROR", sout))){
          stop("When running\n\n  ", x13.bin, 
               "\n\nthe system returned the following message:\n\n", 
               sout,
               "\n\n")
        }
      }
    }
    if (!file.exists(file.path(tdir, "Testairline.html"))){
      if (inherits(sout, "character")){
        stop("X-13 has not produced Testairline.html.\n", 
             "When running\n\n  ", x13.bin, 
             "\n\nthe system returned the following message:\n\n", 
             sout,
             "\n\n")
      } else {
        stop("X-13 has not produced Testairline.html.\n", 
             "When running\n\n  ", x13.bin, 
             "\n\nthe system has returned no message.\n\n")
      }
    }
    if (verbose){
      packageStartupMessage("x13binary is working properly")
    }

  } else {
    ifelse(fail.unsupported, stop, packageStartupMessage)(
      "Unsupported platform: ", Sys.info()["sysname"], Sys.info()["release"],
      "\nFor this platform, there are currently no binaries of X-13ARIMA-SEATS.")
    return(invisible(FALSE))
  }

  invisible(TRUE)
}



