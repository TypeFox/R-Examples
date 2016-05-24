##' R interface to the PHREEQC geochemical modeling program.
##' 
##' An interface to PHREEQC (Version 3). PHREEQC is a geochemical
##' modeling program developed by the US Geological Survey that is designed
##' to perform a wide variety of aqueous geochemical calculations, 
##' including speciation, batch-reaction, one-dimensional 
##' reactive-transport, and inverse geochemical calculations.    
##' 
##' \tabular{ll}{Package: \tab phreeqc\cr Type: \tab Package\cr Version: \tab 3.3.1\cr Date: \tab 2015-09-30\cr License: \tab BSD_3_clause + file LICENSE\cr}
##' 
##' @name phreeqc-package
##' @aliases phreeqc-package phreeqc
##' @docType package
##' @author S.R. Charlton, D.L. Parkhurst, and C.A.J. Appelo, with contributions from D. Gillespie for Chipmunk BASIC (p2c) and S.D. Cohen, A.C. Hindmarsh, R. Serban, D. Shumaker, and A.G. Taylor for CVODE (SUNDIALS) \cr Maintainer: S.R. Charlton \email{charlton@@usgs.gov}
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}\cr \url{http://computation.llnl.gov/casc/sundials/main.html}
##' @keywords package
##' @examples
##'
##' #########################################################################
##' # Run ex2 and plot results
##' #########################################################################
##' 
##' # load the phreeqc.dat database
##' phrLoadDatabaseString(phreeqc.dat)
##' 
##' # run example 2
##' phrRunString(ex2)
##' 
##' # retrieve selected_output as a list of data.frame
##' so <- phrGetSelectedOutput()
##' 
##' # plot the results
##' attach(so$n1)
##' title  <- "Gypsum-Anhydrite Stability"
##' xlabel <- "Temperature, in degrees celcius"
##' ylabel <- "Saturation index"
##' plot(temp.C., si_gypsum, main = title, xlab = xlabel, ylab = ylabel,
##'      col = "darkred", xlim = c(25, 75), ylim = c(-0.4, 0.0))
##' points(temp.C., si_anhydrite, col = "darkgreen")
##' legend("bottomright", c("Gypsum", "Anhydrite"),
##'        col = c("darkred", "darkgreen"), pch = c(1, 1))
NULL



# Package Functions



##' Accumulate line(s) for input to phreeqc.
##' 
##' Appends a line of text to the input buffer in order to be run using
##' \code{\link{phrRunAccumulated}}.
##' 
##' @export phrAccumulateLine
##' @useDynLib phreeqc
##' @param line the line(s) to add for input to phreeqc.
##' @return NULL
##' @family Accumulate
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # this example loads the phreeqc.dat database, accumulates input, and
##' # runs it
##' phrLoadDatabaseString(phreeqc.dat)
##' phrAccumulateLine("TITLE Example 2.--Temperature dependence of solubility")
##' phrAccumulateLine("                  of gypsum and anhydrite")
##' phrAccumulateLine("SOLUTION 1 Pure water")
##' phrAccumulateLine("        pH      7.0")
##' phrAccumulateLine("        temp    25.0")
##' phrAccumulateLine("EQUILIBRIUM_PHASES 1")
##' phrAccumulateLine("        Gypsum          0.0     1.0")
##' phrAccumulateLine("        Anhydrite       0.0     1.0")
##' phrAccumulateLine("REACTION_TEMPERATURE 1")
##' phrAccumulateLine("        25.0 75.0 in 51 steps")
##' phrAccumulateLine("SELECTED_OUTPUT")
##' phrAccumulateLine("        -file   ex2.sel")
##' phrAccumulateLine("        -temperature")
##' phrAccumulateLine("        -si     anhydrite  gypsum")
##' phrAccumulateLine("END")
##' phrSetOutputFileOn(TRUE)
##' if (is.null(phrRunAccumulated())) {
##'   cat(paste("see ", phrGetOutputFileName(), ".\n", sep = ""))
##' }
##' 
phrAccumulateLine <-
function(line) {
  invisible(.Call("accumLineLst", as.character(line), PACKAGE = .packageName))
}



##' Clear the accumulated input buffer.
##' 
##' Clears the accumulated input buffer. The input buffer is accumulated from
##' calls to the \code{\link{phrAccumulateLine}} method.
##' 
##' @export phrClearAccumulatedLines
##' @useDynLib phreeqc
##' @return NULL
##' @family Accumulate
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example loads some keyword input, clears the input, and displays
##' # the results.
##' phrAccumulateLine("SOLUTION 1")
##' phrAccumulateLine("END")
##' cat("The accumulated input is:", phrGetAccumulatedLines(), sep = "\n")
##' phrClearAccumulatedLines()
##' cat("The accumulated input now is:\n", phrGetAccumulatedLines(), sep = "\n")
##' 
phrClearAccumulatedLines <-
function() {
  invisible(.Call("clearAccum", PACKAGE = .packageName))
}



##' Retrieve the accumulated input.
##' 
##' Returns the accumulated input as a character vector.
##' 
##' @export phrGetAccumulatedLines
##' @useDynLib phreeqc
##' @return A character vector containing the accumulated input.
##' @family Accumulate
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example loads some keyword input and displays the contents.
##' phrAccumulateLine("SOLUTION 1")
##' phrAccumulateLine("END")
##' cat("The accumulated input is:", phrGetAccumulatedLines(), sep = "\n")
##' 
phrGetAccumulatedLines <-
function() {
  return(.Call("getAccumLines", PACKAGE = .packageName))
}



##' Retrieve a list containing the current list of components.
##' 
##' @export phrGetComponentList
##' @useDynLib phreeqc
##' @return A list containing the names of the components defined in the current system.
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example runs the ex2 input file and echos the list of components.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrRunString(ex2)
##' cat("components:\n")
##' for (comp in phrGetComponentList()) {
##'   cat(comp, "\n")
##' }
##' 
phrGetComponentList <-
function() {
  return(.Call("listComps", PACKAGE = .packageName))
}



##' Retrieve the name of the dump file.
##' 
##' Retrieves the name of the dump file. This file name is used if not
##' specified within DUMP input. The default value is dump.0.out.
##' 
##' @export phrGetDumpFileName
##' @useDynLib phreeqc
##' @return The name of the dump file as a string.
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetDumpFileOn(TRUE)
##' 
##' input <-              "SOLUTION 1 Pure water     \n"
##' input <- paste(input, "EQUILIBRIUM_PHASES 1      \n")
##' input <- paste(input, "    Calcite 0 10          \n")
##' input <- paste(input, "SAVE solution 1           \n")
##' input <- paste(input, "SAVE equilibrium_phases 1 \n")
##' input <- paste(input, "DUMP                      \n")
##' input <- paste(input, "    -solution 1           \n")
##' input <- paste(input, "    -equilibrium_phases  1\n")
##' 
##' if (!is.null(phrRunString(input))) {
##'   cat(phrGetErrorStrings())
##' }
##' cat(paste("see ", phrGetDumpFileName(), "."))
##' 
##' 
phrGetDumpFileName <-
function() {
  return(.Call("getDumpFileName", PACKAGE = .packageName))
}



##' Retrieve DUMP strings.
##'
##' Retrieves the string buffer containing DUMP output as a character
##' vector.
##' 
##' @export phrGetDumpStrings
##' @useDynLib phreeqc
##' @return The dump output as a character vector.
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetDumpStringsOn(TRUE)
##' 
##' input <-              "SOLUTION 1 Pure water     \n"
##' input <- paste(input, "EQUILIBRIUM_PHASES 1      \n")
##' input <- paste(input, "    Calcite 0 10          \n")
##' input <- paste(input, "SAVE solution 1           \n")
##' input <- paste(input, "SAVE equilibrium_phases 1 \n")
##' input <- paste(input, "DUMP                      \n")
##' input <- paste(input, "    -solution 1           \n")
##' input <- paste(input, "    -equilibrium_phases 1 \n")
##' 
##' if (!is.null(phrRunString(input))) {
##'   cat(phrGetErrorStrings(), sep = "\n")
##' }
##' cat(phrGetDumpStrings(), sep = "\n")
##'
phrGetDumpStrings <-
function() {
  return(.Call("getDumpStrings", PACKAGE = .packageName))
}



##' Retrieve the name of the error file.
##' 
##' Retrieves the name of the error file. The default value is phreeqc.0.err.
##'
##' The error file switch must be set using the \code{\link{phrSetErrorFileOn}} function.
##' 
##' @export phrGetErrorFileName
##' @useDynLib phreeqc
##' @return The name of the error file as a string.
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetErrorFileName <-
function() {
  return(.Call("getErrorFileName", PACKAGE = .packageName))
}



##' Retrieve the current value of the dump file switch.
##' 
##' @export phrGetDumpFileOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetDumpFileOn <-
function() {
  return(.Call("getDumpFileOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the dump strings switch.
##' 
##' @export phrGetDumpStringsOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetDumpStringsOn <-
function() {
  return(.Call("getDumpStringOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the error file switch.
##' 
##' @export phrGetErrorFileOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetErrorFileOn <-
function() {
  return(.Call("getErrorFileOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the error strings switch.
##' 
##' @export phrGetErrorStringsOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetErrorStringsOn <-
function() {
  return(.Call("getErrorStringOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the log file switch.
##' 
##' @export phrGetLogFileOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Log
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetLogFileOn <-
function() {
  return(.Call("getLogFileOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the log strings switch.
##' 
##' @export phrGetLogStringsOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Log
##' @references \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' 
phrGetLogStringsOn <-
function() {
  return(.Call("getLogStringOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the output file switch.
##' 
##' @export phrGetOutputFileOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetOutputFileOn <-
function() {
  return(.Call("getOutputFileOn", PACKAGE = .packageName))
}



##' Retrieve the current value of the output strings switch.
##' 
##' @export phrGetOutputStringsOn
##' @useDynLib phreeqc
##' @return TRUE if errors are currently being written to file.
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' 
phrGetOutputStringsOn <-
function() {
  return(.Call("getOutputStringOn", PACKAGE = .packageName))
}



##' Retrieve error string messages.
##' 
##' Retrieves a character vector containing any error messages that were generated
##' during the last invocation of the following methods:
##' \code{\link{phrAccumulateLine}}, \code{\link{phrLoadDatabase}},
##' \code{\link{phrLoadDatabaseString}}, \code{\link{phrRunAccumulated}},
##' \code{\link{phrRunFile}}, \code{\link{phrRunString}}
##' 
##' This routine is rarely needed when running interactively since the error
##' string is displayed when it occurs.
##' 
##' @export phrGetErrorStrings
##' @useDynLib phreeqc
##' @return The error messages as a character vector.
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # loaddatabase should fail
##' n <- try(phrLoadDatabase("missing.dat"), silent = TRUE)
##' # if n is non-NULL display error string
##' if (!is.null(n)) phrGetErrorStrings()
##' 
phrGetErrorStrings <-
function() {
  return(.Call("getErrorStrings", PACKAGE = .packageName))
}



##' Retrieve the name of the log file.
##' 
##' Retrieves the name of the log file. The default name is phreeqc.0.log.
##' 
##' @export phrGetLogFileName
##' @useDynLib phreeqc
##' @return The name of the log file as a string.
##' @family Log
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # This example checks to see if the log file is turned on
##' # and prints the appropriate message
##' if (phrGetLogFileOn()) {
##'   cat("The name of the log file (is/will be):", phrGetLogFileName(), "\n")
##' } else {
##'   cat("The log file is not turned on\n")
##' }
##' 
phrGetLogFileName <-
function() {
  return(.Call("getLogFileName", PACKAGE = .packageName))
}



##' Retrieve log output.
##' 
##' Retrieves the string buffer containing phreeqc log output.
##' 
##' @export phrGetLogStrings
##' @useDynLib phreeqc
##' @return A character vector containing phreeqc log output.
##' @family Log
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with gypsum with the output file on.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputFileOn(TRUE)
##'
##' input <- vector(mode="character")
##' input <- c(input, "SOLUTION 1 Pure water ")
##' input <- c(input, "EQUILIBRIUM_PHASES 1  ")
##' input <- c(input, "  Gypsum 0 10         ")
##' input <- c(input, "KNOBS                 ")
##' input <- c(input, "  -logfile TRUE       ")
##'
##' if (is.null(phrRunString(input))) {
##'   log <- phrGetLogStrings()
##' }
##' 
phrGetLogStrings <-
function() {
  return(.Call("getLogStrings", PACKAGE = .packageName))
}



##' Retrieve the name of the output file.
##' 
##' Retrieves the name of the output file. The default name is phreeqc.0.out.
##' 
##' @export
##' @useDynLib phreeqc
##' @return The name of the output file as a string.
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with gypsum with the output file on.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputFileOn(TRUE)
##'
##' input <- vector(mode="character")
##' input <- c(input, "SOLUTION 1 Pure water ")
##' input <- c(input, "EQUILIBRIUM_PHASES 1  ")
##' input <- c(input, "  Gypsum 0 10         ")
##'
##' if (is.null(phrRunString(input))) {
##'   output <- readLines(phrGetOutputFileName())
##'   unlink(phrGetOutputFileName())  # tidy up
##' }
##' 
phrGetOutputFileName <-
function() {
  return(.Call("getOutputFileName", PACKAGE = .packageName))
}



##' Retrieve standard phreeqc output.
##' 
##' Retrieves the phreeqc output as a character vector.
##' 
##' A NULL value is returned when there is no selected-output.
##' 
##' @export phrGetOutputStrings
##' @useDynLib phreeqc
##' @return A character vector containing phreeqc output.
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite and displays
##' # the results
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##'
##' input <- vector(mode="character")
##' input <- c(input, "SOLUTION 1 Pure water ")
##' input <- c(input, "EQUILIBRIUM_PHASES 1  ")
##' input <- c(input, "  Gypsum 0 10         ")
##'
##' if (is.null(phrRunString(input))) {
##'   cat(phrGetOutputStrings(), sep = "\n")
##' }
##' 
phrGetOutputStrings <-
function() {
  return(.Call("getOutputStrings", PACKAGE = .packageName))
}



##' Returns the contents of the selected output as a list of data frames.
##' 
##' phrGetSelectedOutput returns a named list containing the resultant
##' selected output blocks.  The names of each data frame are creating by
##' concatenating the letter 'n' and the user number of the selected output
##' block.
##'
##' phrGetSelectedOutput uses the \code{\link{make.names}} function to create
##' valid column names. The allow_ variable is passed to
##' \code{\link{make.names}} and is used for backward compatibility.
##' 
##' @export phrGetSelectedOutput
##' @useDynLib phreeqc
##' @param allow_ used for compatibility with R prior to 1.9.0 (default is TRUE).
##' @return Returns a named list of data frames containing the selected_output from the previous run.
##' @family Selected Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # Load database and run ex2
##' phrLoadDatabaseString(phreeqc.dat)
##' phrRunString(ex2)
##'
##' # display a summary of the results
##' df <- phrGetSelectedOutput()
##' summary(df$n1)
##' 
phrGetSelectedOutput <-
function(allow_ = TRUE) {
  sel_outs <- .Call("getSelOutLst", PACKAGE = .packageName)
  if (!is.null(sel_outs)) {
    for (t in names(sel_outs)) {
      if (!is.null(sel_outs[[t]])) {
        names(sel_outs[[t]]) <- make.names(names(sel_outs[[t]]), unique = TRUE,
                                           allow_ = allow_)
      }
    }
  }
  return(sel_outs)
}



##' Retrieve warning messages.
##' 
##' Returns a character vector containing any warning messages that were
##' generated during the last invocation of the following methods:
##' \code{\link{phrAccumulateLine}}, \code{\link{phrLoadDatabase}},
##' \code{\link{phrLoadDatabaseString}}, \code{\link{phrRunAccumulated}},
##' \code{\link{phrRunFile}}, \code{\link{phrRunString}}.
##' 
##' A NULL value is returned if there are no warnings.
##' 
##' @export phrGetWarningStrings
##' @useDynLib phreeqc
##' @return A character vector containing warning messages or NULL.
##' @family Warning
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example loads the phreeqc.dat database and attempts to use the
##' # DATABASE keyword to set the database to wateq4f.dat.  A warning is
##' # displayed stating that the DATABASE keyword is ignored in the 'R'
##' # implementation.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrAccumulateLine("DATABASE wateq4f.dat")
##' phrAccumulateLine("SOLUTION 1")
##' phrRunAccumulated()
##' if (!is.null(phrGetWarningStrings())) {
##'   cat(phrGetWarningStrings(), sep = "\n")
##' }
##' 
phrGetWarningStrings <-
function() {
  return(.Call("getWarningStrings", PACKAGE = .packageName))
}



##' Load a phreeqc database file
##' 
##' Loads the given phreeqc database into phreeqc.  Returns NULL if successful.
##' 
##' 
##' @export phrLoadDatabase
##' @useDynLib phreeqc
##' @param filename The name of the database file.
##' @return This function returns NULL.
##' @family Load Database
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # create temporary database file
##' tf <- tempfile()
##' writeLines(phreeqc.dat, tf)
##' 
##' if (is.null(phrLoadDatabase(tf))) {
##'   cat("database ok\n")
##' } else {
##'   cat("database contains errors\n")
##' }
##'
##' # delete temporary database file
##' unlink(tf)
##' 
phrLoadDatabase <-
function(filename) {
  invisible(.Call("loadDB", as.character(filename), PACKAGE = .packageName))
}



##' Load a phreeqc database as a string or a list of strings.
##' 
##' Load the specified string(s) as a database into phreeqc. Returns NULL if
##' successful.
##' 
##' All previous definitions are cleared.
##' 
##' @export phrLoadDatabaseString
##' @useDynLib phreeqc
##' @param input String containing data to be used as the phreeqc database.
##' @return This function returns NULL.
##' @family Load Database
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @keywords interface
##' @examples
##' 
##' # this example loads the phreeqc.dat database, turns on the
##' # output file and runs ex2 as a string
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputFileOn(TRUE)
##' if (is.null(phrRunString(ex2))) {
##'   cat(paste("see ", phrGetOutputFileName(), ".\n", sep = ""))
##' }
##' 
phrLoadDatabaseString <-
function(input) {
  invisible(.Call("loadDBLst", as.character(input), PACKAGE = .packageName))
}



##' Runs the accumulated input.
##' 
##' Runs the input buffer as defined by calls to \code{\link{phrAccumulateLine}}.
##' 
##' After calling \code{phrRunAccumulated} \code{\link{phrGetAccumulatedLines}} can
##' be used in case of errors. The accumulated lines will be cleared on the next call
##' to \code{\link{phrAccumulateLine}}.
##'
##' The \code{phrAccumulateLine} method cannot be called until a database
##' has been successfully loaded by calls to either
##' \code{\link{phrLoadDatabase}} or \code{\link{phrLoadDatabaseString}}.
##' 
##' @export phrRunAccumulated
##' @useDynLib phreeqc
##' @return This function returns NULL on success.
##' @family Accumulate
##' @family Run
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # turn on the output file
##' phrSetOutputFileOn(TRUE)
##' 
##' # load the phreeqc.dat database
##' phrLoadDatabaseString(phreeqc.dat)
##' 
##' # accumulate the input
##' phrAccumulateLine("TITLE Example 2.--Temperature dependence of solubility")
##' phrAccumulateLine("                  of gypsum and anhydrite")
##' phrAccumulateLine("SOLUTION 1 Pure water")
##' phrAccumulateLine("        pH      7.0")
##' phrAccumulateLine("        temp    25.0")
##' phrAccumulateLine("EQUILIBRIUM_PHASES 1")
##' phrAccumulateLine("        Gypsum          0.0     1.0")
##' phrAccumulateLine("        Anhydrite       0.0     1.0")
##' phrAccumulateLine("REACTION_TEMPERATURE 1")
##' phrAccumulateLine("        25.0 75.0 in 51 steps")
##' phrAccumulateLine("SELECTED_OUTPUT")
##' phrAccumulateLine("        -file   ex2.sel")
##' phrAccumulateLine("        -temperature")
##' phrAccumulateLine("        -si     anhydrite  gypsum")
##' phrAccumulateLine("END")
##' 
##' # run it and echo the name of the output file
##' if (is.null(phrRunAccumulated())) {
##'   cat(paste("see ", phrGetOutputFileName(), ".\n", sep = ""))
##' }
phrRunAccumulated <-
function() {
  invisible(.Call("runAccum", PACKAGE = .packageName))
}



##' Run phreeqc input file
##' 
##' phrRunFile executes a phreeqc run using a file as input
##' 
##' @export phrRunFile
##' @useDynLib phreeqc
##' @param filename The file name of the phreeqc input file.
##' @return This function returns NULL on success.
##' @family Run
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # load the phreeqc.dat database
##' phrLoadDatabaseString(phreeqc.dat)
##'
##' # create ex2 if it doesn't exist
##' if (!file.exists("ex2")) writeLines(ex2, "ex2")
##' 
##' # run ex2
##' if (is.null(phrRunFile("ex2"))) {
##'   cat("use phrGetSelectedOutput() to see results.")
##' }
##'
##' unlink("ex2")  # tidy up
##' 
phrRunFile <-
function(filename) {
  invisible(.Call("runFile", as.character(filename), PACKAGE = .packageName))
}



##' Runs phreeqc using the given string as input.
##' 
##' Runs phreeqc using the given string as input. Returns the number of
##' errors encountered during the run.
##' 
##' The \code{RunString} method cannot be called until a database has
##' been successfully loaded by one of the following the LoadDatabase
##' methods \code{\link{phrLoadDatabase}}, \code{\link{phrLoadDatabaseString}}.
##' 
##' @export phrRunString
##' @useDynLib phreeqc
##' @param input character vector containing phreeqc input
##' @return This function returns NULL on success.
##' @family Run
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @keywords interface
##' @examples
##'
##' #
##' # This example accumulates phreeqc input into a character vector
##' # and runs it.
##' #
##'
##' # load phreeqc.dat file
##' phrLoadDatabaseString(phreeqc.dat)
##'
##' # create input
##' input <- vector()
##' input <- c(input, "SOLUTION 1") 
##' input <- c(input, "  temp 25.0") 
##' input <- c(input, "  pH    7.0")
##'
##' # turn on output
##' phrSetOutputFileOn(TRUE)
##'
##' # run input
##' phrRunString(input)
##' cat(paste("see", phrGetOutputFileName(), "."))
##' 
phrRunString <-
function(input) {
  invisible(.Call("runStringLst", as.character(input), PACKAGE = .packageName))
}



##' Set the name of the dump file.
##' 
##' Sets the name of the dump file. This file name is used if not specified
##' within the DUMP keyword block. The default value is dump.0.out.
##' 
##' @export phrSetDumpFileName
##' @useDynLib phreeqc
##' @param filename the name of the file.
##' @return NULL
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # This example equilibrates pure water with calcite and writes the 
##' # dump results to file.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetDumpFileOn(TRUE)
##' phrSetDumpFileName("phreeqc.dump")
##' input <- c( 
##'   'SOLUTION 1 Pure water     ',
##'   'EQUILIBRIUM_PHASES 1      ',
##'   '    Calcite 0 10          ',
##'   'SAVE solution 1           ',
##'   'SAVE equilibrium_phases 1 ',
##'   'DUMP                      ',
##'   '    -solution 1           ',
##'   '    -equilibrium_phases 1 '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetDumpFileName(), "\n")
##' }
##' 
phrSetDumpFileName <-
function(filename) {
  invisible(.Call("setDumpFileName", as.character(filename), PACKAGE = .packageName))
}



##' Set the dump file on/off.
##' 
##' Sets the dump file switch on or off. This switch controls whether or
##' not phreeqc writes to the dump file. The initial setting is off.
##' 
##' @export phrSetDumpFileOn
##' @useDynLib phreeqc
##' @param value if TRUE, captures output normally sent to the dump file into a buffer.
##' @return NULL
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite and writes the 
##' # dump results to file.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetDumpFileOn(TRUE)
##' input <- c( 
##'   'SOLUTION 1 Pure water     ',
##'   'EQUILIBRIUM_PHASES 1      ',
##'   '    Calcite 0 10          ',
##'   'SAVE solution 1           ',
##'   'SAVE equilibrium_phases 1 ',
##'   'DUMP                      ',
##'   '    -solution 1           ',
##'   '    -equilibrium_phases 1 '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetDumpFileName(), "\n")
##' }
##' 
phrSetDumpFileOn <-
function(value) {
  invisible(.Call("setDumpFileOn", as.logical(value), PACKAGE = .packageName))
}



##' Set dump strings on/off.
##' 
##' Sets the dump strings switch on or off. This switch controls whether or
##' not the data normally sent to the dump file are stored in a buffer for
##' retrieval. The initial setting is off.
##'
##' @export phrSetDumpStringsOn
##' @useDynLib phreeqc
##' @param value if TRUE, captures output normally sent to the error file into a buffer.
##' @return NULL
##' @family Dump
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # This example equilibrates pure water with calcite and echos the 
##' # dump strings.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetDumpStringsOn(TRUE)
##' input <- c( 
##'   'SOLUTION 1 Pure water     ',
##'   'EQUILIBRIUM_PHASES 1      ',
##'   '    Calcite 0 10          ',
##'   'SAVE solution 1           ',
##'   'SAVE equilibrium_phases 1 ',
##'   'DUMP                      ',
##'   '    -solution 1           ',
##'   '    -equilibrium_phases 1 '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("Dump:\n")
##'   cat(phrGetDumpStrings(), sep = "\n")
##' }
##' 
phrSetDumpStringsOn <-
function(value) {
  invisible(.Call("setDumpStringOn", as.logical(value), PACKAGE = .packageName))
}



##' Set the name of the error file.
##' 
##' Sets the name of the error file. The default value is phreeqc.0.err.
##' 
##' @export phrSetErrorFileName
##' @useDynLib phreeqc
##' @param filename the name of the file.
##' @return NULL
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # This example equilibrates pure water with calcite and displays
##' # the log file name.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetLogFileOn(TRUE)
##' phrSetLogFileName("phreeqc.log")
##' input <- c( 
##'   'SOLUTION 1 Pure water ',
##'   'EQUILIBRIUM_PHASES 1  ',
##'   '    Calcite 0 10      '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetErrorFileName(), "\n")
##' }
##' 
phrSetErrorFileName <-
function(filename) {
  invisible(.Call("setErrorFileName", as.character(filename), PACKAGE = .packageName))
}



##' Set error file on/off.
##' 
##' Sets the error file switch on or off. This switch controls whether
##' or not phreeqc writes to the error file. The initial setting is off.
##' 
##' The try is necessary to keep the error message from displaying immediately.
##'
##' @export phrSetErrorFileOn
##' @useDynLib phreeqc
##' @param value if TRUE, writes output to the the error file.
##' @return NULL
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # This example attempts to run ex1, fails, and writes the error
##' # message to the error file (no database is loaded).
##' phrSetErrorFileOn(TRUE)
##' phrSetErrorFileName("phreeqc.errors")
##' if (!is.null(try(phrRunString(ex1), silent=TRUE))) {
##'   cat("see", phrGetErrorFileName(), "\n")
##' }
##' 
phrSetErrorFileOn <-
function(value) {
  invisible(.Call("setErrorFileOn", as.logical(value), PACKAGE = .packageName))
}



##' Set error strings on/off.
##' 
##' Sets the error strings switch on or off.  This switch controls whether or
##' not the data normally sent to the error file are stored in a buffer for
##' retrieval. The initial setting is on.
##' 
##' The try is necessary to keep the error message from displaying immediately.
##' 
##' @export phrSetErrorStringsOn
##' @useDynLib phreeqc
##' @param value if TRUE, captures output normally sent to the error file into a buffer.
##' @return NULL
##' @family Error
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##'
##' # This example attempts to run ex1, fails, and displays the error message
##' # (no database is loaded).
##' phrSetErrorStringsOn(TRUE)
##' if (!is.null(try(phrRunString(ex1), silent=TRUE))) {
##'   cat(phrGetErrorStrings(), "\n")
##' }
##' 
phrSetErrorStringsOn <-
function(value) {
  invisible(.Call("setErrorStringOn", as.logical(value), PACKAGE = .packageName))
}



##' Set the name of the log file.
##' 
##' Sets the name of the log file. The default value is phreeqc.0.log
##' 
##' Logging must be enabled through the use of the KNOBS -logfile
##' option in order to receive any log messages.
##' 
##' @export phrSetLogFileName
##' @useDynLib phreeqc
##' @param filename the name of the file.
##' @return NULL
##' @family Log
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite and displays
##' # the log file name.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetLogFileOn(TRUE)
##' phrSetLogFileName("phreeqc.log")
##' input <- c( 
##'   'SOLUTION 1 Pure water ',
##'   'EQUILIBRIUM_PHASES 1  ',
##'   '    Calcite 0 10      ',
##'   'KNOBS                 ',
##'   '    -logfile true     '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetLogFileName(), "\n")
##' }
##' 
phrSetLogFileName <-
function(filename) {
  invisible(.Call("setLogFileName", as.character(filename), PACKAGE = .packageName))
}



##' Set log file on/off.
##' 
##' Sets the log file switch on or off. This switch controls whether
##' or not phreeqc writes log messages to the log file. The initial
##' setting is off.
##'
##' Logging must be enabled through the use of the KNOBS -logfile
##' option in order to receive an log messages.
##' 
##' @export phrSetLogFileOn
##' @useDynLib phreeqc
##' @param value if TRUE, writes output to the the log file.
##' @return NULL
##' @family Log
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example runs ex2 with the log file turned on.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetLogStringsOn(TRUE)
##'
##' # turn logging on
##' phrAccumulateLine("KNOBS; -logfile true")
##' phrRunAccumulated()
##' 
##' if (is.null(phrRunString(ex2))) {
##'   cat("see", phrGetLogFileName(), "\n")
##' }
##' 
phrSetLogFileOn <-
function(value) {
  invisible(.Call("setLogFileOn", as.logical(value), PACKAGE = .packageName))
}



##' Set log strings on/off.
##' 
##' Sets the log strings switch on or off.  This switch controls whether or
##' not the data normally sent to the log file are stored in a buffer for
##' retrieval. The initial setting is off.
##' 
##' @export phrSetLogStringsOn
##' @useDynLib phreeqc
##' @param value if TRUE, captures output normally sent to the log file into a buffer.
##' @return NULL
##' @family Log
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example runs ex2 with log strings turned on.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetLogStringsOn(TRUE)
##'
##' # turn logging on
##' phrAccumulateLine("KNOBS; -logfile true")
##' phrRunAccumulated()
##' 
##' if (is.null(phrRunString(ex2))) {
##'   cat(phrGetLogStrings(), sep = "\n")
##' }
##' 
phrSetLogStringsOn <-
function(value) {
  invisible(.Call("setLogStringOn", as.logical(value), PACKAGE = .packageName))
}



##' Set the name of the output file.
##' 
##' Sets the name of the output file. The default value is phreeqc.0.out.
##'
##' The output file must be turned on using the \code{\link{phrSetOutputFileOn}} function.
##' 
##' @export phrSetOutputFileName
##' @useDynLib phreeqc
##' @param filename the name of the file.
##' @return NULL
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite and displays
##' # the resulting file name.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputFileOn(TRUE)
##' phrSetOutputFileName("phreeqc.output")
##' input <- c( 
##'   'SOLUTION 1 Pure water ',
##'   'EQUILIBRIUM_PHASES 1  ',
##'   '    Calcite 0 10      '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetOutputFileName(), "\n")
##' }
##' 
phrSetOutputFileName <-
function(filename) {
  invisible(.Call("setOutputFileName", as.character(filename), PACKAGE = .packageName))
}



##' Set output file on/off.
##' 
##' Sets the output file switch on or off. This switch controls whether
##' or not phreeqc writes to the output file. This is the output normally
##' generated when phreeqc is run. The initial setting is off.
##' 
##' @export phrSetOutputFileOn
##' @useDynLib phreeqc
##' @param value if TRUE, writes output to the the output file.
##' @return NULL
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example runs ex2 with the output file turned on.
##'
##' # write temporary input file
##' tf <- tempfile()
##' writeLines(ex2, tf)
##'
##' # load database and run input file
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputFileOn(TRUE)
##' if (is.null(phrRunFile(tf))) {
##'   cat("see", phrGetOutputFileName(), "\n")
##' }
##'
##' # delete temporary input file
##' unlink(tf)
##' 
phrSetOutputFileOn <-
function(value) {
  invisible(.Call("setOutputFileOn", as.logical(value), PACKAGE = .packageName))
}



##' Set output strings on/off.
##' 
##' Set the output string switch on or off.  This switch controls whether or
##' not the data normally sent to the output file are stored in a buffer for
##' retrieval. The initial setting is off.
##' 
##' The output strings setting is only used by the Run methods:
##' \code{\link{phrRunAccumulated}}, \code{\link{phrRunFile}},
##' \code{\link{phrRunString}}.
##'
##' @export phrSetOutputStringsOn
##' @useDynLib phreeqc
##' @param value if TRUE, captures output normally sent to the output file into a buffer.
##' @return NULL
##' @family Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite and displays
##' # the results.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' input <- c( 
##'   'SOLUTION 1 Pure water ',
##'   'EQUILIBRIUM_PHASES 1  ',
##'   '    Calcite 0 10      '
##'   )
##' 
##' if (is.null(phrRunString(input))) {
##'   cat(phrGetOutputStrings(), sep = "\n")
##' }
##' 
phrSetOutputStringsOn <-
function(value) {
  invisible(.Call("setOutputStringOn", as.logical(value), PACKAGE = .packageName))
}



##' Retrieve the name of the selected_output file.
##' 
##' Retrieves the name of the selected_output file. The default value is selected_{nuser}.0.out.
##'
##' The selected_output file must be turned on using the \code{\link{phrSetSelectedOutputFileOn}} function.
##' 
##' @export phrGetSelectedOutputFileName
##' @useDynLib phreeqc
##' @param nuser the user number specified in the SELECTED_OUTPUT block.
##' @return The name of the selected_output file as a string.
##' @family Selected Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite at various temperatures
##' # and displays the name of the selected_output file.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetSelectedOutputFileOn(1, TRUE)
##' phrSetSelectedOutputFileName(1, "ex2.sel")
##'
##' input <- c( 
##'   'SOLUTION 1 Pure water     ',
##'   'EQUILIBRIUM_PHASES 1      ',
##'   '    Calcite    0.0   1.0  ',
##'   'REACTION_TEMPERATURE 1    ',
##'   '    25.0 75.0 in 51 steps ',
##'   'SELECTED_OUTPUT 1         ',
##'   '    -temperature          ',
##'   '    -si     calcite       '
##'   )
##' 
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetSelectedOutputFileName(1))
##' }
##' 
phrGetSelectedOutputFileName <-
function(nuser) {
  invisible(.Call("getSelectedOutputFileName", as.integer(nuser), PACKAGE = .packageName))
}



##' Set the name of the selected_output file.
##' 
##' Sets the name of the selected_output file. The default value is selected_{nuser}.0.out.
##'
##' The selected_output file must be turned on using the \code{\link{phrSetSelectedOutputFileOn}} function.
##' 
##' @export phrSetSelectedOutputFileName
##' @useDynLib phreeqc
##' @param nuser the user number specified in the SELECTED_OUTPUT block.
##' @param filename the name of the selected_output file.
##' @return NULL
##' @family Selected Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite at various temperatures
##' # and displays the name of the selected_output file.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetSelectedOutputFileOn(1, TRUE)
##' phrSetSelectedOutputFileName(1, "ex2.sel")
##'
##' input <- c( 
##'   'SOLUTION 1 Pure water     ',
##'   'EQUILIBRIUM_PHASES 1      ',
##'   '    Calcite    0.0   1.0  ',
##'   'REACTION_TEMPERATURE 1    ',
##'   '    25.0 75.0 in 51 steps ',
##'   'SELECTED_OUTPUT 1         ',
##'   '    -temperature          ',
##'   '    -si     calcite       '
##'   )
##' 
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetSelectedOutputFileName(1))
##' }
##' 
phrSetSelectedOutputFileName <-
function(nuser, filename) {
  invisible(.Call("setSelectedOutputFileName", as.integer(nuser), as.character(filename), PACKAGE = .packageName))
}



##' Set selected_output file on/off.
##' 
##' Sets the output file switch on or off. This switch controls whether
##' or not phreeqc writes to the output file. This is the output normally
##' generated when phreeqc is run. The initial setting is off.
##' 
##' @export phrSetSelectedOutputFileOn
##' @useDynLib phreeqc
##' @param nuser the user number specified in the SELECTED_OUTPUT block.
##' @param value if TRUE, writes output to the the selected_output file.
##' @return NULL
##' @family Selected Output
##' @references \url{ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/IPhreeqc.pdf}
##' @examples
##' 
##' # This example equilibrates pure water with calcite at various temperatures
##' # and displays the name of the selected_output file.
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetSelectedOutputFileOn(1, TRUE)
##' phrSetSelectedOutputFileName(1, "ex2.sel")
##'
##' input <- c( 
##'   'SOLUTION 1 Pure water     ',
##'   'EQUILIBRIUM_PHASES 1      ',
##'   '    Calcite    0.0   1.0  ',
##'   'REACTION_TEMPERATURE 1    ',
##'   '    25.0 75.0 in 51 steps ',
##'   'SELECTED_OUTPUT 1         ',
##'   '    -temperature          ',
##'   '    -si     calcite       '
##'   )
##' 
##' 
##' if (is.null(phrRunString(input))) {
##'   cat("see", phrGetSelectedOutputFileName(1))
##' }
##' 
phrSetSelectedOutputFileOn <-
function(nuser, value) {
  invisible(.Call("setSelectedOutputFileOn", as.integer(nuser), as.logical(value), PACKAGE = .packageName))
}



##' @name phreeqc.dat
##' @title The phreeqc.dat database
##' @description phreeqc.dat is a phreeqc database file derived from PHREEQE,
##' which is consistent with wateq4f.dat, but has a smaller set of elements and
##' aqueous species. The database has been reformatted for use by
##' \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage phreeqc.dat  # phrLoadDatabaseString(phreeqc.dat)
##' @keywords dataset 
NULL



##' @name ex15.dat
##' @title The ex15.dat database
##' @description ex15.dat is a database used by example 15 (\code{\link{ex15}}).
##' The database has been reformatted for use by
##' \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage ex15.dat  # phrLoadDatabaseString(ex15.dat)
##' @keywords dataset 
NULL



##' @name Amm.dat
##' @title The Amm.dat database.
##' @description Amm.dat is the same as phreeqc.dat, except that ammmonia redox
##' state has been decoupled from the rest of the nitrogen system; that is,
##' ammonia has been defined as a separate component. The database has been
##' reformatted for use by \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage Amm.dat  # phrLoadDatabaseString(Amm.dat)
##' @keywords dataset 
NULL



##' @name wateq4f.dat
##' @title The wateq4f.dat database.
##' @description wateq4f.dat is a database derived from WATEQ4F. The database
##' has been reformatted for use by \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage wateq4f.dat  # phrLoadDatabaseString(wateq4f.dat)
##' @keywords dataset 
NULL



##' @name llnl.dat
##' @title The llnl.dat database.
##' @description llnl.dat is a database derived from databases for EQ3/6 and
##' Geochemist's Workbench that uses thermodynamic data compiled by the
##' Lawrence Livermore National Laboratory. The database has been reformatted
##' for use by \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @usage llnl.dat  # phrLoadDatabaseString(llnl.dat)
##' @keywords dataset 
NULL



##' @name minteq.dat
##' @title The minteq.dat database.
##' @description minteq.dat is a database derived from the databases for the
##' program MINTEQA2. The database has been reformatted for use by
##' \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage minteq.dat  # phrLoadDatabaseString(minteq.dat)
##' @keywords dataset 
NULL



##' @name minteq.v4.dat
##' @title The minteq.v4.dat database.
##' @description minteq.v4.dat is a database derived from MINTEQA2 version 4.
##' The database has been reformatted for use by
##' \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage minteq.v4.dat  # phrLoadDatabaseString(minteq.v4.dat)
##' @keywords dataset 
NULL



##' @name pitzer.dat
##' @title The pitzer.dat database.
##' @description pitzer.dat is a database for the specific-ion-interaction model
##' of Pitzer as implemented in PHRQPITZ. The database has been reformatted for
##' use by \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage pitzer  # phrLoadDatabaseString(pitzer.dat)
##' @keywords dataset 
NULL



##' @name sit.dat
##' @title The sit.dat database.
##' @description sit.dat is a database derived from databases for EQ3/6 and
##' Geochemist's Workbench that uses thermodynamic data compiled by the
##' Lawrence Livermore National Laboratory. The database has been reformatted
##' for use by \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage sit.dat  # phrLoadDatabaseString(sit.dat)
##' @keywords dataset 
NULL



##' @name iso.dat
##' @title The iso.dat database.
##' @description iso.dat is a partial implementation of the individual component
##' approach to isotope calculations as described by Thorstenson and Parkhurst.
##' The database has been reformatted for use by
##' \code{\link{phrLoadDatabaseString}}.
##' @docType data
##' @family Databases
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @usage iso.dat  # phrLoadDatabaseString(iso.dat)
##' @keywords dataset 
NULL



##' @name ex1
##' @title Example 1--Speciation Calculation
##' @description This example calculates the distribution of aqueous species in
##' seawater and the saturation state of seawater relative to a set of minerals.
##' To demonstrate how to expand the model to new elements, uranium is added to
##' the aqueous model defined by \code{\link{phreeqc.dat}}.  The example can be
##' run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex1)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex2
##' @title Example 2--Equilibration With Pure Phases
##' @description This example shows how to calculate the solubility and relative
##' thermodynamic stability of two minerals, gypsum and anhydrite. The example
##' can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex2)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex3
##' @title Example 3--Mixing
##' @description This example demonstrates the capabilities of PHREEQC to
##' perform a series of geochemical simulations, with the final simulations
##' relying on results from previous simulations within the same run. The
##' example investigates diagenetic reactions that may occur in zones where
##' seawater mixes with carbonate groundwater. The example is divided into five
##' simulations, labeled part A through part E. (A) Carbonate groundwater is
##' defined by equilibrating pure water with calcite at a of 10-2.0 atm. (B)
##' Seawater is defined by using the major-ion data given in table.9. (C) The
##' two solutions are mixed together in the proportions 70 percent groundwater
##' and 30 percent seawater. (D) The mixture is equilibrated with calcite and
##' dolomite. (E) The mixture is equilibrated with calcite only to investigate
##' the chemical evolution if dolomite precipitation is assumed to be
##' negligible. The example can be run using the \code{\link{phrRunString}}
##' routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex3)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex4
##' @title Example 4--Evaporation and Homogeneous Redox Reactions
##' @description Evaporation is accomplished by removing water from the chemical
##' system. Water can be removed by several methods: (1) water can be specified
##' as an irreversible reactant with a negative reaction coefficient in the
##' REACTION keyword input, (2) the solution can be mixed with pure water which
##' is given a negative mixing fraction in MIX, or (3) "H2O" can be specified as
##' the alternative reaction in EQUILIBRIUM_PHASES keyword input, in which case
##' water is removed or added to the aqueous phase to attain equilibrium with a
##' specified phase. This example uses the first method; the REACTION data block
##' is used to simulate concentration of rainwater by approximately 20-fold by
##' removing 95 percent of the water. The resulting solution contains only about
##' 0.05 kg of water. In a subsequent simulation, the MIX keyword is used to
##' generate a solution that has the same concentrations as the evaporated
##' solution, but has a total mass of water of approximately 1 kg. The example
##' can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex4)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex5
##' @title Example 5--Irreversible Reactions
##' @description This example demonstrates the irreversible reaction
##' capabilities of PHREEQC in modeling the oxidation of pyrite. Oxygen (O2) and
##' NaCl are added irreversibly to pure water in six amounts (0.0, 0.001, 0.005,
##' 0.01, 0.03, and 0.05 mol); the relative proportion of O2 to NaCl in the
##' irreversible reaction is 1.0 to 0.5. Pyrite, calcite, and goethite are
##' allowed to dissolve to equilibrium and the carbon dioxide partial pressure
##' is maintained at 10-3.5 (atmospheric partial pressure). In addition, gypsum
##' is allowed to precipitate if it becomes supersaturated. The example can be
##' run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex5)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex6
##' @title Example 6--Reaction-Path Calculations
##' @description In this example, the precipitation of phases as a result of
##' incongruent dissolution of K-feldspar (microcline) is investigated. Only the
##' four phases originally addressed by Helgeson and others (1969)--K-feldspar,
##' gibbsite, kaolinite, and K-mica (muscovite)--are considered. The
##' thermodynamic data for the phases (PHASES keyword) are derived from Robie
##' and others (1978) and are the same as for test problem 5 in the PHREEQE
##' manual (Parkhurst and others, 1980). The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex6)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex7
##' @title Example 7--Gas-Phase Calculations
##' @description This example demonstrates the capabilities of PHREEQC to model
##' the evolution of gas compositions in equilibrium with a solution with a
##' fixed (total) pressure or a fixed volume of the gas phase. In the case of a
##' fixed-pressure gas phase, a gas bubble forms as soon as the sum of the
##' partial pressures of the component gases exceeds the specified pressure of
##' the gas phase. Once the bubble forms, its volume and composition will vary
##' with the extent of reactions. This case applies to gas bubbles forming in
##' surface water or groundwater at a given depth, where the total pressure is
##' constant. With a fixed-volume gas phase, the aqueous solution is in contact
##' with a head space of a fixed volume, which is typical for a laboratory
##' experiment with a closed bottle. The gas phase always exists in this head
##' space, but its pressure and composition will vary with the reactions.
##' Another way to model gas-liquid reactions in PHREEQC is to maintain a fixed
##' partial pressure by using the EQUILIBRIUM_PHASES data block. This
##' fixed-partial-pressure approach is illustrated in this example by fixing
##' the CO2 pressure for a SOLUTION. The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##' 
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex7)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex8
##' @title Example 8--Surface Complexation
##' @description In all surface complexation models, sorption is a function of
##' both chemical and electrostatic energy as described by the free energy
##' relationship. Sorption is stronger when the Gibbs energy decreases. Thus, a
##' counter-ion that carries a charge opposite to the surface charge tends to be
##' sorbed electrostatically, while a co-ion that carries a charge with the same
##' sign as the surface tends to be rejected. The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' # example 8 requires the selected_output file to be turned on
##' phrSetSelectedOutputFileOn(1, TRUE)
##' phrRunString(ex8)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex9
##' @title Example 9--Kinetic Oxidation of Dissolved Ferrous Iron With Oxygen
##' @description Kinetic rate expressions can be defined in a completely general
##' way in PHREEQC by using Basic statements in the RATES data block. The rate
##' expressions can be used in batch-reaction and transport calculations with
##' the KINETICS data block. For transport calculations (ADVECTION or TRANSPORT),
##' kinetic reactions can be defined cell by cell by the number range following
##' the KINETICS keyword (KINETICS m-n). The rate expressions are integrated with
##' an embedded (up to) 5th-order Runge-Kutta-Fehlberg algorithm, or with a
##' stiff, variable-order, variable-step multistep solver (Cohen and Hindmarsh,
##' 1996). Equilibrium is calculated before a kinetic calculation is initiated
##' and again when a kinetic reaction increment is added. Equilibrium includes
##' solution species equilibrium; exchange-, equilibrium-phase-, solid-solution-,
##' and surface-assemblage equilibrium; and gas-phase equilibrium. A check is
##' performed to ensure that the difference between estimates of the integrated
##' rate over a time interval is smaller than a user-defined tolerance. If the
##' tolerance is not satisfied, then the integration over the time interval is
##' automatically restarted with a smaller time interval. The example can be run
##' using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex9)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex10
##' @title Example 10--Aragonite-Strontianite Solid Solution
##' @description PHREEQC has the capability to model multicomponent ideal and
##' binary nonideal solid solutions. For ideal solid solutions, the activity of
##' each end member solid is equal to its mole fraction. For nonideal solid
##' solutions, the activity of each end member is the product of the mole
##' fraction and an activity coefficient, which is determined from the mole
##' fraction and Guggenheim excess free-energy parameters. Example 10 considers
##' an aragonite (CaCO3)-strontianite (SrCO3) solid solution and demonstrates
##' how the composition of the solid solution and the aqueous phase change as
##' strontium carbonate is added to an initially pure calcium carbonate system.
##' The example can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex10)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex11
##' @title Example 11--Transport and Cation Exchange
##' @description The following example simulates the chemical composition of the
##' effluent from a column containing a cation exchanger (Appelo and Postma,
##' 2005). Initially, the column contains a sodium-potassium-nitrate solution
##' in equilibrium with the exchanger. The column is flushed with three pore
##' volumes of calcium chloride solution. Calcium, potassium, and sodium react
##' to equilibrium with the exchanger at all times. The problem is run two
##' ways--by using the ADVECTION data block, which models only advection, and by
##' using the TRANSPORT data block, which simulates advection and dispersive
##' mixing. The example can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex11)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex12
##' @title Example 12--Advective and Diffusive Flux of Heat and Solutes
##' @description This example demonstrates the capability of PHREEQC to
##' calculate transient transport of heat and solutes in a column or along a 1D
##' flowline. A column is initially filled with a dilute KCl solution at 25
##' degrees C in equilibrium with a cation exchanger. A KNO3 solution then
##' advects into the column and establishes a new temperature of 0 degrees C.
##' Subsequently, a sodium chloride solution at 24 degrees C is allowed to
##' diffuse from both ends of the column, assuming no heat is lost through the
##' column walls. At one end, a constant boundary condition is imposed, and at
##' the other end, the final cell is filled with the sodium chloride solution
##' and a closed boundary condition is prescribed. For the column end with a
##' constant boundary condition, an analytical solution is compared with PHREEQC
##' results, for unretarded Cl- (R = 1.0) and retarded Na+ and temperature
##' (R = 3.0). Finally, the second-order accuracy of the numerical method is
##' verified by increasing the number of cells by a factor of three and
##' demonstrating a decrease in the error of the numerical solution by
##' approximately one order of magnitude relative to the analytical solution.
##' The example can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex12)
##' phrGetOutputStrings()
##' 
NULL



##' @rdname ex13
##' @name ex13a
##' @aliases ex13a ex13b ex13c
##' @title Example 13--Aragonite-Strontianite Solid Solution
##' @description PHREEQC has the capability to model multicomponent ideal and
##' binary nonideal solid solutions. For ideal solid solutions, the activity of
##' each end member solid is equal to its mole fraction. For nonideal solid
##' solutions, the activity of each end member is the product of the mole
##' fraction and an activity coefficient, which is determined from the mole
##' fraction and Guggenheim excess free-energy parameters. Example 10 considers
##' an aragonite (CaCO3)-strontianite (SrCO3) solid solution and demonstrates
##' how the composition of the solid solution and the aqueous phase change as
##' strontium carbonate is added to an initially pure calcium carbonate system.
##' The example can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex13a)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex14
##' @title Example 14--Advective Transport, Cation Exchange, Surface
##' Complexation, and Mineral Equilibria
##' @description This example uses the phase-equilibrium, cation-exchange, and
##' surface-complexation reaction capabilities of PHREEQC in combination with
##' advective-transport capabilities to model the evolution of water in the
##' Central Oklahoma aquifer. The geochemistry of the aquifer has been described
##' by Parkhurst and others (1996). Two predominant water types occur in the
##' aquifer: a calcium magnesium bicarbonate water with pH in the range of 7.0
##' to 7.5 in the unconfined part of the aquifer and a sodium bicarbonate water
##' with pH in the range of 8.5 to 9.2 in the confined part of the aquifer. In
##' addition, marine-derived sodium chloride brines exist below the aquifer and
##' presumably in fluid inclusions and dead-end pore spaces within the aquifer.
##' Large concentrations of arsenic, selenium, chromium, and uranium occur
##' naturally within the aquifer. Arsenic is associated almost exclusively with
##' the high-pH, sodium bicarbonate water type. The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex14)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex15
##' @title Example 15--1D Transport: Kinetic Biodegradation, Cell Growth, and
##' Sorption
##' @description A test problem for advective-dispersive-reactive transport was
##' developed by Tebes-Stevens and Valocchi (1997) and Tebes-Stevens and others
##' (1998). Although based on relatively simple speciation chemistry, the
##' solution to the problem demonstrates several interacting chemical processes
##' that are common to many environmental problems: bacterially mediated
##' degradation of an organic substrate; bacterial cell growth and decay; metal
##' sorption; and aqueous speciation, including metal-ligand complexation. In
##' this example, the test problem is solved with PHREEQC, which produces
##' results almost identical to those of Tebes-Stevens and Valocchi (1997) and
##' Tebes-Stevens and others (1998). The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' # this example takes longer than 5 seconds
##' phrLoadDatabaseString(ex15.dat)
##' phrSetOutputStringsOn(TRUE)
##' \dontrun{phrRunString(ex15)}
##' phrGetOutputStrings()
##' 
NULL



##' @name ex16
##' @title Example 16--Inverse Modeling of Sierra Spring Waters
##' @description This example repeats the inverse modeling calculations of the
##' chemical evolution of spring-water compositions in the Sierra Nevada that
##' are described in a classic paper by Garrels and Mackenzie (1967). The same
##' example is described in the manual for the inverse-modeling program NETPATH
##' (Plummer and others, 1991 and 1994). The example uses two spring-water
##' compositions, one from an ephemeral spring, which is less chemically
##' evolved, and one from a perennial spring, which probably has had a longer
##' residence time in the subsoil. The differences in composition between the
##' ephemeral spring and the perennial spring are assumed to be caused by
##' reactions between the water and the minerals and gases it contacts. The
##' object of inverse modeling in this example is to find sets of minerals and
##' gases that, when reacted in appropriate amounts, account for the differences
##' in composition between the two solutions. The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex16)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex17
##' @title Example 17--Inverse Modeling With Evaporation
##' @description Evaporation is handled in the same manner as other
##' heterogeneous reactions for inverse modeling. To model evaporation (or
##' dilution), it is necessary to include a phase with the composition "H2O".
##' The important concept in modeling evaporation is the water mole-balance
##' equation (see Parkhurst and Appelo, 1999, "Equations and Numerical Method
##' for Inverse Modeling"). The moles of water in the initial solutions times
##' their mixing fractions, plus water gained or lost by dissolution or
##' precipitation of phases, plus water gained or lost through redox reactions,
##' must equal the moles of water in the final solution. The equation is still
##' approximate because it does not include the moles of water gained or lost in
##' hydrolysis and complexation reactions in the solutions. The results of
##' inverse modeling are compared with a forward model using Pitzer equations to
##' calculate the sequence of salts that precipitate during evaporation. The
##' example can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(pitzer.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex17)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex18
##' @title Example 18--Inverse Modeling of the Madison Aquifer
##' @description In this example, inverse modeling, including isotope
##' mole-balance modeling, is applied to the evolution of water in the Madison
##' aquifer in Montana. Plummer and others (1990) used mole-balance modeling to
##' quantify the extent of dedolomitization at locations throughout the aquifer.
##' In the dedolomitization process, anhydrite dissolution causes the
##' precipitation of calcite and dissolution of dolomite. Additional reactions
##' identified by mole-balance modeling include sulfate reduction, cation
##' exchange, and halite and sylvite dissolution (Plummer and others, 1990). Del
##' 13C and del 34S data were used to corroborate the mole-balance models and
##' carbon-14 was used to estimate groundwater ages (Plummer and others, 1990).
##' Initial and final water samples were selected from a flow path that extends
##' from north-central Wyoming northeast across Montana (Plummer and others,
##' 1990, flow path 3). This pair of water samples was selected specifically
##' because it was one of the few pairs that showed a relatively large
##' discrepancy between previous mole-balance approaches and the mole-balance
##' approach of PHREEQC, which includes uncertainties; results for most sample
##' pairs were not significantly different between the two approaches. In
##' addition, this pair of samples was selected because it was modeled in detail
##' in Plummer and others (1990) to determine the sensitivity of mole-balance
##' results to various model assumptions and was used as an example in the
##' NETPATH manual (Plummer and others, 1994, example 6). Results of PHREEQC
##' calculations are compared to NETPATH calculations. This example is also
##' discussed in Parkhurst (1997). The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex18)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex19
##' @title Example 19--Modeling Cd+2 Sorption With Linear, Freundlich, and
##' Langmuir Isotherms, and With a Deterministic Distribution of Sorption Sites
##' for Organic Matter, Clay Minerals, and Iron Oxyhydroxides
##' @description Sorption of heavy metals and organic pollutants on natural
##' materials can be described by linear, Freundlich, or Langmuir isotherms. All
##' three isotherms can be calculated by PHREEQC, as shown in this example for
##' Cd+2 sorbing on a loamy soil (Christensen, 1984; Appelo and Postma, 2005). A
##' more mechanistic approach, also illustrated here, is to model the
##' distribution of Cd+2 over the sorbing components in the soil, in this case,
##' in and on organic matter, clay minerals, and iron oxyhydroxides. The example
##' can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex19)
##' phrGetOutputStrings()
##' 
NULL



##' @rdname ex20
##' @name ex20a
##' @aliases ex20a ex20b
##' @title Example 20--Distribution of Isotopes Between Water and Calcite
##' @description The database \code{\link{iso.dat}} implements the approach to
##' isotope reactions described by Thorstenson and Parkhurst (2002, 2004), in
##' which minor isotopes are treated as individual thermodynamic components. The
##' aqueous and solid species of minor isotopes have slightly different
##' equilibrium constants than those of the major isotopes, which account for
##' fractionation processes. The treatment of isotopes in gases requires a
##' separate species for each isotopic variant of a gas; for example, the
##' isotopic variants of carbon dioxide are CO2, C18OO, C18O2, 13CO2, 13C18OO,
##' and 13C18O2. Similarly, every isotopic variant of a mineral must be included
##' as a component of a solid solution to represent completely the isotopic
##' composition of the solid. The equilibrium constants in iso.dat are derived
##' from empirical fractionation factors, from statistical mechanical theory,
##' or, where no data are available (the most common case), by assuming no
##' fractionation. However, the database is a framework that can be expanded as
##' additional isotopic thermodynamic data become available. The example can be
##' run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(iso.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex20a)
##' phrGetOutputStrings()
##' 
NULL



##' @name ex21
##' @title Example 21--Modeling Diffusion of HTO, 36Cl-, 22Na+, and Cs+ in a
##' Radial Diffusion Cell
##' @description This example illustrates how PHREEQC version 3 can simulate a
##' diffusion experiment, as is now often performed for assessing the properties
##' of a repository for nuclear waste in a clay formation. A sample is cut from
##' a core of clay, enveloped in filters, and placed in a diffusion cell (see
##' Van Loon and others, 2004, for details). Solutions with tracers are
##' circulated at the surfaces of the filters, the tracers diffuse into and out
##' of the clay, and the solutions are sampled and analyzed regularly in time.
##' The concentration changes are interpreted with Fick's diffusion equations to
##' obtain transport parameters for modeling the rates of migration of elements
##' away from a waste repository. Transport in clays is mainly diffusive because
##' of the low hydraulic conductivity, and solutes are further retarded by
##' sorption (cations) and by exclusion from part of the pore space (anions).
##' The example can be run using the \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' # example 21 requires the selected_output file to be turned on
##' phrSetSelectedOutputFileOn(1, TRUE)
##' phrSetOutputStringsOn(TRUE)
##' # this takes longer than 5 seconds
##' \dontrun{phrRunString(ex21)}
##' phrGetOutputStrings()
##' 
NULL



##' @name ex22
##' @title Example 22--Modeling Gas Solubilities: CO2 at High Pressures
##' @description PHREEQC calculates the fugacity coefficient with the
##' Peng-Robinson equation of state (Peng and Robinson, 1976) from the critical
##' pressure and temperature, and the acentric factor of the gas in a gas
##' mixture to obtain the limiting volume and the attraction factor in the Van
##' der Waals equation. The fugacity coefficient is close to 1 when the total
##' pressure of the gas phase is less than about 10 atm, and it can be neglected
##' in the solubility calculation. At higher pressures, the effect can be
##' substantial. At low pressures, the concentration of CO2 increases
##' near-linearly with pressure. At 25 degrees C and pressures higher than 62
##' atm, the concentration increases more gradually because the fugacity
##' coefficient drops rapidly. The example can be run using the
##' \code{\link{phrRunString}} routine.
##' @docType data
##' @family Examples
##' @references \url{http://pubs.usgs.gov/tm/06/a43/pdf/tm6-A43.pdf}
##' @source \url{http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc}
##' @keywords dataset 
##' @examples
##'
##' phrLoadDatabaseString(phreeqc.dat)
##' phrSetOutputStringsOn(TRUE)
##' phrRunString(ex22)
##' phrGetOutputStrings()
##' 
NULL
