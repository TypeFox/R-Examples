# Generic function for bringing up the visualization tool
showDSA <- function(x, javacmd=getOption("javacmd"), quietly=FALSE, ...) {
  UseMethod('showDSA')
}

# Visualize a partDSA/dsa object
showDSA.dsa <- function(x, javacmd=getOption("javacmd"), quietly=FALSE, ...) {
  # Dump visualization information to a temp file
  tfile <- tempfile(pattern="dsa")
  on.exit(unlink(tfile))
  dumpDSA(x, file=tfile)
  showDSA(tfile, javacmd=javacmd, quietly=quietly, wait=TRUE, ...)
}

# Visualize a partDSA/dsa object that has been dumped to a file via dumpDSA
showDSA.character <- function(x, javacmd=getOption("javacmd"),
                              quietly=FALSE, wait=FALSE, ...) {
  if (is.null(javacmd))
    javacmd <- "java"

  # Get the path to the jar file included in the "java" directory
  jardir <- system.file("java", package="partDSA")
  jarpath <- Sys.glob(file.path(jardir, "dsa-*.jar"))
  n <- length(jarpath)
  if (n == 0) {
    stop("unable to find the jar file")
  } else if (n > 1) {
    warning("found multiple jar files: using ", jarpath[n])
  }

  # Execute the Java interpreter to display the visualization information
  cmdv <- c(javacmd, '-jar', jarpath[n], x)
  cmd <- qcmd(cmdv)
  if (! quietly) {
    message(sprintf("Executing command: %s", cmd))
    if (wait) {
      message("The R session will block until the Java program exits")
    }
  }

  # Run the Java visualization program
  stat <- system(cmd, wait=wait)

  if (stat != 0) {
    warning("Got non-zero status from 'system' function when executing Java program")
    warning("Java must be installed to use the showDSA function")
    if (javacmd == "java") {
      warning("The 'java' command must be in a directory in PATH",
              " or its path specified via the 'javacmd' argument")
    }
  } else if (! quietly) {
    if (wait) {
      message("The Java program has exited")
    } else {
      message("The Java program is running in another process")
    }
  }
}

# Turn an argv into a command string
qcmd <- function(cmdv) {
  win.quote <- function(arg) {
    if (! nzchar(arg)) {
      '""'
    } else if (length(grep('[[:space:]"]', arg)) == 0) {
      arg
    } else {
      q <- '"'
      nbs <- 0
      v <- strsplit(arg, split='')[[1]]
      for (c in v) {
        if (c == '\\') {
          q <- paste(q, c, sep='')
          nbs <- nbs + 1
        } else if (c == '"') {
          q <- paste(q, paste(rep('\\', nbs + 1), collapse=''), c, sep='')
          nbs <- 0
        } else {
          q <- paste(q, c, sep='')
          nbs <- 0
        }
      }
      paste(q, paste(rep('\\', nbs), collapse=''), '"', sep='')
    }
  }

  unix.quote <- function(arg) {
    q <- "'"
    v <- strsplit(arg, split='')[[1]]
    for (c in v) {
      if (c == "'") {
        c <- "'\\''"
      }
      q <- paste(q, c, sep='')
    }
    paste(q, "'", sep='')
  }

  fun <- if (.Platform$OS.type == 'windows') win.quote else unix.quote
  paste(lapply(cmdv, fun), collapse=' ')
}
