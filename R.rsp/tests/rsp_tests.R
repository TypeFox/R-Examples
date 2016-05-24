library("R.rsp")
library("R.utils") # Arguments

verbose <- Arguments$getVerbose(TRUE)

# Change option (used in one or more of the tests)
oopts <- options(papersize="letter")

path <- system.file(package="R.rsp")
path <- file.path(path, "rsp_tests")

# Find all RSP files with matching "truth" files
pathnames <- list.files(path=path, pattern="[.]rsp$", full.names=TRUE)
pathnames <- pathnames[file_test("-f", pathnames)]

# Reference ("truth") files
pathnamesR <- gsub("[.]rsp$", "", pathnames)
keep <- file_test("-f", pathnamesR)
pathnames <- pathnames[keep]
pathnamesR <- pathnamesR[keep]

for (kk in seq_along(pathnames)) {
  pathname <- pathnames[kk]
  pathnameR <- pathnamesR[kk]
  verbose && enter(verbose, sprintf("RSP file #%d ('%s') of %d", kk, pathname, length(pathnames)))

  rc <- rcode(file=pathname)

  rs <- rstring(file=pathname)
  s <- as.character(rs)
  sR <- readChar(pathnameR, nchars=1e6)
  if (length(sR) == 0L) sR <- ""
  sR <- gsub("\r\n", "\n", sR, fixed=TRUE)
  # If there is a '<EOF>' string, drop it and everything beyond
  sR <- sub("<EOF>.*", "", sR)

  # Compare
  res <- all.equal(s, sR, check.attributes=FALSE)
  discrepancy <- !isTRUE(res);
  if (discrepancy) {
    if (verbose) {
      enter(verbose, "Detected discrepancy of output and expected output")
      cat(verbose, "Output:")
      printf(verbose, ">>>%s<<<\n", gsub("\n", "\\n", s, fixed=TRUE))
      ruler(verbose)
      cat(verbose, "Expected output:")
      printf(verbose, ">>>%s<<<\n", gsub("\n", "\\n", sR, fixed=TRUE))
      ruler(verbose)
      cat(verbose, "Difference:")
      fs <- tempfile("s")
      fsR <- tempfile("sR")
      cat(s, file=fs)
      cat(sR, file=fsR)
      res <- system2("diff", args=c("-bw", shQuote(fs), shQuote(fsR)),
                                             stdout=TRUE, stderr=TRUE)
      cat(verbose, res)
      file.remove(c(fs, fsR))
      ruler(verbose)
      cat(verbose, "Reason for not being equal:")
      print(verbose, res)
      ruler(verbose)
      exit(verbose)
    }
  }

  verbose && exit(verbose)
}

# Restore options
options(oopts)

if (discrepancy) {
  stop("Detected a discrepancy above.")
}
