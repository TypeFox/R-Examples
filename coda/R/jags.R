"read.jags" <- function (file = "jags.out", start, end, thin, quiet=FALSE) 
{
  if (!is.R()) {
    stop("This function is not yet available in S-PLUS")
  }
  read.coda(file, start, end, thin, quiet)
}

bugs2jags <- function(infile, outfile)
{
  if (!is.R()) {
    stop("This function is not yet available in S-PLUS")
  }
  
  ## Convert S-style data for WinBUGS into the R dump format
  ## used by JAGS.
  bugs.dat <- dget(infile)
  for (bugs.variable.name in names(bugs.dat)) {
    if(!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
      ## Manually reverse order of dimensions of arrays
      dim(bugs.dat[[bugs.variable.name]]) <-
        rev(dim(bugs.dat[[bugs.variable.name]]))
      ## Then transpose
      bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
    }
    assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
  }
    dump(names(bugs.dat), file=outfile)
}
