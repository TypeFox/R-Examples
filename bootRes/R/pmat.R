pmat <- function(x, start = -6, end = 9, vnames = NULL) {
  years <- unique(x[, 1])
  n <- length(years)
  no.vars <- dim(x)[2] - 2 # number of variables
  months <- paste(c(rep("prev.", 12), rep("curr.", 12)), rep(c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), 2), sep = "")
  month.ids <- c(-1:-12, 1:12)
  used.months <- months[which(month.ids == start):which(month.ids == end)]
  no.months <- length(used.months)

  ## check for specified variable names, else default to colnames of
  ## x, if not present to V1, V2 etc.
  if (is.null(vnames) || length(unique(vnames)) != no.vars) {
    if (!is.null(colnames(x))) {
      vnames <- colnames(x)[3:(no.vars+2)]		
    } else {
      vnames <- paste(rep("V", no.vars), 1:no.vars, sep = "")
    }	
  }

  ## create unique names for variables
  vnames.mat <- matrix(NA, nrow = no.months, ncol = no.vars) 
  for (i in 1:no.vars) {
    vnames.mat[, i] <- paste(vnames[i], ".", used.months, sep = "")
  }
  
  m <- matrix(NA, nrow = no.months*no.vars, ncol = n - 1)
  colnames(m) <- years[-1]
  rownames(m) <- as.vector(vnames.mat)
  
  for (i in 2:n) {
    if (start < 0) {
      start.with <- which(x[, 1] == years[i - 1])[abs(start)] # start month in previous year
    } else {
      start.with <- which(x[, 1] == years[i])[start] # start month in current year
    }
    for (k in 1:no.vars) { # loop through variables
      for (j in 1:no.months) { # loop through months
        m[(j + (no.months*(k-1))), (i - 1)] <- x[(start.with + j - 1), 2+k]
      }
    }
  }
  
  pmatrix <- as.data.frame(t(m))
  attributes(pmatrix)$npar <- no.vars
  
  pmatrix
}

