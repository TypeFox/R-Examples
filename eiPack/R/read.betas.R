read.betas <- function(rows, columns, units, dir = NULL, ret.mcmc = TRUE) {
  "%w/o%" <- function(x,y) x[!x %in% y] #--  x without y

  idx <- expand.grid(rows, columns, units)
  idx <- apply(idx, 1, paste, collapse = ".")
  files <- paste(idx, ".txt.gz", sep = "")
  
  if (is.null(dir)) {
    dir <- getwd()
  }
  else { 
    current.dir <- getwd()
    setwd(dir)
  }
  check <- system("ls", intern = TRUE)
  if (!all(files %in% check)) 
    stop(paste(paste(paste(files %w/o% check, ".txt.gz", sep = ""), collapse = " "),
         " files not found in ", dir, sep = ""))
  readin <- list()
  for (ii in files)
    readin[[ii]] <- read.table(gzfile(files[ii]), as.is = TRUE, header=TRUE)[[1]]
  
  vidx <- names(readin)
  if (ret.mcmc) {
    dat <- as.mcmc(readin)
    class(dat) <- c(class(dat), "eiMD.beta")
  }
  else { 
    idx1 <- list(rows = rows, columns = columns, units = units) 
    idx1[[4]] <- 1:length(readin[[1]])
    names(idx1) <- c("rows", "columns", "units", "sims")
    dat <- array(t(readin), dim = sapply(idx, length),
                 dimnames = idx)
    class(dat) <- c("array", "eiMD.beta")
  }
  dat
}
