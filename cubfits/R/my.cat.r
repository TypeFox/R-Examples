### These are for cat within some iterations of long MCMC.

### Get the specific function according to the option.
get.my.cat <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel is not found.")
  }
  ret <- eval(parse(text = paste("my.cat.",
                                 parallel[1], sep = "")))
  assign("my.cat", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.cat().


### For lapply.
my.cat.lapply <- function(..., file = "", sep = " ", fill = FALSE,
    labels = NULL, append = FALSE){
  cat(..., file = file, sep = sep, fill = fill, labels = labels,
      append = append)
} # End of my.cat.lapply().

### For mclapply.
my.cat.mclapply <- my.cat.lapply

### For task pull.
my.cat.task.pull <- function(..., file = "", sep = " ", fill = FALSE,
    labels = NULL, append = FALSE){
  if(pbdMPI::comm.rank() == 0L){
    cat(..., file = file, sep = sep, fill = fill, labels = labels,
        append = append)
  }
  invisible()
} # End of my.cat.task.pull().

### For pbdLapply. 
my.cat.pbdLapply <- function(..., file = "", sep = " ", fill = FALSE,
    labels = NULL, append = FALSE){
  if(pbdMPI::comm.rank() == 0L){
    cat(..., file = file, sep = sep, fill = fill, labels = labels,
        append = append)
  }
  invisible()
} # End of my.cat.pbdLapply().
