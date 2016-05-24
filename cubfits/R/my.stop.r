### These are for stop within some iterations of long MCMC.

### Get the specific function according to the option.
get.my.stop <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel is not found.")
  }
  ret <- eval(parse(text = paste("my.stop.",
                                 parallel[1], sep = "")))
  assign("my.stop", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.stop().


### For lapply.
my.stop.lapply <- function(..., call. = TRUE, domain = NULL){
  stop(..., call. = call., domain = domain)
} # End of my.stop.lapply().

### For mclapply.
my.stop.mclapply <- my.stop.lapply

### For task pull.
my.stop.task.pull <- function(..., call. = TRUE, domain = NULL){
  pbdMPI::comm.stop(..., call. = call., domain = domain)
  invisible()
} # End of my.stop.task.pull().

### For pbdLapply. 
my.stop.pbdLapply <- function(..., call. = TRUE, domain = NULL){
  pbdMPI::comm.stop(..., call. = call., domain = domain)
  invisible()
} # End of my.stop.pbdLapply().

