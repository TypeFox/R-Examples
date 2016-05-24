### These are for print within some iterations of long MCMC.

### Get the specific function according to the option.
get.my.print <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel is not found.")
  }
  ret <- eval(parse(text = paste("my.print.",
                                 parallel[1], sep = "")))
  assign("my.print", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.print().


### For lapply.
my.print.lapply <- function(x, ...){
  print(x, ...)
} # End of my.print.lapply().

### For mclapply.
my.print.mclapply <- my.print.lapply

### For task pull.
my.print.task.pull <- function(x, ...){
  if(pbdMPI::comm.rank() == 0L){
    print(x, ...)
  }
  invisible()
} # End of my.print.task.pull().

### For pbdLapply. 
my.print.pbdLapply <- function(x, ...){
  if(pbdMPI::comm.rank() == 0L){
    print(x, ...)
  }
  invisible()
} # End of my.print.pbdLapply().
