### These are for dumpping temporary within some iterations of long MCMC.

### Get the specific function according to the options.
get.my.dump <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel is not found.")
  }
  ret <- eval(parse(text = paste("my.dump.",
                                 parallel[1], sep = "")))
  assign("my.dump", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.dump().


### For lapply.
my.dump.lapply <- function(iter, list = NULL, envir = parent.frame()){
  if(.CF.DP$dump){
    if((iter %% .CF.DP$iter) == 0){
      if(is.null(list)){
        list <- ls(envir = envir)
      }

      cat("dump start:", iter, "\t", date(), "\n")

      file.data <- paste(.CF.DP$prefix.dump, "data_",
                         Sys.getpid(), ".rda", sep = "")
      save(list = list, file = file.data, envir = envir)

      .cubfitsEnv$iter <- iter
      if(exists(".Random.seed", envir = .GlobalEnv)){
        .cubfitsEnv$Random.seed <- .GlobalEnv$.Random.seed
      }
      file.env <- paste(.CF.DP$prefix.dump, "env_",
                        Sys.getpid(), ".rda", sep = "")
      save(list = ls(all.names = TRUE, envir = .cubfitsEnv),
           file = file.env, envir = .cubfitsEnv)

      cat("dump end:", iter, "\t", date(), "\n")
    }
  }

  invisible()
} # End of my.dump.lapply().

### For mclapply.
my.dump.mclapply <- my.dump.lapply

### For task pull.
my.dump.task.pull <- function(iter, list = NULL, envir = parent.frame()){
  if(pbdMPI::comm.rank() == 0L){
    my.dump.lapply(iter, list = list, envir = envir)
  }
  invisible()
} # End of my.dump.task.pull().

### For pbdLapply.
my.dump.pbdLapply <- function(iter, list = NULL, envir = parent.frame()){
  if(pbdMPI::comm.rank() == 0L){
    my.dump.lapply(iter, list = list, envir = envir)
  }
  invisible()
} # End of my.dump.pbdLapply().

