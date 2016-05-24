### For print some information within MCMC.
my.verbose <- function(verbose, iter, report){
  if(verbose){
    if(iter != 0){
      if((iter %% report) == 0){
        .cubfitsEnv$my.cat("iter:", iter, "\t", date(), "\n")
      }
      if((iter %% .CF.DP$report.proc) == 0){
        total.time <- proc.time() - .cubfitsEnv$start.time
        total.time <- rbind(total.time, total.time / iter)
        rownames(total.time)[2] <- "avg.time"
        .cubfitsEnv$my.print(total.time)
      }
    } else{
      .cubfitsEnv$my.cat("pid:\t\t", Sys.getpid(), "\n")
      .cubfitsEnv$my.cat("start:\t\t", date(), "\n")
      .cubfitsEnv$start.time <- proc.time()
    }
  }
  invisible()
} # End of my.verbose().
