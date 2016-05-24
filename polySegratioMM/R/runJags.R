`runJags` <-
function(jags.control, jags="jags", quiet = FALSE,
                     cmd.file=paste(jags.control$stem,".cmd",sep=""),
                     timing=TRUE)
{
  if (class(jags.control) != "jagsControl")
    stop("'jags.control' must be of class 'jagsControl'")
  
  if (! file.exists(cmd.file))  ## maybe a bit redundant
    writeControlFile(jags.control, cmd.file)
  
  ## adapted from JAGScall in 'bayesmix'

  start.time <- date()
  if(timing) {
    ptm <- proc.time()
  }
  
  if (.Platform$OS.type == "windows")
  ##  exit <- system(paste(jags, cmd.file))
    exit <- shell(paste(jags, cmd.file))  # not sure why but system not working
  else exit <- system(paste(jags, "< ", cmd.file ," > /dev/null"),
                      ignore.stderr = quiet)
  if (exit)
    stop("System call not successful")
  if (! file.exists((paste(jags.control$stem,"CODAchain1.txt",sep=""))))
    exit <- 1

  res <- list(jags.control=jags.control, exit=exit, cmd.file=cmd.file,
              start.time=start.time, end.time=date(), call=match.call())
  if(timing) {
    res$elapsed.time <- proc.time() - ptm
  }

  class(res) <- "runJags"
  return(res)
}

