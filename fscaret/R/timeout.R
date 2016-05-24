# Function eval_fork() of opencpu package modified
# Original code by Jeroen Ooms <jeroen.ooms at stat.ucla.edu> of OpenCPU package
# 

timeout <- function(..., seconds) {
  
  fork_to_check <- parallel::mcparallel(
    {eval(...)},
    silent = FALSE)
  
  # call mccollect to wait "seconds" for returning result of mcparallel.
  my_result <- parallel::mccollect(fork_to_check, wait = FALSE, timeout = seconds)
  # If my_result is returned kill fork
  tools::pskill(fork_to_check$pid, tools::SIGKILL)
  tools::pskill(-1 * fork_to_check$pid, tools::SIGKILL)
  
  # kill the fork of forks if they were spawned
  parallel::mccollect(fork_to_check, wait = FALSE)
  # If the function mccollect had NULL (timedout), make stop
  if (is.null(my_result))
    stop("Time limit has reached!")
  
  # remove list format
  my_result <- my_result[[1]]
  
  # return result
  return(my_result)
}