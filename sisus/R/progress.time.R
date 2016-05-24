progress.time <-
function# returns the time taken since the time started
### internal function for sisus
(time.start
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  time.sofar = proc.time()[3] - time.start;   # calculate time elapsed
  return(time.sofar);
  ### internal variable
}
