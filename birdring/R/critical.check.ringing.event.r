
critical.check.ringing.event <- function(dat, id="birdID"){
  # dat EURING data read by read.euring2000plus
  # id = name of the variable containing the individual identifier (normally a combination of ringing scheme and ring number)
  if (sum(is.element(names(dat), "metal.ring.info"))==0)  {
    stop('The variable metal.ring.info is missing.')
  }
  if (sum(is.element(names(dat), id))==0)  {
    stop(paste('The variable', id,   'is missing.'))
  }
  
  eventbybird <- tapply(dat$metal.ring.info, dat[,id],  function(x) ifelse(sum(is.element(x, c(1:3)))!=1, FALSE, TRUE))
  dat$ringing.event <- eventbybird[match(dat[,id], levels(factor(dat[,id])))]
  return(dat)
}
