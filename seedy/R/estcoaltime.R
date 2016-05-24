estcoaltime <-
function(bottlesize, popsize, bottletimes, obstime) {
  
  backwardbottles <- sort(obstime-bottletimes, decreasing=TRUE)
  summer <- 0
  for (j in 1:(obstime-1)) {
    part <- (1-1/popsize)^(j-1-sum(backwardbottles<j))*(1-1/bottlesize)^sum(backwardbottles<j)*j
    if (j %in% backwardbottles) {
      part <- part/bottlesize
    } else {
      part <- part/popsize
    }
    summer <- summer+part
  }
  summer <- summer + (1-1/popsize)^(j-1-sum(backwardbottles<j))*(1-1/bottlesize)^sum(backwardbottles<j)*obstime
  return(summer)
}
