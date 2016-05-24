# DRAW RANDOM VALUES FOR THE PARAMETERS OF THE PSTEP
# TIME 

psteppar<- function(pstep){

psteptime <- runif(1, pstep$time_min, pstep$time_max)

return(psteptime)

} # END FUNCTION