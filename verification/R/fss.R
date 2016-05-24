
fss <- function(obs, pred, w = 0, FUN = mean, ...){
### compare matrixes of forecast of observed values and forecast.
### values can be calcuated using different windows.

### with a window size of 0, obs is returned.
obs.matrix <- matrix.func(DAT = obs, w = w, FUN = FUN) 
  
### with a window size of 0, obs is returned.
frcs.matrix <- matrix.func(DAT = pred, w = w, FUN = FUN)   

if(nrow(obs)!= nrow(pred) &  ncol(obs)!=  nrow(obs) ) stop("Observation matrix and forecast matrix different sizes")
  
n   <- prod(dim(obs.matrix))  ### number of gridpoints

N   <- sum((obs.matrix-frcs.matrix)^2, na.rm = TRUE)/n ### numerator
D   <- (sum(obs.matrix^2, na.rm = TRUE) +sum(frcs.matrix^2, na.rm = TRUE))/n ### denominator

FSS <- 1 - N/D  
return(FSS)
}

matrix.func <- function(DAT, w = 0, FUN = mean, ...){

### w is the '' radius'' of window.  eg. w = 2, defines a 5 by 5 square
  
### define function
FUN <- match.fun(FUN)

### define output dimension
II <- nrow(DAT) - 2*w   ## output row dimension
JJ <- ncol(DAT) - 2*w  

if(JJ<=0|II <= 0) {stop("The window exceeds the size of the observation" ) } 
OUT <- matrix(NA, nrow= II, ncol = JJ)

for(i in 1:II){
for(j in 1:JJ){
sub <- DAT[ i :(i + 2*w ),
             j :(j + 2*w ) ] # subset data

OUT[i,j] <- FUN(sub,...)
}  ## close J
}  ## close I

return(OUT)
}  ## close function

 
