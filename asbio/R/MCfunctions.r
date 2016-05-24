#-------------------------------------------------------------------------#
# Arguments:
# mat = A symmetric matrix.
# pow = The power that the matrix is to be raised to.
#
mat.pow <- function(mat, pow){
res <- mat 
for(i in 2 : pow){
res <- res %*% mat
}
res
}
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# Arguments:
# start = starting value. 
# length = length of chain.
# T = transition matrix.
# states = Number of states.

MC <- function(T, start, length){ 
T <- as.matrix(T)
states <- ncol(T)
m <- seq(1 : length)   # The vector m; will hold the MCMC results
m[1] <- start         # The 1st element in m; the value specified in "start". 
for(i in 2 : length){
m[i] <- sample(1 : states, size = 1, prob = T[m[i - 1],])  
# randomly acquire a new number, 1 through 4, based on                             
# the probabilities in the current row of T.
}
m			  # MCMC results
}
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# Arguments:
# res = A result from MC. 
# states = Number of states. 

Rf <- function(res){
states <- nlevels(as.factor(res))
M <- seq(1, states); M1 <- M  
for(i in 1 : states){ 
M1[i] <- length(res[res == M[i]]) 
}
M1/length(res)			 
}
#------------------------------------------------------------------------#
