log_normalize <-
function(vals,maxi = max(vals), ID = which.max(vals)[1]){
### Normalization of a vector of log probabilities.
### K is the length of the vector
### The maximum log-p value and the corresponding element ID (R ID first element has ID 1) are past as arguments
### It returns the corresponding vector of probabilities in the first element
### the log probabilities in the second element
### and the log of the sum of the probabilities in the third element
  ID <- ID -1
  .Call("log_norm",vals=as.numeric(vals),max=as.numeric(maxi),ID=as.integer(ID))

}
