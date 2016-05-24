#SOME UTILITY FUNCTIONS FOR USE IN THIS PACKAGE

#Function to determine the actual window length from possibly decimal (fractional) window length
.determineWindowLength <- function(x,w){
  if(!is.numeric(x) | !is.numeric(w)){stop("Arguments 'x' and 'w' must be numeric")}
  
  #Vars
  l      = length(x)
  w.orig = w
  
  #If Fractional
  if(w > 0 && w < 1){w = w*l} 
  
  #Positive Integer not less than 1
  ret = as.integer(max(abs(w[1]),1))
  
  #Just in case.
  if(length(x) <= ret | ret < 1){
    stop(paste("Resultant window length is out of range (",ret,"), must be >= 1",sep=""),call.=FALSE)
  }
  
  #Done
  return(ret)
}

#Ifthenelse function
.ifthenelse <- function(t,y,n){if(t){y}else{n}}

#Hidden Normalization Function, sums weights to 1.
.normalize <- function(x){
  if(!is.numeric(x)){stop("argument 'x' must be numeric")}
  
  #If Ony One Value, return unity.
  if(length(x) == 1){return(1)} 
  
  #Else determine the total
  total = sum(x)
  
  #Divide by Zero Error Check
  if(total == 0){stop("argument 'x' must not sum to zero")}
  
  #Done
  return(x/total)
}