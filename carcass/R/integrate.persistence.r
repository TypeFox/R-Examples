# Function to integrate persistence probabilities to account for 
# constant (continuous) arrival of carcasses
#-----------------------------------------------------------------------
integrate.persistence <- function(s, n=20, d=1){
  

if(length(s)>1){                 # integration over the last time intervval
  pprop <- s/c(1,s[-length(s)])
  s1 <- s
  s <- c(1,s1[-length(s1)])*((pprop-1)/log(pprop))
  s[pprop==1] <- s1[pprop==1]       # in cases where survivor function does not decrease or persistence prob is equal
}

if(length(s)==1){
  s1 <- s
  ns <- ifelse(is.na(n*d), 100, n*d)
  s <- s1^c(1:ns)/log(s1) - s1^c(0:(ns-1))/log(s1)
}
return(s)
}

