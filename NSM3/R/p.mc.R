function (E, n, min.reps=100, max.reps=1000, delta=10^-3) 
{
  # returns the monte carlo estimate for the probability.  
  dsn = numeric() #initialize
  for(i in 1:min.reps){
    dsn = c(dsn,find.Epstein(rexp(n)))
  }
  reps = min.reps
  while(reps<=max.reps){
    p = length(dsn[dsn>E])/reps
    dsn = c(dsn,find.Epstein(rexp(n)))
    if(abs(p-length(dsn[dsn>E])/reps)<=delta){
      return(p) #if p converges to be w/i delta, then return
    }
    reps = reps +1
  }
  print("Warning: reached maximum reps without converging within delta")
  return(p)
  
}
