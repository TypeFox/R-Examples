ST_Unordered <-
function(N,x,p_x)  {
  F_x <- cumsum(p_x)
  disc_sim <- numeric(length=N)
  for(i in 1:N){
    temp <- runif(1)
    disc_sim[i] <- x[min(which(F_x>temp))]
  }
  return(disc_sim)
}
