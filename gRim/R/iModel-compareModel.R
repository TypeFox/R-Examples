
.compareModel <- function(x, basestat=c(0,0), crit="aic", k=2,direction="drop", ...){

  switch(crit,
    "aic"={
      ans <- c(-x[1]+basestat[1] + k*(x[2]-basestat[2]))
      if (direction=="add")
        ans <- -ans
    },
    "test"={
      ans <- 1-pchisq(abs(x[1]-basestat[1]), df=abs(x[2]-basestat[2]))
    } 
  )
  return(ans)
}
