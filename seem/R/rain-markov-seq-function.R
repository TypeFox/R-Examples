# rain sequence daily
markov.rain.seq <- function(P, ndays){

# arguments: markov matrix P & number of days

# define array with all 0 
x <- rep(0,ndays); wet <- 0

# start first day with rain at random
y <- runif(1,0,1)
if(y > 0.5) {x[1] <- 1; wet<- wet+1}

# loop for remaining days
for(i in 2:ndays) {

# apply markov
  y <- runif(1,0,1)
  if(x[i-1]==0) {
   if(y > P[1,1]) x[i] <- 1
  }
  else {
   if(y > P[1,2]) x[i] <- 1
  }
 if(x[i] >0) wet <- wet+1

} # end of days loop
expec.wet.days <- ndays*P[2,1]/(P[1,2]+P[2,1])
expec.dry.days <- ndays - ndays*P[2,1]/(P[1,2]+P[2,1])
dry <- ndays-wet

return(list(x=round(x,2),wet.days=wet,expected.wet.days=expec.wet.days,
            dry.days=dry,expected.dry.days=expec.dry.days))
}




