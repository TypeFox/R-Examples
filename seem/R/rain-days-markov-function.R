# weather generator daily

markov.rain <- function(P, ndays, amount.param,plot.out=T){

# arguments: markov matrix P & number of days
# and amount stats parameters 

mu <- amount.param[[1]]

# define array with all 0 
x <- rep(0,ndays); wet <- 0

# start first day with rain at random
y <- runif(1,0,1)
if(y > 0.5) {x[1] <- rain.day(amount.param); wet<- wet+1}

# loop for remaining days
for(i in 2:ndays) {

# apply markov
  y <- runif(1,0,1)
  if(x[i-1]==0) {
   if(y > P[1,1]) x[i] <- rain.day(amount.param)
  }
  else {
   if(y > P[1,2]) x[i] <- rain.day(amount.param)
  }
 if(x[i] >0) wet <- wet+1

} # end of days loop
expec.wet.days <- ndays*P[2,1]/(P[1,2]+P[2,1])
expec.dry.days <- ndays - ndays*P[2,1]/(P[1,2]+P[2,1])
dry <- ndays-wet
rain.tot<-round(sum(x),2);expec.rain.tot <- expec.wet.days*mu
rain.avg=round(mean(x),2); expec.avg <- expec.rain.tot/ndays
rain.wet.avg=round(sum(x)/wet,2); expec.wet.avg <- mu

if(plot.out==T){
 mat<- matrix(1:2,2,1,byrow=T)
 layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
 plot(x,type="s",xlab="Day",ylab="Rain (mm/day)")
 Rain <- x
 hist(Rain,prob=T,main="",xlab="Rain (mm/day)")
}
return(list(x=round(x,2),wet.days=wet,expected.wet.days=expec.wet.days,
            dry.days=dry,expected.dry.days=expec.dry.days,
            rain.tot=rain.tot, expec.rain.tot=expec.rain.tot,
            rain.avg=rain.avg,expected.rain.avg = expec.avg,
            rain.wet.avg=rain.wet.avg,expected.wet.avg=expec.wet.avg))
}




