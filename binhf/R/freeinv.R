`freeinv` <-
function(y,n){

c<-n+1

a<-(n/2)*(1-sign(cos(y))*(1-(sin(y)+(sin(y)-1/sin(y))/n)^2)^.5)

a
}

