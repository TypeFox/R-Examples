SimulateAR1 <-
function(n, phi){
sde <- sqrt(1-phi^2)
e<-rnorm(n, sd=sde) #variance of time series will be 1
u<-numeric(n)
u[1]<-rnorm(1)
for (i in 2:n)
    u[i] <- phi*u[i-1] + e[i]
u
}

