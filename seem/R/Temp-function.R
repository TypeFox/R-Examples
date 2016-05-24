# parabolic function
parab <- function(Temp,a,b){
 nT <- length(Temp)
 y <- 4*(b-Temp)*(Temp-a)/(b-a)^2	
 for(i in 1:nT) if(y[i] <0) y[i] <- 0
 return(y)
}

# rate q10
kT.rate <- function(Temp,k0,q10){
   kT <- k0*q10^(Temp/10)
   return(kT)
} 

# Q10 function

Q10 <- function(Temp, param=c(Topt,Tmax,q10)){ 

Topt <- param[1];Tmax<-param[2];q10<- param[3]
 d <- (Tmax-Topt)*(log(q10)/10)		 
 a <- (Tmax - Temp)/(Tmax-Topt)
 y <- a^d*exp(d-a*d)


return(y)
}

