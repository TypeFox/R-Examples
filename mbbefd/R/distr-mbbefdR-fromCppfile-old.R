

#exposure curve

.G<-function(x, a, b, g)
{
  stop("deprecated.")
}

#its derivative
dG<-function(x,a,b,g)
{
  stop("deprecated.")
}

#the survival function
.Sx<-function(x,a,b,g)
{
  stop("deprecated.")
}

#the function to compute the exposure function
#could add the following line mbbefdExposure <- ecmbbefdR in mbbefdR-1stparam.R ?
mbbefdExposure<-function(x, a, b, g)
{
 stop("deprecated.")
}



####################################
#classical functions
#distribution function
pmbbefd2<-function(q,a,b,g)
{
  stop("deprecated.")
}

#random generation function (now using the Rcpp version)
#.f4Random<-function(x,a,b) ifelse( ( x>= 1-(a+1)*b/(a+b) ),1,log( (a*(1-x)) / (a+x) ) /log(b))  




#density functio

dmbbefd2<-function(x,a,b,g)
{
  stop("deprecated.")
}


