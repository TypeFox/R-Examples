# weather generator - daily amount

rain.day <- function(param){
 # mu, std and skew of daily rain used for skew
 # shape for gamma, exp, and Weibull
 # model.pdf "s" skewed, "w" weibull and "e" exponential, "g" gamma
  mu <- param[[1]]; std <-param[[2]]; skew <- param[[3]]
  shape <- param[[4]]; model.pdf <- param[[5]]

 # calc rain
 if(model.pdf=="e") {
   u <- runif(1,0,1) # generate uniform
   scale <- mu	
   y <- scale*(-log(u))
 }
 if(model.pdf=="w") {
   u <- runif(1,0,1) # generate uniform
   scale <- mu/gamma(shape+1)	
   y <- scale*(-log(u))^shape
 }
 if(model.pdf=="g") {
   scale <- mu/shape	
   y <- rgamma(1,scale,shape)
 }
 if (model.pdf == "s") {
   z <- rnorm(1,0,1) # generate standard normal
   y <- mu+ 2*(std/skew)*( ( ((skew/6)*(z-skew/6)+1))^3 -1)
 } 

 return(y)
}

