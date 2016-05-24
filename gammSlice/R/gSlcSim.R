gSlcSim <-
function(eg, family, numRep = 10, numGrp = 100, randomFactor)
{
# Set seed

  set.seed(49032)

# Find the number of observations
   
   if (randomFactor == FALSE) {
       numRep = 1
   }   

   num.obs <- numGrp * numRep

# Define true link of mean function

   f3 <- function(x) return(2*x + cos(4*pi*x))
   f4 <- function(x) return(sin(2*pi*x^2))
   f5 <- function(x) return(x + cos(x*pi))

# Create the link function

getobs <- function(family,lipred) {
   if (family=="binomial"){
      mu <- 1/(1+exp(-(lipred)))
      y <- rbinom(num.obs,1,mu)
   }

   if (family=="poisson"){
      mu <- exp(lipred)
      y <- rpois(num.obs,mu)
   }
   return(y)
}

# Set true value of offset parameter
# and random intercept variance component.
   beta0.true <- 0.2
   beta1.true <- 0.5
   beta2.true <- 0.7
   sigsq.u.true <- 1.0


# Create random factor and id number

   Udata<- rep(0,num.obs)
   Unorm<-rnorm(numGrp,0,sqrt(sigsq.u.true))

   id.num <- matrix(0,num.obs,1)
   for (id in 1:numGrp) {
       id.num[(numRep *(id-1)+1):(numRep * id),1] <- id 
       Udata[(numRep *(id-1)+1):(numRep * id)]  <- Unorm[id]
   }
   if (randomFactor == FALSE) {
      Udata <- 0
      id.num <- c(1:num.obs)
   }
 
# Construct the data matrix.
   
   if (eg == 1) {
       x1 <- sample(c(0,1), num.obs, replace = TRUE)
       mu1 <- beta0.true + beta1.true * x1  + Udata
       y <- getobs(family,mu1)
       dataYX <- data.frame(y=y,x1=x1) 
       if (randomFactor) {
          dataYX <- data.frame(idnum = id.num,y=y,x1=x1)
       } else {dataYX <- data.frame(y=y,x1=x1)}
   }

   if (eg == 2) {
       x1 <- sample(c(0,1), num.obs, replace = TRUE)
       x2 <- sample(c(0,1), num.obs, replace = TRUE)
       lipred2 <- beta0.true + beta1.true * x1 + beta2.true * x2 + Udata
       y <- getobs(family,lipred2)
       if (randomFactor) {
           dataYX <- data.frame(idnum=id.num, y=y, x1=x1, x2=x2)
       } else {dataYX <- data.frame(y=y, x1=x1,x2=x2)}   
   }

   if (eg == 3) {
       x1 <- runif(num.obs)
       mu3 <- f3(x1) + Udata
       y <- getobs(family,mu3)
       if (randomFactor) {
           dataYX <- data.frame(idnum = id.num,y=y, x1=x1)
       } else {dataYX <- data.frame(y=y,x1=x1)}
   }

   if (eg == 4) {
      x1 <- rnorm(num.obs)
      x2 <- runif(num.obs)
      lipred4 <- beta1.true * x1 + f4(x2) + Udata
      y <- getobs(family,lipred4)
      if (randomFactor) {
          dataYX <- data.frame(idnum = id.num, y=y, x1=x1, x2=x2)
      } else {dataYX <- data.frame(y=y,x1=x1,x2=x2)} 
   }

   if (eg == 5) {
      x1 <- runif(num.obs,-1,1)
      x2 <- runif(num.obs,-1,1)
      lipred5 <- f3(x1) + f4(x2) + Udata
      y <- getobs(family, lipred5)
      if (randomFactor) {
         dataYX <- data.frame(idnum = id.num, y=y, x1=x1, x2=x2)
      } else {dataYX <- data.frame(y=y,x1=x1,x2=x2)
      } 
   }

   if (eg == 6) {
      x1 <- sample(c(0,1), num.obs, replace = TRUE)
      x2 <- sample(c(0,1), num.obs, replace = TRUE)
      x3 <- runif(num.obs)
      x4 <- runif(num.obs)
      
      lipred6 <- beta0.true + beta1.true * x1 + beta2.true * x2 +
             f3(x3) + f4(x4) + Udata
      y <- getobs(family,lipred6)
      if (randomFactor) {
          dataYX <- data.frame(idnum = id.num, y=y, x1 = x1,
                               x2 = x2,x3 = x3, x4 = x4) 
      } else {dataYX <- data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4)
      }
   }

   return(dataYX)
}
