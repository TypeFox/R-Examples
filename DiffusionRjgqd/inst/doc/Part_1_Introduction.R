## ----eval=FALSE----------------------------------------------------------
#  G1 <- function(t){-sin(2*pi*t)}
#  Q0 <- function(t){1}

## ----eval=FALSE----------------------------------------------------------
#  G1   <- function(t){-sin(2*pi*t)}
#  Q0   <- function(t){1}
#  Lam1 <- function(t){0.15}
#  Jmu  <- function(t){0.0}
#  Jsig <- function(t){0.5}

## ----fig.align = 'center'------------------------------------------------
mu      <- 0.1
sigma.x <- 0.5
lam     <- 0.5
mu.z    <- 0.1
sigma.z <- 0.5

true.density=function(X0,Xt,t,order =10)
{
   dens <- 1/sqrt(2*pi*sigma.x^2*t)*exp(-(Xt-X0)^2/(2*sigma.x^2*t)+mu/sigma.x^2*(Xt-X0))*exp(-(mu^2/(2*sigma.x^2)+lam)*t)
   for(j in 1:order)
   {
      dens <-  dens +exp(-lam*t)*lam^j/(factorial(j)*sqrt(2*pi)*sqrt(sigma.x^2*t+j*sigma.z^2))*exp(-(Xt-X0-mu*t-j*mu.z)^2/(2*(sigma.x^2*t+j*sigma.z^2)))*t^j
   }
   return(list(density=dens,Xt=Xt))
}
Xt  <- seq(-2,2,1/10)
res <- true.density(0,Xt,1)
plot(res$density~Xt,type='l',main='Transition Density',xlab='X_t',ylab='density')

## ----fig.align = 'center'------------------------------------------------
library(DiffusionRjgqd)
# Remove any existing coefficients:
JGQD.remove()

# Define the model coefficients:
G0   <- function(t){mu}
Q0   <- function(t){sigma.x^2}
Lam0 <- function(t){lam}
Jmu  <- function(t){mu.z}
Jsig <- function(t){sigma.z}

# Calculate the transitional density:
BM <- JGQD.density(0,Xt,0,1,factorize=TRUE,Dtype='Normal.A')

# Plot the transitional density:
plot(res$density~Xt,type='l',main='Transition Density',xlab='X_t',ylab='density')
lines(BM$density[,100]~BM$Xt,col='blue',lty='dashed',lwd=2)

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRjgqd')

