# This file shows the calibrator package used on a simple 1-d test
# case.  It generates some data randomly, then attempts to infer what
# the parameters used to generate that data are.  One can then compare
# the estimates with the true values.  

# Many of the computationally intensive steps (eg stage2()) are
# executed with only a small number of evaluations, to save time. 


# First, a Boolean specifying whether or not to do the hyperparamter
# optimization (which takes a long time):
do.hyper.opt <- TRUE


#Load the libraries:
library(calibrator)
library(emulator)

# First, how many code observations:
n1 <- 21

# and how many field observations:
n2 <- 30



# Now, source some files that set up the problem:
source("extractor_maker_1d.R")
source("design_maker_1d.R")
source("h1_h2_maker_1d.R")
source("params_maker_1d.R")
source("phi_fun_1d.R")
source("phi_maker_1d.R")
source("Ed_theta_1d.R")
source("data_maker_1d.R")


# Right, now can we estimate the coefficients?  The true values are
# given in params_maker_1d.R, named beta1.TRUE and beta2.TRUE.  The
# thing estimated by betahat.fun.koh is c(beta1,beta2):
betahat.fun.koh(theta=theta.TRUE, d=d.1d, D1=D1.1d, D2=D2.1d, H1=H1.1d, H2=H2.1d, phi=phi.TRUE)

# Ok, not too bad; the true answer is (0,1,1,1,1), so there is an
# error of about 2%


# Now what happens if we use a very wrong  value for theta:
betahat.fun.koh(theta=100, d=d.1d, D1=D1.1d, D2=D2.1d, H1=H1.1d,
                H2=H2.1d, phi=phi.TRUE)

# Terrible!  the estimate is waaaay off!



# Now, can we calculated the likelihood of some thetas?  I define two
# wrapper functions, f(), which gives a log value (ie the support),
# and g() that returns the actual value (ie the likelihood).

f.single <-
  function(x){p.eqn8.supp(theta=x, D1=D1.1d, D2=D2.1d, H1=H1.1d,
                          H2=H2.1d, d=d.1d,
                          phi=phi.TRUE,return.log=TRUE)}
f.single_log <-
  function(x){p.eqn8.supp(theta=x, D1=D1.1d, D2=D2.1d, H1=H1.1d,
                          H2=H2.1d, d=d.1d,
                          phi=phi.TRUE,return.log=FALSE)}
f <- function(x){sapply(x,f.single)}
g <- function(x){sapply(x,f.single_log)}

# Thus f() and g() differ in that f() returns the log likelihood and
# g() returns the likelihood.


# We can plot these, remembering that the true value of theta is
# specified in params_maker_1d.R as theta.TRUE, which is 0.5:

par(ask=TRUE)
x <- seq(from=0,to=1,len=63)
plot(x,g(x),type= "b",
     main="likelihood function for theta; true value is 0.5", 
     xlab=expression(theta)
     )
abline(v=0.5)





# Now use Metropolis-Hastings.  We need a wrapper to give the PDF for theta:


pi.1d <- function(x){p.eqn8.supp(theta=x, D1=D1.1d, D2=D2.1d, H1=H1.1d,
                             H2=H2.1d, d=d.1d,
                             phi=phi.TRUE,return.log=FALSE)}


#So to sample from the posterior for theta we need MH() with a suitable kernel:

theta.sample.1d <- MH(n=10, start=0.5,sigma=diag(1)*0.1 ,pi=pi.1d)

# Is it Gaussian?
if(max(theta.sample.1d)>min(theta.sample.1d)){
  shapiro.test(theta.sample.1d)
} else {
  print("shapiro test not applicable because MH() returned a sample of identical values")
}

# Bear in mind that there is no "true" value for the process at x=0.1.
# The best we can do is to determine the posterior distribution of the
# process at x=0.1 and this is done using Ez.eqn9.supp():

Ez.eqn9.supp(x=0.1,  theta=theta.sample.1d,  d=d.1d, D1=D1.1d,  D2=D2.1d, H1=H1.1d, H2=H2.1d, phi=phi.TRUE)



# We can now evaluate the covariance at, say, x=0.4 for a range of thetas:
jj.theta <- as.matrix( (1:10)/10)
cov.p5.supp(x=0.4,xdash=NULL,theta=jj.theta,d=d.1d, D1=D1.1d, D2=D2.1d, H1=H1.1d, H2=H2.1d , phi=phi.1d)






# Now, all the above was conditional on the hyperparameters being
# correct.  The next few lines will try to estimate the
# hyperparameters.  Bear in mind that this is hard.

# Also bear in mind that the number of iterations is small here, to
# save time.  If you want to see how the optimizations perform with
# more iterations, set variable "maximum.iterations" to a larger value:

maximum.iterations <- 2




# First, try stage1():

if(do.hyper.opt){
phi.stage1 <- stage1(D1=D1.1d, y=y.1d, H1=H1.1d, maxit=maximum.iterations,
 method="SANN", trace=0, do.print=FALSE, phi.fun=phi.fun.1d,
 phi=phi.1d)

# The only thing that this changes is psi1:
phi.stage1$psi1


# (although changing psi1 has a trickle-down effect on some other
# elements of phi:
phi.stage1$sigma1squared
phi.stage1$omega_x
phi.stage1$omega_t



# OK, now stage 2.  Because this is hard, we can extract a subset of
# the observations and use them.  

  use1 <- 1:10
  use2 <- 1:11
  phi.stage2 <- stage2(D1=D1.1d[use1,,drop=FALSE], D2=D2.1d[use2,,drop=FALSE], H1=H1.1d, H2=H2.1d,
      y=y.1d[use1], z=z.1d[use2], extractor=extractor.1d,
     phi.fun=phi.fun.1d, E.theta=E.theta.1d, Edash.theta=Edash.theta.1d,
     maxit=maximum.iterations, method="SANN", phi=phi.stage1)


# This will change psi2, _and_ rho _and_ lambda:
phi.stage2$rho
phi.stage2$lambda
phi.stage2$sigma2squared
phi.stage2$psi2


 # Right, now calibrated prediction.  First, we need to specify a
 # distribution for X:
jj.xdist.mean <- rep(0.5,1)
names(jj.xdist.mean) <- c("x")
jj.xdist.var <- diag(c(0.1),nrow=1)
rownames(jj.xdist.var) <- c("x")
colnames(jj.xdist.var) <- c("x")

X.dist.1d <- list(mean=jj.xdist.mean,var=jj.xdist.var)

# And specify a hbar function:
hbar.fun.1d <- 
function (theta, X.dist, phi) 
{
    if (is.vector(theta)) {
        theta <- t(theta)
    }
    first.bit <- phi$rho * H1.1d(D1.fun(X.dist$mean, theta))
    second.bit <- H2.1d(X.dist$mean)
    jj.names <- colnames(second.bit)
    second.bit <- kronecker(second.bit, rep(1, nrow(first.bit)))
    colnames(second.bit) <- jj.names
    return(t(cbind(first.bit, second.bit)))
}


# And finally we can call EK.eqn10.supp():
EK.eqn10.supp(X.dist=X.dist.1d, D1=D1.1d, D2=D2.1d,
  H1=H1.1d, H2=H2.1d, d=d.1d, hbar.fun=hbar.fun.1d,
  lower.theta=c(-3), upper.theta=c(3),extractor=extractor.1d,
  phi=phi.stage2)

} else {
   print("optimization not performed because variable do.hyper.opt set to FALSE.  To do the optimization, set variable do.hyper.opt to TRUE, near the beginning of this file, and source() it again")
}






