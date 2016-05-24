########
# 1D dynamic densities
########

# Simulate the data: 100 time points with 10 obsv each
nx <- 10
total <- nx*100
x <- c()
times <- c()
sd <- 13
xx <- seq(-120,120,length=100)
dd <- c()
for(i in 1:10)
 { r <- rbinom(nx, 1, 0.5)  
   x <- c(x, rnorm(nx, 80, sd)*r + rnorm(nx, -80, sd)*(1-r) )
   times <- c(times, rep(i-1,nx))
   dd <- rbind(dd, dnorm(xx, 80, sd)/2 + dnorm(xx, -80, sd)/2) }
for(i in 1:40)
 { r <- rbinom(nx, 1, 0.5)  
   x <- c(x, rnorm(nx, 80-2*i, sd+i/4)*r + rnorm(nx, -80+2*i, sd+i/4)*(1-r) )
   times <- c(times, rep(10+i-1,nx))
   dd <- rbind(dd, dnorm(xx, 80-2*i, sd+i/4)/2 + dnorm(xx, -80+2*i, sd+i/4)/2) }
for(i in 1:40)
 { r <- rbinom(nx, 1, 0.5)  
   x <- c(x, rnorm(nx, 2*i, sd+(40-i)/4)*r + rnorm(nx, -2*i, sd+(40-i)/4)*(1-r) )
   times <- c(times, rep(50+i-1,nx))
   dd <- rbind(dd, dnorm(xx, 2*i, sd+(40-i)/4)/2 + dnorm(xx, -2*i, sd+(40-i)/4)/2) }
for(i in 1:10)
 { r <- rbinom(nx, 1, 0.5)  
   x <- c(x, rnorm(nx, 80, sd)*r + rnorm(nx, -80, sd)*(1-r) )
   times <- c(times, rep(90+i-1,nx))
   dd <- rbind(dd, dnorm(xx, 80, sd)/2 + dnorm(xx, -80, sd)/2) }


alpha <- 4
params <- c(0, #gamma
            .2, #kappa
            3, #nu
            3, #gam0
            50 #psi0
            )
N <- 50 # very small number of particles! You'll notice markov error in repeated runs

# independent DP for each time
l0 <- mix(x, alpha=alpha, g0params=params,
          times=times, rho=0, cat=0,
          N=N, niter=0, read=0, print=1)

# BAR stick-breaking with rho=1/2
l1 <- mix(x, alpha=alpha, g0params=params,
          times=times, rho=0.5, cat=0,
          N=N, niter=0, read=0, print=1) 

# Plot the Bayes factor for rho=.5 vs independence
bf <- l1$logprob-l0$logprob
par(mai=c(.7,.7,0.4,0.4), mfrow=c(1,1))
plot(c(-100:(total+100)), rep(0,total+201), type="l", col=grey(.5), xlim=c(10,total+10), ylim=range(bf), 
		xlab="", ylab="", main="", cex.axis=.8)
mtext("Log Bayes Factor", side=2, font=3, cex=1.1, line=2.3)
lines(bf, col=6)
text(x=total+20, y=bf[total], label="0.5", cex=.8, font=3) 
mtext("Observation", side=1, font=3, cex=1.1, line=-1.25, outer=TRUE)

# Extract mean pdfs and compare the filtered densities
dens <- function(prt)
  { pdf <- rep(0,100)
    for(j in 1:nrow(prt))
      { pdf <- pdf + prt$p[j]*dt( (xx-prt[j,]$a.1)/sqrt(prt[j,]$B.1), df = prt$c[j] )/sqrt( prt[j,]$B.1 ) }
    return(pdf) }

prts1 <- vector(mode="list", length=0)
prts0 <- vector(mode="list", length=0)
for(t in 1:99){
  prt <-  vector(mode="list", length=N) 
  for(i in 1:N) prt[[i]] <- particle(i, l0, t, 0)
  prts0 <- cbind(prts0, prt)
  for(i in 1:N) prt[[i]] <- particle(i, l1, t, 0.5)
  prts1 <- cbind(prts1, prt) }

post0 <- lapply(prts0,dens)
post1 <- lapply(prts1,dens)

pdfs0 <- array( unlist(post0), dim=c(100,N,99) )
pdfs1 <- array( unlist(post1), dim=c(100,N,99) )

mf0 <- apply(pdfs0, c(1,3), mean)
mf1 <- apply(pdfs1, c(1,3), mean)

rl <- readline("press RETURN to continue: ")

# plot
cols <- rainbow(99)
par(mfrow=c(1,3))
pmat <- persp(x=xx, y=1:100, z=t(dd), theta=20, phi=40, expand=.6, ticktype="detailed", r=100, tcl=.1,
              xlab="x", ylab="time", zlab="",  border=0, col=0,  zlim=range(dd))
text(trans3d(x=-115, y=0, z=.025, pmat=pmat), label="f(x)", cex=1, font=3)
mtext("Filtered AR Fit", side=3, font=3)
for(i in 99:1){ lines(trans3d(x=xx, y=i, z=mf1[,i], pmat=pmat), col=cols[i]) }
pmat <- persp(x=xx, y=1:100, z=t(dd), theta=20, phi=40, expand=.6, ticktype="detailed", r=100,
              xlab="x", ylab="time", zlab="",  border=NA, col=matrix(rep(cols,99), ncol=99, byrow=TRUE), zlim=range(dd) )
text(trans3d(x=-115, y=0, z=.025, pmat=pmat), label="f(x)", cex=1, font=3)
mtext("The Truth", side=3, font=3)
pmat <- persp(x=xx, y=1:100, z=t(dd), theta=20, phi=40, expand=.6, ticktype="detailed", r=100,
              xlab="x", ylab="time", zlab="", border=0, col=0, zlim=range(dd))
text(trans3d(x=-115, y=0, z=.025, pmat=pmat), label="f(x)",  font=3)
for(i in 99:1){ lines(trans3d(x=xx, y=i, z=mf0[,i], pmat=pmat), col=cols[i]) }
mtext("Independent Fit", side=3, font=3)

