## beta distributed radii
lam <-3000
theta <- list("shape1"=2,"shape2"=4)
S <- simSphereSystem(theta,lam, rdist="rbeta", box=list(c(0,5)),pl=101)

sp <- planarSection(S,d=2.5)
ret <- unfold(sp,nclass=20)
 
## Point process intensity
cat("Intensities: ", sum(ret$N_V)/25, "vs.",lam,"\n")
 
## original diameters
r3d <- unlist(lapply(S,function(x) 2.0*x$r))
rest3d <- unlist(lapply(2:(length(ret$breaks)),
            function(i) rep(ret$breaks[i],sum(ret$N_V[i-1]))))
 
op <- par(mfrow = c(1, 2))
hist(r3d[r3d<=max(ret$breaks)], breaks=ret$breaks, main="Radius 3d",
     freq=FALSE, col="gray",xlab="r")
hist(rest3d, breaks=ret$breaks,main="Radius estimated",
     freq=FALSE, col="gray", xlab="r")
par(op)
