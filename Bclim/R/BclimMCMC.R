BclimMCMC <-
function(Bclimdata,chron.loc,nchron=10000,control.mcmc=list(iterations=100000,burnin=20000,thinby=40,report=100),control.chains=list(v.mh.sd=2,phi1.mh.sd=1,phi2.mh.sd=10,vstart=statmod::rinvgauss(Bclimdata$n-1,2,1),Zstart=sample(1:Bclimdata$G,Bclimdata$n,replace=TRUE),phi1start=rep(3,Bclimdata$m),phi2start=rep(20,Bclimdata$m)),control.priors=list(phi1dlmean=rep(1.275,Bclimdata$m),phi1dlsd=rep(0.076,Bclimdata$m),phi2dlmean=rep(4.231,Bclimdata$m),phi2dlsd=rep(0.271,Bclimdata$m))) {

# Create output matrices
remaining <- (control.mcmc$iterations-control.mcmc$burnin)/control.mcmc$thinby
if(remaining!=as.integer(remaining)) stop("Iterations minus burnin divided by thinby must be an integer")

vout <- rep(0,length=Bclimdata$m*(Bclimdata$n-1)*remaining)
zout <- rep(0,length=Bclimdata$m*Bclimdata$n*remaining)
chronout <- rep(0,length=Bclimdata$n*remaining)
cout <- rep(0,length=Bclimdata$m*(Bclimdata$n)*remaining)
phi1out = phi2out = rep(0,length=Bclimdata$m*remaining)

# Re-dim the precisions matrix
Bclimprec <- Bclimdata$tau.mat[,1,]

cat("\n")
cat("Running MCMC...\n")

# Run C code
out <- .C("BclimMCMC3D", 
        as.integer(Bclimdata$G),
        as.integer(Bclimdata$n),
        as.integer(Bclimdata$m),
        as.integer(nchron),
        as.double(Bclimdata$p.mat),
        as.double(Bclimdata$mu.mat),
        as.double(Bclimprec),
        as.character(chron.loc),
        as.integer(control.mcmc$iterations),
        as.integer(control.mcmc$burnin),
        as.integer(control.mcmc$thinby),
        as.integer(control.mcmc$report),
        as.double(control.chains$v.mh.sd),
        as.double(control.chains$vstart),
        as.integer(control.chains$Zstart),
        as.double(control.chains$phi1start),
        as.double(control.chains$phi2start),
        as.double(vout),
        as.integer(zout),
        as.double(chronout),
        as.double(cout),
        as.double(phi1out),
        as.double(phi2out),
        as.double(control.priors$phi1dlmean),
        as.double(control.priors$phi1dlsd),
        as.double(control.priors$phi2dlmean),
        as.double(control.priors$phi2dlsd),
        as.double(control.chains$phi1.mh.sd),
        as.double(control.chains$phi2.mh.sd)
        )
          #,PACKAGE="Bclim")

#browser()
vout  <- array(NA,dim=c(remaining,Bclimdata$n-1,Bclimdata$m))
cout  <- array(NA,dim=c(remaining,Bclimdata$n,Bclimdata$m))
for(i in 1:remaining) {
    for(j in 1:Bclimdata$m) {
        vout[i,,j] <- out[[18]][seq(1,Bclimdata$n-1)+(j-1)*(Bclimdata$n-1)+(i-1)*(Bclimdata$n-1)*Bclimdata$m]
        cout[i,,j] <- out[[21]][seq(1,Bclimdata$n)+(j-1)*(Bclimdata$n)+(i-1)*(Bclimdata$n)*Bclimdata$m]
    }   
}
chronout <- matrix(out[[20]],ncol=Bclimdata$n,nrow=remaining,byrow=TRUE)
zout <- matrix(out[[19]],ncol=Bclimdata$n,nrow=remaining,byrow=TRUE)
phi1out = matrix(out[[22]],ncol=Bclimdata$m,nrow=remaining,byrow=TRUE)
phi2out = matrix(out[[23]],ncol=Bclimdata$m,nrow=remaining,byrow=TRUE)

                           
return(list(v.store=vout,chron.store=chronout,c.store=cout,z.store=zout,zout=zout,chron.loc=chron.loc,nchron=nchron,phi1.store=phi1out,phi2.store=phi2out))

}
