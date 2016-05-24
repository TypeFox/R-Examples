library(ConnMatTools)
data(chile.loco)

# Get appropriate collapse slope
# critical.FLEP=0.2 is just an example
slope <- settlerRecruitSlopeCorrection(chile.loco,critical.FLEP=0.2)

# Make the middle 20 sites a reserve
# All other sites: scorched earth
n <- dim(chile.loco)[2]
LEP <- rep(0,n)
nn <- round(n/2)-9
LEP[nn:(nn+19)] <- 1

Rmax <- 1

recruits0 <- rep(Rmax,n)

# Use DPR model
ret <- DispersalPerRecruitModel(LEP,chile.loco,recruits0,1:20,slope=slope,Rmax=Rmax,
                                settler.recruit.func=BevertonHolt)
image(1:n,1:20,ret$settlers,xlab="sites",ylab="timesteps",
      main=c("Settlement","click to proceed"))
locator(1)
plot(ret$settlers[,20],xlab="sites",ylab="equilibrium settlement",
     main="click to proceed")
locator(1)

# Same, but with a uniform Laplacian dispersal matrix and hockeyStick
cm <- laplacianConnMat(n,10,0,"circular")
ret <- DispersalPerRecruitModel(LEP,cm,recruits0,1:20,slope=1/0.35,Rmax=Rmax)
image(1:n,1:20,ret$settlers,xlab="sites",ylab="timesteps",
      main=c("Settlement","click to proceed"))
locator(1)
plot(ret$settlers[,20],xlab="sites",ylab="equilibrium settlement")
