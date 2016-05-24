BchronDensity <-
function(ages,ageSds,calCurves,pathToCalCurves=system.file('data',package='Bchron'),dfs=rep(100,length(ages)),numMix=30,iterations=10000,burn=2000,thin=8,updateAges=FALSE) {

if(length(ages)!=length(ageSds)) stop("ages and 1-sigma errors must be same length")
if(length(ages)!=length(calCurves)) stop("ages and Calibration curves must be same length")
  
# Calibrate ages
x = BchronCalibrate(ages=ages,ageSds=ageSds,calCurves=calCurves,pathToCalCurves=pathToCalCurves,eps=0,dfs=rep(100,length(ages)))
xSmall = BchronCalibrate(ages=ages,ageSds=ageSds,calCurves=calCurves,pathToCalCurves=pathToCalCurves,dfs=rep(100,length(ages)))

# Get thetaRange to calculate mu values
n = length(x)
thetaRange = range(xSmall[[1]]$ageGrid)
for(i in 2:n) thetaRange = range(c(thetaRange,xSmall[[i]]$ageGrid))

# Put in offset for normal calibration curve (enables faster lookup)
offset=vector(length=n)
for(i in 1:n) {
  offset[i] = ifelse(x[[i]]$calCurve == 'normal',61,0)
}

# Create some Gaussian basis functions to use
gauss <- function(x, mu, sig) {
  # Gaussian-shaped function
  u <- (x - mu) / sig
  y <- exp(- u * u / 2)
  y }

gbase <- function(x, mus) {
  # Construct Gaussian basis
  sig <- (mus[2] - mus[1]) / 2
  G <- outer(x, mus, gauss, sig)
  G }

clrInv = function(phi) {
  return(exp(phi)/sum(exp(phi)))  
}

# Starting values
J = numMix
mu = seq(thetaRange[1],thetaRange[2],length=numMix)
theta = vector(length=n)
for(j in 1:n) theta[j] = round(stats::rnorm(1,mean=x[[j]]$ageGrid[match(max(x[[j]]$densities),x[[j]]$densities)],sd=ageSds[j]),3)
phi = c(stats::runif(J-1,-10,10),0)
p = as.numeric(clrInv(phi))
G = gbase(theta,mu)
  
# Storage
remaining = (iterations-burn)/thin
thetaStore = matrix(ncol=length(theta),nrow=remaining)
pStore = matrix(ncol=J,nrow=remaining)

# Get thetas up front for faster results
thetaAll = matrix(NA,ncol=n,nrow=iterations)
for(j in 1:n) thetaAll[,j] = sample(xSmall[[j]]$ageGrid,size=iterations,prob=xSmall[[j]]$densities,replace=TRUE)

# Create function for quick calling of mixture density
mu2 = mu
sigma2 = (mu[2] - mu[1]) / 2
my_dnorm = function(x) stats::dnorm(x,mean=mu2,sd=sigma2)

# Loop through iterations
pb = utils::txtProgressBar(min = 1, max = iterations, style = 3,width=60,title='Running BchronDensity')
for(i in 1:iterations) {
  utils::setTxtProgressBar(pb, i)

  # Store stuff
  if(i>burn & i%%thin==0) {
    ind = (i-burn)/thin
    thetaStore[ind,] = theta
    pStore[ind,] = p
  }
  
  # Update theta
  if(updateAges) {
    for(j in 1:n) {
      thetaNew = round(stats::rnorm(1,theta[j],0.5),3)
      thetaNewMatch = as.integer(thetaNew+offset[j])+1      
      thetaNewLogDens = max(log(x[[j]]$densities[thetaNewMatch]),-1000000)
      priorNew.dens = sum(p*stats::dnorm(thetaNew,mean=mu2,sd=sigma2))
      thetaMatch = as.integer(theta[j]+offset[j])+1
      thetaLogDens = max(log(x[[j]]$densities[thetaMatch]),-1000000)
      priorDens = sum(p*stats::dnorm(theta[j],mean=mu2,sd=sigma2))
        
      logRtheta = thetaNewLogDens - thetaLogDens + log(priorNew.dens) - log(priorDens)
      if(stats::runif(1)<exp(logRtheta)) theta[j] = thetaNew      
    }    
  } else {
    theta = thetaAll[i,]  
  }
  
  # Update phi
  for(j in 1:(J-1)) {
    phiNew = stats::rnorm(1,phi[j],1)
    phiAllNew = phi
    phiAllNew[j] = phiNew
    pNew = as.numeric(clrInv(phiAllNew))
    phiNewLogDens = sum(log(G%*%pNew))
    phiLogDens = sum(log(G%*%p))
    logRphi = phiNewLogDens - phiLogDens + stats::dunif(phiNew,-10,10,log=TRUE) - stats::dunif(phi[j],-10,10,log=TRUE)
    if(stats::runif(1)<exp(logRphi)) {
      phi[j] = phiNew
      p = as.numeric(clrInv(phi))
    } 
  }
}

output = list(theta=thetaStore,p=pStore,mu=mu,calAges=xSmall,G=G)
class(output) = 'BchronDensityRun'
return(output)
  
}
