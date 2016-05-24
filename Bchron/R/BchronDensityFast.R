BchronDensityFast <-
function(ages,ageSds,calCurves,pathToCalCurves=system.file('data',package='Bchron'),dfs=rep(100,length(ages)),samples=2000,G=30) {

if(length(ages)!=length(ageSds)) stop("ages and 1-sigma errors must be same length")
if(length(ages)!=length(calCurves)) stop("ages and Calibration curves must be same length")
  
# Calibrate ages
x = BchronCalibrate(ages=ages,ageSds=ageSds,calCurves=calCurves,pathToCalCurves=pathToCalCurves,dfs=rep(100,length(ages)))

# Get number of dates
n = length(x)

# Get a huge load of samples from the posteriors here
thetaBig = vector(length=n*samples)
for(i in 1:n) thetaBig[((i-1)*samples+1):(i*samples)] = sample(x[[i]]$ageGrid,size=samples,prob=x[[i]]$densities,replace=TRUE)

# Now run mclust
mclustOutput = mclust::densityMclust(data=thetaBig,G=G)

output = list(out=mclustOutput,calAges=x)
class(output) = 'BchronDensityRunFast'
return(output)
  
}
