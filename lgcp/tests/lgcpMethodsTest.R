library(lgcp)
library(spatstat)
library(sp)
library(raster)

set.seed(1)

sd <- lgcpSimSpatial()

xyt <- lgcpSim()

# check plotting

plot(sd)
plot(xyt)

# now check fitting algorithms

sdsave <- file.path(tempdir(),"lgm")
xytsave <- file.path(tempdir(),"lg")

exceed <- exceedProbs(c(1.5,2))

lgm <- lgcpPredictSpatial(  sd=sd,
                            model.parameters=lgcppars(sigma=2,phi=0.1),
        					spatial.covmodel="exponential",
        					cellwidth=0.1,
        					spatial.intensity=density(sd),				
        					mcmc.control=mcmcpars(mala.length=200,burnin=20,
                                retain=20,adaptivescheme=andrieuthomsh(inith=1,alpha=0.5,C=1,
                                targetacceptance=0.574)),
        					output.control=setoutput(gridfunction=
                                dump2dir(dirname=sdsave,forceSave=TRUE),
                                gridmeans=MonteCarloAverage("exceed")),
        					gradtrunc=Inf, # no gradient truncation
        					ext=2)

lg <- lgcpPredict(  xyt=xyt,
                    T=8,
				    laglength=4,
				    model.parameters=lgcppars(sigma=2,phi=0.1,theta=1),
				    spatial.covmodel="exponential",
				    cellwidth=0.1,
				    spatial.intensity=density(xyt),
				    temporal.intensity=function(x){return(100)},					
				    mcmc.control=mcmcpars(mala.length=200,burnin=20,
				        retain=20,adaptivescheme=andrieuthomsh(inith=0.01,alpha=0.5,C=1,targetacceptance=0.574)),
				    output.control=setoutput(gridfunction=
				                dump2dir(dirname=xytsave,forceSave=TRUE),
                                gridmeans=MonteCarloAverage("exceed")),	
				    autorotate=FALSE,
				    gradtrunc=Inf)

# check extraction

meanfield(lgm)
meanfield(lg)

varfield(lgm)
varfield(lg)

rr(lgm)
rr(lg)

serr(lgm)
serr(lg)

intens(lgm)
intens(lg)

seintens(lgm)
seintens(lg)				    
				    
				    
# now check conversion

as.array(lgm$y.mean)
as.array(lg$y.mean)		

raster(meanfield(lgm))
raster(varfield(lg))

as.SpatialPixelsDataFrame(lg$y.mean)
as.SpatialPixelsDataFrame(lgm$y.var)

# check quantiles

quantile(lgm,c(0,0.1,0.5))
quantile(lg,c(0,0.1,0.5))
			   
t1 <- all.equal(expectation(lgm,function(x){return(x)})[[1]],lgm$y.mean$grid[[1]])
if(!t1){stop("error in computing expectation in lgcpMethodsTest.R")}
t2 <- all.equal(expectation(lg,function(x){return(x)})[[1]],lg$y.mean$grid[[length(lg$y.mean$grid)]])
if(!t2){stop("error in computing expectation in lgcpMethodsTest.R")}

				    