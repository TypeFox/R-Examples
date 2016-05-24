havePackages = c(
    'INLA' = requireNamespace('INLA', quietly=TRUE)
)

print(havePackages)

# number of cells... smaller is faster but less interesting
Ncell = 25

# as in example
require('geostatsp')
 
data('swissRain')
swissRain$lograin = log(swissRain$rain)

if(all(havePackages)) {
  swissFit =  glgm("lograin", swissRain, Ncell, 
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE)
)

swissFit$parameters$summary
if(!interactive()) pdf("swissGlgmExc.pdf")
swissExc = excProb(swissFit$inla$marginals.random$space, 0, template=swissFit$raster)
plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		
if(!interactive()) dev.off()
if(!interactive()) pdf("swissGlgmExc2.pdf")
swissExcP = excProb(swissFit$inla$marginals.predict, 3, template=swissFit$raster)
plot(swissExcP, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		
if(!interactive()) dev.off()

# intercept only, add a covariate just to confuse it
swissFit =  glgm(formula=lograin~1, 
		data=swissRain, grid=Ncell, 
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
)

swissFit$parameters$summary

if(!interactive()) pdf("swissGlgmExc3.pdf")
	swissExc = excProb(swissFit$inla$marginals.random$space, 0, template=swissFit$raster)
	plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
	plot(swissBorder, add=TRUE)		
if(!interactive()) dev.off()


# now with formula
swissFit =  glgm(lograin~ CHE_alt,
		swissRain, 
		Ncell, 
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFit$parameters$summary



# covariates are in data
newdat = swissRain
newdat$elev = extract(swissAltitude, swissRain)
swissFit =  glgm(lograin~ elev + land,
		newdat, Ncell, 
		covariates=list(land=swissLandType),
		family="gaussian", buffer=40000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFit$parameters$summary
if(!interactive()) pdf("swissCovInData.pdf")
plot(swissFit$raster[['predict.mean']])
if(!interactive()) dev.off()
# formula, named list elements
swissFit =  glgm(lograin~ elev,
		swissRain, Ncell, 
		covariates=list(elev=swissAltitude), 
		family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFit$parameters$summary

# categorical covariates
swissFit =  glgm(
		formula = lograin ~ elev + factor(land),
		data = swissRain, grid = Ncell, 
covariates=list(elev=swissAltitude,land=swissLandType), 
family="gaussian", buffer=20000,
priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
control.family=list(hyper=list(prec=list(prior="loggamma", 
						param=c(.1, .1))))
)
swissFit$parameters$summary
table(swissFit$inla$.args$data$land)
if(!interactive()) pdf("swissFitCategorical.pdf")
plot(swissFit$raster[['predict.mean']])
if(!interactive()) dev.off()
# put some missing values in covaritates
# also don't put factor() in formula
temp = values(swissAltitude)
temp[seq(10000,12000)] = NA
values(swissAltitude) = temp
swissFitMissing =  glgm(rain ~ elev + land,swissRain,  Ncell, 
		covariates=list(elev=swissAltitude,land=swissLandType), 
		family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFitMissing$parameters$summary

}

# these tests are time consuming, so only patrick will do them
if(all(havePackages) & Sys.info()['user'] =='patrick') {

data('loaloa')
rcl = rbind(
		# wedlands and mixed forests to forest
		c(5,2),c(11,2),
# savannas to woody savannas
		c(9,8),
		# croplands and urban changed to crop/natural mosaid
		c(12,14),c(13,14))
ltLoaR = reclassify(ltLoa, rcl)
levels(ltLoaR) = levels(ltLoa)

 
elevationLoa = elevationLoa - 750
elevLow = reclassify(elevationLoa, c(0, Inf, 0))
elevHigh = reclassify(elevationLoa, c(-Inf, 0, 0))

 covList = list(elLow = elevLow, elHigh = elevHigh, 
		land = ltLoaR, evi=eviLoa)

 loaFit = glgm(
		 y ~ land + evi + elHigh + elLow, #+ f(villageID,model="iid"),
		 loaloa,
		  Ncell, 
		  covariates=covList, 
		family="binomial", Ntrials = loaloa$N,
		shape=2, buffer=25000,
		priorCI = list(sd=c(0.2, 4), range=c(20000,500000)))

loaFit$par$summary

if(!interactive()) png("loaFitted.png")
plot(loaFit$raster[['predict.exp']])
if(!interactive()) dev.off()

# prior for observation standard deviation
swissFit =  glgm( formula="lograin",data=swissRain, grid=Ncell,
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.1, 2), range=c(50000,500000), 
				sdNugget=c(0.1, 2)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE)
)



# a model with little data, posterior should be same as prior

data2 = SpatialPointsDataFrame(cbind(c(1,0), c(0,1)),
		data=data.frame(y=c(0,0), offset=c(-50,-50), x=c(-1,1)))

res = glgm(data=data2, grid=20, formula=y~1 + x+offset(offset), 
    covariates=NULL,
		priorCI = list(sd=c(0.3,0.5), range=c(0.25, 0.4)),
		family="poisson",buffer=0.5, 
    control.fixed=list(mean.intercept=0, prec.intercept=1,
        mean=0,prec=4),
    control.mode=list(theta=c(2, 2),restart=TRUE)
)

if(!interactive()) pdf("nodata.pdf")

if(!interactive()) par(mfrow=c(3,1))

# intercept
plot(res$inla$marginals.fixed[['(Intercept)']], col='blue', type='l',
    xlab='intercept',lwd=3)
xseq = res$inla$marginals.fixed[['(Intercept)']][,'x']
lines(xseq, dnorm(xseq, 0, 1),col='red',lty=2,lwd=3)
legend("topright", col=c("blue","red"),lty=1,legend=c("prior","post'r"))

# beta
plot(res$inla$marginals.fixed[['x']], col='blue', type='l',
    xlab='beta',lwd=3)
xseq = res$inla$marginals.fixed[['x']][,'x']
lines(xseq, dnorm(xseq, 0, 1/2),col='red',lty=2,lwd=3)
legend("topright", col=c("blue","red"),lty=1,legend=c("prior","post'r"))


# sd
plot(res$parameters$sd$prior,type='l', col='blue',xlab='sd',lwd=3)
lines(res$parameters$sd$post,col='red',lty=2,lwd=3)
legend("topright", col=c("blue","red"),lty=1,legend=c("prior","post'r"))


# range
plot(res$parameters$range$prior,type='l', col='blue', xlab='range',lwd=3)
lines(res$parameters$range$post,col='red',lty=2,lwd=3)
legend("topright", col=c("blue","red"),lty=1,legend=c("prior","post'r"))
if(!interactive()) dev.off()


# covariates are in data, interactions
newdat = swissRain
newdat$elev = extract(swissAltitude, swissRain)
swissFit =  glgm(
    formula = lograin~ elev : land,
    data=newdat, 
    grid=squareRaster(swissRain,50), 
    covariates=list(land=swissLandType),
    family="gaussian", buffer=0,
    priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
    control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
    control.family=list(hyper=list(prec=list(prior="loggamma", 
                param=c(.1, .1))))
)
swissFit$parameters$summary
if(!interactive()) pdf("swissCovInDataInteraction.pdf")
plot(swissFit$raster[['predict.mean']])
if(!interactive()) dev.off()



}
