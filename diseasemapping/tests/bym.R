havePackages = c(
  'INLA' = requireNamespace('INLA', quietly=TRUE),
  "spdep" = requireNamespace('spdep', quietly=TRUE)
)

print(havePackages)

library('diseasemapping')
data('kentucky')

if(FALSE) {
	# must have an internet connection to do the following
	larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
	# get rid of under 10's
	larynxRates = larynxRates[-grep("_(0|5)$",names(larynxRates))]
	dput(larynxRates)
} else {
	larynxRates = structure(c(0, 0, 0, 0, 1e-06, 6e-06, 2.3e-05, 4.5e-05, 9.9e-05, 
					0.000163, 0.000243, 0.000299, 0.000343, 0.000308, 0.000291, 0.000217, 
					0, 0, 0, 1e-06, 1e-06, 3e-06, 8e-06, 1.3e-05, 2.3e-05, 3.5e-05, 
					5.8e-05, 6.8e-05, 7.5e-05, 5.5e-05, 4.1e-05, 3e-05), .Names = c("M_10", 
					"M_15", "M_20", "M_25", "M_30", "M_35", "M_40", "M_45", "M_50", 
					"M_55", "M_60", "M_65", "M_70", "M_75", "M_80", "M_85", "F_10", 
					"F_15", "F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", 
					"F_55", "F_60", "F_65", "F_70", "F_75", "F_80", "F_85"))
	
}

kentucky = getSMR(kentucky, larynxRates, larynx,
		regionCode="County")

if(all(havePackages)){

  kBYM = bym(
			formula = observed ~ offset(logExpected) + poverty,
      data=kentucky@data,
			adjMat = spdep::poly2nb(kentucky, row.names=kentucky$County),
      priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)),
			region.id="County"
  )
	kBYM$parameters$summary

	pdf("priorPostKentucky.pdf")
	plot(kBYM$parameters$sdSpatial$posterior, type='l', 
			xlim=c(0,max(kBYM$parameters$sdSpatial$priorCI)))
	lines(kBYM$parameters$sdSpatial$prior, col='blue')
	legend('topright', lty=1, col=c('black','blue'), legend=c('posterior','prior'))
	dev.off()

	
	


# also try no covariate or prior

kBYM = bym(
		formula = observed ~ offset(logExpected),
		data=kentucky)


if(require('geostatsp', quietly=TRUE)) {
 	kBYM$data$exc1 = geostatsp::excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
} else {
	kBYM$data$exc1 = rep(NA, length(kBYM$data))
}

kBYM$par$summary

if(require('mapmisc', quietly=TRUE)) {

colFit = colourScale(kBYM$data$fitted.exp,
		breaks=6, dec=3)
	
plot(kBYM$data, col=colFit$plot)
legendBreaks('topleft', colFit)

colExc = colourScale(kBYM$data$exc1 ,
		style='fixed',
		breaks=c(0, 0.2, 0.8,0.9, 1), 
		col=rainbow, rev=TRUE
	)

	plot(kBYM$data, col=colExc$plot)
	legendBreaks('topleft', colExc)
 		
}
# and try passing a data frame and adjacency matrix

	
adjMat = spdep::poly2nb(kentucky, row.names =as.character(kentucky$County) )
kBYM = bym(data=kentucky@data, formula=observed ~ offset(logExpected) + poverty,
		adjMat = adjMat, region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

kBYM$par$summary

# subtract a few regions

kBYM = bym(
    formula=observed ~ offset(logExpected) + poverty,
    data=kentucky@data[-(1:4),],  
 	  adjMat = adjMat, region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))
 

kBYM$par$summary

# intercept only, no offset


kBYM = bym(data=kentucky,  formula=observed ~ 1,
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

kBYM$par$summary


if(require('mapmisc', quietly=TRUE)) {
	
	colFit = colourScale(kBYM$data$fitted.exp,
			breaks=6, dec=1)
	
	plot(kBYM$data, col=colFit$plot)
	legendBreaks('topleft', colFit)
	
}

 

# give spdf but some regions have no data
# but keep the 'county' column as is
kentucky@data[1:2,-grep("County", names(kentucky))] = NA 

kBYM = bym(observed ~ offset(logExpected) + poverty,
		kentucky, 
		region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)))

 
kBYM$par$summary


# missing value in a categorical variable

pCuts = quantile(kentucky$poverty, na.rm=TRUE)
kentucky$povertyFac = cut(kentucky$poverty, 
		breaks = pCuts,
		labels = letters[seq(1,length(pCuts)-1)])
kentucky$povertyFac[c(2,34,100)] = NA

kBYM = bym(
		formula = observed ~ offset(logExpected) + povertyFac,
		data = kentucky, 
		region.id="County",
		priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5))
)


kBYM$par$summary
}


