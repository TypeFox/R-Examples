library("diseasemapping")

larynxRates= try(cancerRates("USA", year=1998:2002,site="Larynx"), silent=TRUE)

if(class(larynxRates)=='try-error') {
	warning("iarc web site appears to be down")
} else {

	larynxRates

data("kentucky")

kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

if(require('mapmisc', quietly=TRUE)) {
	mycol = colourScale(kentucky2$SMR, breaks=9, 
			dec=-log10(0.5), style='equal', transform='sqrt')
	plot(kentucky2, col=mycol$plot)
	legendBreaks('topleft', mycol)
}
}