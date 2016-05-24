GreatTitMalaria
str(GreatTitMalaria)
xtabs(~ response + treatment, GreatTitMalaria)
malariaT <- as.data.frame(xtabs(~ response + treatment, GreatTitMalaria))
malariaT
barchart(Freq ~ treatment, groups=response, malariaT, col=c('red','orange'))
barchart(Freq ~ response | treatment, malariaT, col=c('red','orange'))
barchart(Freq ~ treatment, groups=response, malariaT, stack=TRUE, col=c('red','orange'))

if (require(vcd)) {
	mosaic(~treatment + response, GreatTitMalaria, direction = 'v')
	mosaic(~treatment + response, GreatTitMalaria, shade = TRUE,
			gp=gpar(fill=c('red','orange')),
			labeling=labeling_values,
			direction='v')
	# reorder the response factor 
	GreatTitMalaria$response <- ordered(GreatTitMalaria$response, 
                                      levels=rev(levels(GreatTitMalaria$response)))
	mosaic(~treatment + response, GreatTitMalaria,shade=TRUE,
			gp=gpar(fill=c('red','orange')),
			labeling=labeling_values,
			direction='v')
} 
