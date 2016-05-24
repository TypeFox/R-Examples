# Table 3.4-1
with(SticklebackPlates,
	data.frame( genotype = levels(genotype),
	            frequency = as.vector(table(genotype)),
				proportion = as.vector(table(genotype)) / nrow(SticklebackPlates)
		)
	)

# Proportion is a mean
mean( SticklebackPlates$genotype == 'mm' )
mean( SticklebackPlates$genotype == 'Mm' )
mean( SticklebackPlates$genotype == 'MM' )
