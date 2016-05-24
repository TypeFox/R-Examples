## ---- eval=FALSE---------------------------------------------------------
#  install.packages("biotic")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github('robbriers/biotic')

## ------------------------------------------------------------------------
# load library
library(biotic)

# show the format of the built-in almond dataset
head(almond)

## ------------------------------------------------------------------------
# calculate the BMWP index for the River Almond dataset
# 'index' and 'type' do not have to specified as defaults are used
# ("BMWP" and "num")

calcindex(almond)

# calculate the PSI index for the almond samples
# 'type' argument again not needed as the data are numeric abundances
calcindex(almond, "PSI")

## ------------------------------------------------------------------------
# example of processing data in alphabetic log abundance categories
# using the 'type' argument

# 'braidburn' dataset contains alphabetic log category data
# see ?braidburn for details

# calculate the Whalley revised BMWP index (including N-taxa and ASPT)

calcindex(braidburn, "Whalley", "alpha")

# example of processing data in numeric log abundance categories
# using the 'type' argument

# 'greenburn' dataset contains numeric log category data
# see ?greenburn for details

# calculate the LIFE index for this dataset

calcindex(greenburn, "LIFE", "log")

## ------------------------------------------------------------------------
# calculate the BMWP index for almond samples

calcBMWP(almond)

# calculate the AWIC index for almond samples

calcAWIC(almond)

# calculate the WHPT abundance-weighted index for Green Burn samples
# (numeric log abundance categories)

calcWHPT_AB(greenburn, "log")

# calculate LIFE index for Braid Burn samples (alphabetic log categories)

calcLIFE(braidburn, type="alpha")

