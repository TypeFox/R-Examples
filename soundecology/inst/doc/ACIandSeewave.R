## ----results='show', message=FALSE, warning=FALSE------------------------
library(seewave)
library(soundecology)

data(tropicalsound)
duration <- length(tropicalsound@left)/tropicalsound@samp.rate # 20 seconds

duration

ACI(tropicalsound, nbwindows=(duration/5))

acoustic_complexity(tropicalsound) # j is set to 5 by default


