## ----results='show', message=FALSE, warning=FALSE------------------------
library(soundecology)
data(tropicalsound)

#Using Shannon's Diversity Index:
acoustic_diversity(tropicalsound)

#Using the original code:
acoustic_diversity(tropicalsound, shannon = FALSE)

