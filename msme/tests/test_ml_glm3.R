library(msme)

# library(msme, lib.loc="lib")

data(medpar)

med.nb.g <- ml_glm3(los ~ hmo + white,
                   family = "gNegBinomial", 
                   link = "log",
                   group = medpar$provnum, 
                   data = medpar)

summary(med.nb.g)

