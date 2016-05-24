#Loading package
library(R0)

## Woodall reported an attack rate of 0.31 in a population of 1732 during
## the 1957 H2N2 influenza pandemic ('Age and Asian Influenza, 1957', BMJ, 1958)

est.R0.AR(pop.size=1732, AR=0.31)
# Reproduction number estimate using Attack Rate method
# R :  1.19698[ 1.179606 , 1.215077 ]

est.R0.AR(AR=0.31)
# Reproduction number estimate using  Attack Rate  method.
# R :  1.19698

est.R0.AR(pop.size=1732, incid=31)
# Reproduction number estimate using Attack Rate method
# R :  1.009057[ 1.005873 , 1.012269 ]

est.R0.AR(pop.size=1732, incid=c(2,3,4,7,4,2,4,5))
# Reproduction number estimate using Attack Rate method
# R :  1.009057[ 1.005873 , 1.012269 ]

est.R0.AR(pop.size=1732, incid=c(2,3,0,7,4,2,0,5))
# Reproduction number estimate using Attack Rate method
# R :  1.006699[ 1.003965 , 1.009453 ]
