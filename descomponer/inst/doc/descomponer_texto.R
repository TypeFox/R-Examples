### R code from vignette source 'descomponer_texto.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: descomponer_texto
###################################################
options(keep.source = TRUE, width = 60)
descomponer_texto <- packageDescription("descomponer")


###################################################
### code chunk number 2: descomponer
###################################################
library(descomponer)
data(ipi)
descomponer(ipi,12,1)$datos
gdescomponer(ipi,12,1,2002,1)


