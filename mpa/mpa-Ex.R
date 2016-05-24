pkgname <- "mpa"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mpa')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("diagram.mpa")
### * diagram.mpa

flush(stderr()); flush(stdout())

### Name: diagram.mpa
### Title: Strategic diagram
### Aliases: diagram.mpa
### Keywords: multivariate

### ** Examples

data(revista)
mat <- matriz.mpa(revista, fmin=3, cmin=1)
clas <- mpa(mat$Matriza,10,mat$Palabras)
diagram.mpa(clas,tmin=3)



cleanEx()
nameEx("leer.mpa")
### * leer.mpa

flush(stderr()); flush(stdout())

### Name: leer.mpa
### Title: Reading corpus
### Aliases: leer.mpa
### Keywords: multivariate

### ** Examples

#revista <- leer.mpa("revista,txt",encoding="latin1")
data(revista)
revista



cleanEx()
nameEx("matriz.mpa")
### * matriz.mpa

flush(stderr()); flush(stdout())

### Name: matriz.mpa
### Title: Calculation of co-occurrences matrix and matrix associations
### Aliases: matriz.mpa
### Keywords: multivariate

### ** Examples

data(revista)
mat <- matriz.mpa(revista, fmin=3, cmin=1)
mat$Matriza
mat$Matrizc
diag(mat$Matrizc)



cleanEx()
nameEx("mpa")
### * mpa

flush(stderr()); flush(stdout())

### Name: mpa
### Title: CoWords method
### Aliases: mpa contar.si reemplazar.si
### Keywords: multivariate

### ** Examples

#revista <- leer.mpa("revista.txt",encoding="latin1")
data(revista)
mat <- matriz.mpa(revista, fmin=3, cmin=1)
clas <- mpa(mat$Matriza,10,mat$Palabras)
clas



cleanEx()
nameEx("plotmpa")
### * plotmpa

flush(stderr()); flush(stdout())

### Name: plotmpa
### Title: Network group's internal associations
### Aliases: plotmpa
### Keywords: multivariate

### ** Examples

data(revista)
mat <- matriz.mpa(revista, fmin=3, cmin=1)
clas <- mpa(mat$Matriza,10,mat$Palabras)
clas
plotmpa(1,mat$Matriza,clas)
plotmpa(6,mat$Matriza,clas)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
