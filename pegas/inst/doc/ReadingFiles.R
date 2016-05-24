### R code from vignette source 'ReadingFiles.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ReadingFiles.Rnw:31-32
###################################################
options(width=60)


###################################################
### code chunk number 2: ReadingFiles.Rnw:106-111
###################################################
library(pegas)
x <- read.loci("toto", header = FALSE)
x
print(x, details = TRUE)
class(x)


###################################################
### code chunk number 3: ReadingFiles.Rnw:129-132
###################################################
y <- read.loci("titi", header = FALSE, allele.sep = "-")
print(y, details = TRUE)
identical(x, y)


###################################################
### code chunk number 4: ReadingFiles.Rnw:137-138
###################################################
args(read.loci)


###################################################
### code chunk number 5: ReadingFiles.Rnw:165-166
###################################################
print(read.loci("tutu", FALSE), TRUE)


###################################################
### code chunk number 6: ReadingFiles.Rnw:182-185
###################################################
X <- read.loci("tyty")
print(X, TRUE)
summary(X)


###################################################
### code chunk number 7: ReadingFiles.Rnw:203-206
###################################################
z <- read.loci("tata", loci.sep = "\t", col.loci = 2:3, col.pop = 4, row.names = 1)
z
print(z, details = TRUE)


###################################################
### code chunk number 8: ReadingFiles.Rnw:212-213
###################################################
getAlleles(z)


###################################################
### code chunk number 9: ReadingFiles.Rnw:218-219
###################################################
attr(z, "locicol")


###################################################
### code chunk number 10: ReadingFiles.Rnw:224-225
###################################################
str(z)


###################################################
### code chunk number 11: ReadingFiles.Rnw:233-234
###################################################
args(read.vcf)


