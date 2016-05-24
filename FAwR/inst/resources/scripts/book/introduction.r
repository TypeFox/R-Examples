### R code from vignette source 'introduction.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")
if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(xtable, quietly=TRUE)) install.packages("xtable")
#rm( list=ls() )


###################################################
### code chunk number 2: introduction.rnw:208-211
###################################################
a <- rnorm( 1000 )
b <- mean( a )
b


###################################################
### code chunk number 3: introduction.rnw:271-272
###################################################
1 + 2


###################################################
### code chunk number 4: introduction.rnw:279-280
###################################################
a.out <- 1 + 2


###################################################
### code chunk number 5: introduction.rnw:291-296
###################################################
hi.there <- function() { 
  a.out <- 1 + 2
  cat("Hello World!\n")
  return(a.out)
}


###################################################
### code chunk number 6: introduction.rnw:303-304
###################################################
hi.there()


###################################################
### code chunk number 7: introduction.rnw:311-317
###################################################
hi.there <- function(arbitrary.number) { 
  a.out <- 1 + arbitrary.number
  cat("Hello World!\n")
  return(a.out)
}
hi.there(pi)


###################################################
### code chunk number 8: introduction.rnw:331-332
###################################################
ls()


###################################################
### code chunk number 9: introduction.rnw:337-339
###################################################
rm(a.out)
ls()


###################################################
### code chunk number 10: introduction.rnw:347-348
###################################################
rm(list=ls())


###################################################
### code chunk number 11: introduction.rnw:541-542
###################################################
mean(c(1,2,3))


###################################################
### code chunk number 12: introduction.rnw:551-552
###################################################
library(boot)


###################################################
### code chunk number 13: Setup
###################################################
rm(list=ls())                


###################################################
### code chunk number 14: Setup
###################################################
sprintf("%.1f", 3.01)


###################################################
### code chunk number 15: introduction.rnw:643-646
###################################################
class(mean)
class(c(1,2,3))
class(mean(c(1,2,3)))


###################################################
### code chunk number 16: introduction.rnw:651-654
###################################################
not.really.numeric <- c("1","2","3")
class(not.really.numeric)
mean(not.really.numeric)


###################################################
### code chunk number 17: introduction.rnw:659-662
###################################################
not.really.numeric <- as.numeric(not.really.numeric)
class(not.really.numeric)
mean(not.really.numeric)


###################################################
### code chunk number 18: introduction.rnw:673-674
###################################################
str(not.really.numeric)


