### R code from vignette source 'RGIFT.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: RGIFT.Rnw:118-120
###################################################
library(RGIFT)
source("MCQExample.R")


###################################################
### code chunk number 2: RGIFT.Rnw:149-150
###################################################
source("TFQExample.R")


###################################################
### code chunk number 3: RGIFT.Rnw:160-161
###################################################
source("SAQExample.R")


###################################################
### code chunk number 4: RGIFT.Rnw:173-174
###################################################
source("MQExample.R")


###################################################
### code chunk number 5: RGIFT.Rnw:185-186
###################################################
source("MWQExample.R")


###################################################
### code chunk number 6: RGIFT.Rnw:201-202
###################################################
source("NQExample.R")


###################################################
### code chunk number 7: RGIFT.Rnw:219-220
###################################################
source("EExample.R")


###################################################
### code chunk number 8: RGIFT.Rnw:232-233
###################################################
source("DExample.R")


###################################################
### code chunk number 9: RGIFT.Rnw:274-275
###################################################
GIFTparse("Take care with $, {, } and ~.")


###################################################
### code chunk number 10: RGIFT.Rnw:287-288
###################################################
source(file("UTF8.R", encoding="ISO-8859-1"))


###################################################
### code chunk number 11: RGIFT.Rnw:311-337
###################################################
#Image repository
repos<-"http://www.uclm.es/profesorado/vgomez/imgtests/"

#Create image
imgfile<-"q1.png"
png(imgfile, width=480, height=240)

par(mfrow=c(1,4))
#Sample several random variables and plot histogram
hist(rbinom(100, 10, .25), main="(a)", xlab="x")
hist(rnorm(100, 10*.25, sqrt(10*.25*.75)), main="(b)", xlab="x")
hist(rexp(100, 1/(10*.25)), main="(c)", xlab="x")
hist(runif(100, 8, 12), main="(d)", xlab="x")

dev.off()

#Create question in HTML format
qtxt<-paste("[html]What\'s the histogram for a Binomial ",
   "variable with n=10 and prob=0.25?\n<p>\n<img src=\'", 
   repos, imgfile, "\'>\n<p>", sep="")

#Print question
sink(file="imgqtion.txt", type="output")
GIFTMC(qtxt, c("(a)", "(b)", "(c)", "(d)"), rightans=1, wwrong="-33.333")
sink()



###################################################
### code chunk number 12: RGIFT.Rnw:356-361
###################################################
sink(file="eqqtion.txt")
qtxt<-paste("[html]The expectation of a random variable which is \n",
   "Binomial with $$n=10$$ and $$\\pi=0.7$$ is ", sep="")
GIFTMC(qtxt, c("0.7", "1", "7", "0.07"), rightans=3)
sink()


