#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
## Time-stamp: <[createData.R] by DSB Die 10/03/2009 11:06 (GMT) on daniel@puc-home>
##
## Description:
## Create the data for the example BayesX run.
##
## History:
## 10/03/2009   file creation
#####################################################################################

### setup
nObs <- 300
set.seed(991)


### create covariates and their effects

## smooth functions in x1 and x2
x1 <- round(runif(n=nObs,
                  min=-pi, max=pi),
            2)
x1Effect <- sin(x1) * cos(x1)
plot(x1Effect[order(x1)] ~ x1[order(x1)],
     type="l")

x2 <- round(runif(n=nObs),
            2)
x2Effect <- (x2 - 0.5)^2
plot(x2Effect[order(x2)] ~ x2[order(x2)],
     type="l")

## linear functions in x3 and x4
x3 <- rnorm(n=nObs)
x3Effect <- 9.2 * x3

x4 <- rexp(n=nObs)
x4Effect <- 5.1 * x4

## spatial effect from the district in Tanzania
library(BayesX)
tanzania <- read.bnd(file="tanzania.bnd")

tanzaniaEffects <- rnorm(n=length(tanzania))
names(tanzaniaEffects) <- names(tanzania)

drawmap(map=tanzania,
        data=
        data.frame(x=names(tanzaniaEffects),
                   y=tanzaniaEffects),
        regionvar="x",
        plotvar="y")

district <- sample(x=names(tanzania),
                   size=nObs,
                   replace=TRUE)
districtEffect <- tanzaniaEffects[district]


### now generate the response

linearPredictor <- x1Effect + x2Effect + x3Effect + x4Effect + districtEffect
y <- linearPredictor + rnorm(n=nObs)


### write data into text file

data <- data.frame(x1=x1,
                   x2=x2,
                   x3=x3,
                   x4=x4,
                   district=district,
                   y=y)
write.table(x=data, file="data.txt",
            quote=FALSE, col.names=TRUE, row.names=FALSE)


