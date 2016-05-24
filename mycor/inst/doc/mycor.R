## ----,echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment=NA)

## ------------------------------------------------------------------------
cor(mtcars[1:5])

## ----,message=FALSE------------------------------------------------------
require(lattice)
require(mycor)
mycor(iris)

## ------------------------------------------------------------------------
out=mycor(iris,alternative="greater", method="kendall",digits=2)
out1=mycor(~mpg+disp+hp+wt,data=mtcars)
summary(out1)

## ----,fig.width=6, fig.height=6------------------------------------------
plot(out)

## ----,fig.width=6, fig.height=6------------------------------------------
plot(out,groups=species,main="Test of mycor::plot")

## ----,fig.width=6, fig.height=6------------------------------------------
plot(out,type=2,groups=spe)

## ----,fig.width=6, fig.height=6------------------------------------------
plot(out1,type=3)

## ----,fig.width=6, fig.height=6------------------------------------------
plot(out,type=4,groups=spe)
plot(out1,type=4,groups=cyl)

