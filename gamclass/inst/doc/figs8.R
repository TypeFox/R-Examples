## ----setup, cache=FALSE, echo=FALSE, include=FALSE--------------------
library(knitr)
options(replace.assign=FALSE,width=72)
opts_chunk$set(fig.path='figs/tree-', cache.path='cache/tree-',
               fig.align='center', dev='pdf', fig.width=3.5,
               fig.height=3.5, out.width="0.8\\textwidth",
               fig.show='hold', par=TRUE,
               tidy=FALSE,  comment=NA)
knit_hooks$set(pars=function(before, options, envir){
if (before && options$fig.show!='none') par(mar=c(4,4,1.6,.1),
              font.main=1, cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),
              tcl=-.3)
              par(options$pars)
}, crop=hook_pdfcrop)
pdf.options(pointsize=12)
oldopt <- options(digits=4)

## ----fig8_1A, eval=TRUE, echo=TRUE------------------------------------
fig8.1A <- function(){
  if(!exists('car90.rpart'))
  car90.rpart <- rpart(Mileage ~ tonsWt, data=Car90)
  plot(car90.rpart)
text(car90.rpart, xpd=TRUE, digits=3)
mtext(side=3, line=1.25, "A: Regression tree", adj=0)
}

## ----fig8_1B, eval=TRUE, echo=TRUE------------------------------------
fig8.1B <- function(){
if(!exists('car90.rpart'))
  car90.rpart <- rpart(Mileage ~ tonsWt, data=Car90)
plot(Mileage ~ tonsWt, data=Car90)
wt <- with(Car90, tonsWt)
hat <- predict(car90.rpart)
addhlines(wt, hat, lwd=2, col="gray")
mtext(side=3, line=1.25, "B: Predicted values from tree", adj=0)
}

## ----fig8_2, eval=TRUE, echo=TRUE-------------------------------------
fig8.2 <- function(){
BSS <- bssBYcut(tonsWt, Mileage, Car90)
with(BSS, plot(xOrd, bss, xlab="Cutpoint",
               ylab="Between groups sum of squares"))
abline(v=1.218, lty=2)
}

## ----fig8_3A, eval=TRUE, echo=TRUE------------------------------------
fig8.3A <- function(){
opar <- par(mar=c(4,4,2.6,1.6))
if(!exists('car90x.rpart'))
  car90x.rpart <- rpart(Mileage ~ tonsWt, data=Car90,
                        minbucket=5, minsplit=10,
                        cp=0.001)
plot(car90x.rpart, uniform=TRUE)
text(car90x.rpart, digits=3, xpd=TRUE)
mtext(side=3, line=0.75, "A: Decision tree", adj=0)
par(opar)
}

## ----fig8_3B, eval=TRUE, echo=TRUE------------------------------------
opar <- par(mar=c(4,4,2.6,1.6))
fig8.3B <- function(){
if(!exists('car90x.rpart'))
  car90x.rpart <- rpart(Mileage ~ tonsWt, data=Car90,
                        minbucket=5, minsplit=10,
                        cp=0.001)
plot(Mileage ~ tonsWt, data=Car90)
hat <- predict(car90x.rpart)
wt <- with(Car90, tonsWt)
addhlines(wt, hat, lwd=2, col="gray")
mtext(side=3, line=0.75, "B: Mileage vs tonsWt", adj=0)
par(opar)
}

## ----fig8_4, eval=TRUE, echo=TRUE-------------------------------------
fig8.4 <- function(){
if(!exists('car90x.rpart'))
  car90x.rpart <- rpart(Mileage ~ tonsWt, data=Car90,
                        minbucket=5, minsplit=10,
                        cp=0.001)
plotcp(car90x.rpart)
}

## ----fig8_5, eval=TRUE, echo=TRUE-------------------------------------
fig8.5 <- function(){
if(!exists('car90.rf'))
  car90.rf <- randomForest(Mileage ~ tonsWt,
                          data=Car90)
plot(Mileage ~ tonsWt, data=Car90, type="n")
with(Car90, points(Mileage ~ tonsWt, cex=0.8))
hat <- predict(car90.rf)
with(Car90, points(hat ~ tonsWt, pch="-"))
}

## ----fig8_6, eval=TRUE, echo=TRUE-------------------------------------
fig8.6 <- function(){
ran <- range(errsmat)
at <- round(ran+c(0.02,-0.02)*diff(ran),2)
lis <- list(limits=ran, at=at, labels=format(at, digits=2))
lims=list(lis,lis,lis,lis,lis,lis)
library(lattice)
splom(errsmat,
      pscales=lims,
      par.settings=simpleTheme(cex=0.75),
      col=adjustcolor("black", alpha=0.5),
      panel=function(x,y,...){lpoints(x,y,...)
      panel.abline(0,1,col="gray")}
)
}

## ----pkgs-figs8, eval=TRUE, message=FALSE, warning=FALSE--------------
pkgs <- c("rpart","mgcv","randomForest","gamclass")
z <- sapply(pkgs, require, character.only=TRUE, warn.conflicts=FALSE)
if(any(!z)){
  notAvail <- paste(names(z)[!z], collapse=", ")
  print(paste("The following packages should be installed:", notAvail))
}

## ----Car90------------------------------------------------------------
if(!exists('Car90'))
Car90 <- na.omit(car90[, c("Mileage","Weight")])
## Express weight in metric tonnes
Car90 <- within(Car90, tonsWt <- Weight/2240)

## ----meuse------------------------------------------------------------
getmeuse <- function(){
if(require('sp', quietly=TRUE)){
data("meuse", package="sp", envir = environment())
meuse <- within(meuse, {levels(soil) <- c("1","2","2")
                        ffreq <- as.numeric(ffreq)
                        loglead <- log(lead)
                        })
invisible(meuse)
} else if(!exists("meuse"))
  print("Dataset 'meuse' was not found, get from package 'sp'")
}

## ----meuse-rf---------------------------------------------------------
cfRF <- function(nrep=50){
form1 <- ~ dist + elev + soil + ffreq
form3 <- ~ s(dist, k=3) + s(elev,k=3) + soil +ffreq
form3x <- ~ s(dist, k=3) + s(elev,k=3) + s(x, k=3) + soil+ffreq
form8x <- ~ s(dist, k=8) + s(elev,k=8) + s(x, k=8) + soil+ffreq
formlist <- list("Hybrid1"=form1, "Hybrid3"=form3,
                 "Hybrid3x"=form3x, "Hybrid8x"=form8x)
## ----rfgam-setup----
rfVars <- c("dist", "elev", "soil", "ffreq", "x", "y")
errsmat <- matrix(0, nrep, length(formlist)+2)
dimnames(errsmat)[[2]] <- c(names(formlist), "rfTest", "rfOOB")
n <- 95
for(i in 1:nrep){
sub <- sample(1:nrow(meuse), n)
meuseOut <- meuse[-sub,]
meuseIn <- meuse[sub,]
errsmat[i, ] <- gamRF(formlist=formlist, yvar="loglead",
                      rfVars=rfVars,
                      data=meuseIn, newdata=meuseOut,
                      printit=FALSE)
}
invisible(errsmat)
}

## ----car90-plot12, fig.width=3.35, fig.height=3, out.width="0.47\\textwidth", echo=FALSE----
fig8.1A()
fig8.1B()

## ----bss-plot, echo=FALSE, fig.width=3.35, fig.height=3, out.width="0.8\\textwidth"----
fig8.2()

## ----Car90-loosen, echo=FALSE, fig.width=3.5, fig.height=2.75, out.width="0.47\\textwidth"----
fig8.3A()
fig8.3B()

## ----car90x-plotcp, fig.width=4.0, fig.height=3.25, out.width="0.92\\textwidth"----
fig8.4()

## ----Car90-rf, eval=TRUE, echo=FALSE, fig.width=4.0, fig.height=3.5----
fig8.5()

## ----cf-models-rfgam, echo=FALSE, fig.width=7, fig.height=7, out.width="0.8\\textwidth"----
nrep <- NA
meuse <- getmeuse()
if(!exists('errsmat'))errsmat <- cfRF(nrep=25)
fig8.6()
if(exists('errsmat'))nrep <- nrow(errsmat)

