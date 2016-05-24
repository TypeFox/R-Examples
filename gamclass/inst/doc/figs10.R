## ----setup, cache=FALSE, echo=FALSE-----------------------------------
library(knitr)
options(replace.assign=FALSE,width=72)
opts_chunk$set(fig.path='figs/class-', cache.path='cache/class-',
               fig.align='center', dev='pdf', fig.width=7.5,
               fig.height=6, out.width="0.97\\textwidth",
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

## ----fig10_1, eval=TRUE, echo=TRUE------------------------------------
fig10.1 <- function(){
 ## ---- xWITHerr ----
tau <- (0:5)/2.5; m <- length(tau); n <- 200; SD <- 2
x0 <- rnorm(n, mean=12.5, sd=SD)  # Generate x-values
df <- data.frame(sapply(tau, function(xtau)x0+rnorm(n, sd=SD*xtau)))
  # Columns after the first are x-values with added error
df$y = 15+2.5*x0 + rnorm(n, sd=1.5)
names(df) <- c(paste("X", tau, sep=""), "y")
lab <- c(list("0"),
         lapply(tau[-1], function(x)substitute(A*s[z], list(A=x))))
form <- formula(paste("y ~ ", paste(paste("X", tau, sep=""),
                                  collapse="+")))
library(latticeExtra)
xlabel <- expression(italic(x)*' ('*italic(z)*' with error)')
striplabel <- strip.custom(strip.names=TRUE,
                           var.name="SD(added err)",
                           sep=expression(" = "),
                           factor.levels=as.expression(lab))
gph <- xyplot(form, data=df, outer=TRUE, xlab=xlabel, strip=striplabel,
               type=c("p", "r"))
gph+layer(panel.abline(15, 2.5, lty=2))
}

## ----fig10_2, eval=TRUE, echo=TRUE------------------------------------
fig10.2 <- function(){
set.seed(31)         # Reproduce graph shown
## Use function errorsINx(), from DAAG
errorsINx(gpdiff=4, timesSDx=1.25, SDyerr=2.5, n=80, plotit=FALSE)[["gph"]]
}

## ----pkgs-figs10, eval=TRUE, message=FALSE, warning=FALSE-------------
pkgs <- c("latticeExtra", "DAAG")
z <- sapply(pkgs, require, character.only=TRUE, warn.conflicts=FALSE)
if(any(!z)){
  notAvail <- paste(names(z)[!z], collapse=", ")
  print(paste("The following package requires to be installed:", notAvail))
}

## ----fig10_1x, eval=TRUE, echo=TRUE, fig.width=6.75, fig.height=5-----
if(require('latticeExtra')) fig10.1() else
  print("Package 'latticeExtra' is not available, cannot do graph")

## ----fig10_2x, eval=TRUE, echo=TRUE,  fig.width=6.0, fig.height=3.5----
if(require('latticeExtra')&require('DAAG')) fig10.2() else
  print("Packages 'latticeExtra' &/or 'DAAG' is/are not available, cannot do graph")

