## ----echo=FALSE----------------------------------------------------------
library(highlight)
library(knitr)
##homedir <- '/home/rmcd/tex/d67/Rtutorial/'
options(digits=4)
figsize <- 4.5
opts_chunk$set(size='footnotesize', prompt=FALSE, comment=NA, 
               fig.align='center', fig.width = figsize, 
               fig.height=figsize, out.width='3.75in')

#              , fig.width=4.5*3.75/3.25, fig.height=4.5, 
#              , out.width='3.75in', out.height='3.25in'
#               )
opts_knit$set(highlight = TRUE, 
              eval.after='fig.cap',
              prompt=TRUE,
              renderer=renderer_latex(document=FALSE),
              size='footnotesize') 

## ----echo=FALSE----------------------------------------------------------
library(knitr)
library(derivmkts)
library(markdown)

opts_chunk$set(collapse=TRUE)

## ------------------------------------------------------------------------
s <- 100; k <- 100; r <- 0.08; v <- 0.30; tt <- 2; d <- 0
bscall(s, k, v, r, tt, d)
bsput(s, k, v, r, tt, d)
bsput(s, c(95, 100, 105), v, r, tt, d)
bsopt(s, c(95, 100, 105), v, r, tt, d)$Call

## ------------------------------------------------------------------------
s <- 100; k <- 100; r <- 0.08; v <- 0.30; tt <- 2; d <- 0.03
binomopt(s, k, v, r, tt, d, nstep=4)
binomopt(s, k, v, r, tt, d, nstep=4, returnparams=TRUE)
binomopt(s, k, v, r, tt, d, nstep=4, putopt=TRUE)
binomopt(s, k, v, r, tt, d, nstep=4, returntrees=TRUE, putopt=TRUE)

## ------------------------------------------------------------------------
H <- 105
uicall(c(95, 100, 105), k, v, r, tt, d, H)
bscall(c(95, 100, 105), k, v, r, tt, d)

## ------------------------------------------------------------------------
H <- 105
greeks(uicall(s, k, v, r, tt, d, H))
greeks2(uicall, s=s, k=k, v=v, r=r, tt=tt, d=d, H=H)


## ----binomplot1, fig.cap='Basic option plot showing stock prices and nodes at which the option is exercised.\\label{fig:binomplot1}'----
binomplot(s, k, v, r, tt, d, nstep=6, american=TRUE, putopt=TRUE)


## ----binomplot2, fig.cap='Same plot as Figure \\ref{fig:binomplot1} except that values and arrows are added to the plot.\\label{fig:binomplot2}'----
binomplot(s, k, v, r, tt, d, nstep=6, american=TRUE, putopt=TRUE,
    plotvalues=TRUE, plotarrows=TRUE)

## ----binomplot3, fig.cap="Binomial plot when nstep is 40.\\label{fig:binomplot3}"----
d <- 0.06
binomplot(s, k, v, r, tt, d, nstep=40, american=TRUE)

## ----binomplot4, fig.cap="Binomial plot when nstep is 40 using the argument ylimval to focus on a subset.\\label{fig:binomplot4}"----
d <- 0.06
binomplot(s, k, v, r, tt, d, nstep=40, american=TRUE, ylimval=c(75, 225))

