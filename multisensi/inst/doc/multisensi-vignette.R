## ----setup, include=FALSE, cache=FALSE----------------------------
library(knitr)
opts_chunk$set(fig.path='figure/multisensi-vignette-')
options(formatR.arrow=TRUE, width=68)
knit_hooks$set(close.dev = function(before, options, envir) {
    if (!before) graphics.off()})

## ----ME, echo=TRUE, tidy=TRUE-------------------------------------
verhulst <- function(K, Y0, a, t){
    output <- K / (1 + (K/Y0-1)*exp(-a*t))
    return(output) }

## ----MEmulti, echo=TRUE, tidy=TRUE--------------------------------
T <- seq(from=5, to=100, by=5)
verhulst2 <- function(X, t=T){
    out <- matrix(nrow=nrow(X), ncol=length(t), NA)
    for(i in 1:nrow(X)){
        out[i,] <- verhulst(X$K[i], X$Y0[i], X$a[i], t) }
    out <- as.data.frame(out) ; names(out) <- paste("t",t,sep="")
    return(out)	}

## ----MEplot, echo=TRUE, tidy=TRUE, include=TRUE, dev='pdf', close.dev=TRUE, fig.align='center',fig.width=5, fig.height=3.5, out.width='0.9\\linewidth'----
n <- 10 ; set.seed(1234)
X <- data.frame(K=runif(n, min=100,max=1000), Y0=runif(n, min=1, max=40), 
                a=runif(n, min=0.05,max=0.2))
Y <- verhulst2(X)
par(cex.axis=0.7, cex.lab=0.8)
plot(T, Y[1,], type="l", xlab= "Time", ylab="Population size", ylim=c(0,1000))
for(i in 2:n){ lines(T, Y[i,], type="l", col=i) }

## ----MEdynsi1, tidy=FALSE, echo=TRUE------------------------------
library(multisensi)
verhulst.seq <- multisensi(model=verhulst2, reduction=NULL, center=FALSE, 
  design.args = list( K=c(100,400,1000), Y0=c(1,20,40), a=c(0.05,0.1,0.2)))

## ----<MEdynsi2, tidy=TRUE, eval=TRUE------------------------------
print(verhulst.seq, digits=2)

## ----MEdynsi3, tidy=TRUE, echo=-(1),include=TRUE,  dev='pdf',close.dev=TRUE, fig.show='hold', fig.align='center',fig.width=6, fig.height=6, out.width='0.45\\linewidth'----
par(cex.axis=1, cex.lab=0.9)
# color palettes: rainbow, heat.colors, terrain.colors, topo.colors, cm.colors
plot(verhulst.seq, normalized=TRUE, color=terrain.colors)
title(xlab="Time in half-decades.")
plot(verhulst.seq, normalized=FALSE, color=terrain.colors)
title(xlab="Time in half-decades.")

## ----MEalteruni, echo=TRUE, tidy=FALSE----------------------------
X <- expand.grid(K=c(100,400,1000), Y0=c(1,20,40), a=c(0.05,0.1,0.2))
Y <- verhulst2(X) ## this part can be performed outside R if necessary
verhulst.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE)

## ----MEgsi, tidy=FALSE--------------------------------------------
## Note that this old syntax of multisensi still works:
## verhulst.gsi <- gsi(formula=2, Y, X) 
verhulst.pca <- multisensi(design=X, model=Y, reduction=basis.ACP, scale=FALSE)

## ----MEgsisummary, tidy=TRUE--------------------------------------
summary(verhulst.pca, digits=2)

## ----MEgsiplot1, tidy=TRUE, echo=-(1), include=TRUE, dev='pdf', close.dev=TRUE,  fig.width=6, fig.height=6, out.width='0.9\\linewidth', fig.align='center'----
par(cex.axis=0.8, cex.lab=0.9)
plot(verhulst.pca, graph=1)

## ----MEgsiplot2, tidy=TRUE, include=TRUE, dev='pdf', close.dev=TRUE,  fig.width=6, fig.height=6, out.width='0.6\\linewidth', fig.align='center'----
plot(verhulst.pca, graph=2)

## ----MEgsiplot3, tidy=TRUE, include=TRUE, dev='pdf', close.dev=TRUE,  fig.width=6, fig.height=6, out.width='0.6\\linewidth', fig.align='center'----
plot(verhulst.pca, graph=3)

## ----MEgsipoly2, echo=TRUE, tidy=FALSE----------------------------
verhulst.poly <- multisensi(design = X, model = Y, reduction = basis.poly,
      dimension = 0.99, center = FALSE, scale = FALSE, cumul = FALSE,
      basis.args = list(degree=6, x.coord=T), analysis = analysis.anoasg, 
      analysis.args = list(formula=2, keep.outputs=FALSE))
summary(verhulst.poly, digits=2)

## ----MEgsipolyplot, tidy=TRUE, echo=-(1), include=TRUE, dev='pdf', close.dev=TRUE,  fig.align='center',fig.width=9, fig.height=6, out.width='0.9\\linewidth'----
par(cex.axis=0.9, cex.lab=1.0)
plot(verhulst.poly, nb.comp=3,graph=1)

## ----MEgsi2.bsplines, tidy=FALSE----------------------------------
## bsplines
verhulst.bspl <- multisensi(design=X, model=Y, reduction=basis.bsplines, 
                        dimension=NULL, center=FALSE, scale=FALSE, 
                        basis.args=list(knots=10, mdegree=3), cumul=FALSE, 
                        analysis=analysis.anoasg, 
                        analysis.args=list(formula=2, keep.outputs=FALSE))

## ----MEgsibsplinesplot, tidy=TRUE, echo=-(1), include=TRUE, dev='pdf', close.dev=TRUE, fig.align='center',fig.width=15, fig.height=6, out.width='0.9\\linewidth'----
#par(cex.axis=0.9, cex.lab=1.0)
plot(verhulst.bspl, nb.comp=5,graph=1)

## ----MEsensitivity, tidy=TRUE, echo=TRUE, include=TRUE------------
library(sensitivity)

## ----MEdynsiSobol, tidy=TRUE, echo=TRUE, include=TRUE-------------
m <- 10000
Xb <- data.frame(K=runif(m, min=100,max=1000), Y0=runif(m, min=1, max=40), 
                   a=runif(m, min=0.05,max=0.2))
verhulst.seq.sobol<-
    multisensi(design=sobol2007, model=verhulst2,
               reduction=NULL, analysis=analysis.sensitivity, center=TRUE, 
               design.args=list(X1=Xb[1:(m/2),], X2=Xb[(1+m/2):m,], nboot=100), 
               analysis.args=list(keep.outputs=FALSE))
print(verhulst.seq.sobol,digits=2)

## ----MEdynsiSobolplot, tidy=TRUE, echo=-(1), include=TRUE, dev='pdf', close.dev=TRUE, fig.align='center',fig.width=6, fig.height=6, out.width='0.7\\linewidth'----
  par(cex.axis=0.8, cex.lab=0.9)
  plot(verhulst.seq.sobol, normalized=TRUE, color=terrain.colors)
  title(xlab="Time in half-decades")

## ----MEdynsiFast, tidy=FALSE, echo=TRUE, include=TRUE-------------
verhulst.seq.fast <- multisensi(design = fast99, model = verhulst2,
      center = FALSE, reduction = NULL, analysis = analysis.sensitivity, 
      design.args=list( factors=c("K","Y0","a"), n=1000, q = "qunif", 
        q.arg = list(list(min=100, max=1000), list(min=1, max=40), 
          list(min = 0.05, max = 0.2))),
      analysis.args=list(keep.outputs=FALSE))
  print(verhulst.seq.fast,digits=2)

## ----MEdynsiFastplot, tidy=TRUE, echo=-(1), include=TRUE, dev='pdf', close.dev=TRUE, fig.align='center',fig.width=6, fig.height=6, out.width='0.7\\linewidth'----
  par(cex.axis=0.8, cex.lab=0.9)
  plot(verhulst.seq.fast, normalized=TRUE, color=terrain.colors)
  title(xlab="Time in half-decades")

