# ---------------------------------------------------------------- #

plotpos <- function(x, a=0, orient="xF", ...) {

  # INPUT
  # x = colonna
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta cartesiana della serie

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=plotpos
  if (orient=="xF") plot(x,F,...)
  else if (orient=="Fx") plot(F,x,...)
  else stop("plotpos(x, a, orient): orient unknown") 
  grid()
}


# ---------------------------------------------------------------- #

pointspos <- function(x, a=0, orient="xF", ...) {

  # INPUT
  # x = colonna
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta cartesiana della serie

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)

  x=ordinato
  F=plotpos
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("pointspos(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

loglogplot <- function(x, a=0, orient="xF", ...) {

  # INPUT
  # x = colonna
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione log-log

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)

  x=ordinato
  F=plotpos
  if (orient=="xF") plot(log(x),log(1-F),...)
  else if (orient=="Fx") plot(log(1-F),log(x),...)
  else stop("loglogplot(x, a, orient): orient unknown") 
  grid()
}


# ---------------------------------------------------------------- #

loglogpoints <- function(x, a=0, orient="xF", ...) {

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)

  x=ordinato
  F=plotpos
  if (orient=="xF") points(log(x),log(1-F),...)
  else if (orient=="Fx") points(log(1-F),log(x),...)
  else stop("loglogpoints(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

normplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica normale della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  s <- sd(ordinato)
  u <- (ordinato-m)/s
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qnorm(plotpos)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   grid(nx=NULL, ny=NA)
   abline(h=qq, col="gray", lty="dotted")
   if(line==TRUE) lines(ordinato,u,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   grid(ny=NULL, nx=NA)
   abline(v=qq, col="gray", lty="dotted")
   if(line==TRUE) lines(u,ordinato,lty=2)
  }
  else stop("normplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

normpoints <- function(x, a=0, orient="xF", ...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica normale della serie

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qnorm(plotpos)

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("normpoints(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

lognormplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica lognormale della serie

  ordinato <- sort(x)
  logordinato <- log(ordinato)
  m <- mean(logordinato)
  s <- sd(logordinato)
  u <- (logordinato-m)/s
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qnorm(plotpos)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,log="x",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
   if (line==TRUE) lines(ordinato,u,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
   if (line==TRUE) lines(u, ordinato, lty=2)
  }
  else stop("lognormplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

studentplot <- function(x, df, a=0, orient="xF", line=FALSE,...) {

  # INPUT
  # df = degrees of freedom

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qt(plotpos, df)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qt(pp, df)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   grid(nx=NULL, ny=NA)
   abline(h=qq, col="gray", lty="dotted")
   #if(line==TRUE) lines(ordinato,u,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   grid(ny=NULL, nx=NA)
   abline(v=qq, col="gray", lty="dotted")
   #if(line==TRUE) lines(u,ordinato,lty=2)
  }
  else stop("studentplot(x, df, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

studentpoints <- function(x, df, a=0, orient="xF", ...) {

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qt(plotpos, df)

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("studentpoints(x, df, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

logisplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  ordinato <- sort(x)
  mu <- mean(ordinato)
  s <- sqrt(3)*sd(ordinato)/pi
  adim <- (ordinato - mu)/s

  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qlogis(plotpos)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qlogis(pp)

  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   grid(nx=NULL, ny=NA)
   abline(h=qq, col="gray", lty="dotted")
   if(line==TRUE) lines(ordinato,adim,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   grid(ny=NULL, nx=NA)
   abline(v=qq, col="gray", lty="dotted")
   if(line==TRUE) lines(adim,ordinato,lty=2)
  }
  else stop("logisplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

logispoints <- function(x, a=0, orient="xF", ...) {

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qlogis(plotpos)

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("logispoints(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

gammaplot <- function(x, shape, a=0, orient="xF", line=FALSE,...) {

  # INPUT
  # shape of the gamma distribution

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qgamma(plotpos, shape)
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qgamma(pp, shape)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   grid(nx=NULL, ny=NA)
   abline(h=qq, col="gray", lty="dotted")
   #if(line==TRUE) lines(ordinato,u,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   grid(ny=NULL, nx=NA)
   abline(v=qq, col="gray", lty="dotted")
   #if(line==TRUE) lines(u,ordinato,lty=2)
  }
  else stop("gammaplot(x, shape, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

gammapoints <- function(x, shape, a=0, orient="xF", ...) {

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- qgamma(plotpos, shape)

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("gammapoints(x, shape, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

gumbelplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica di Gumbel della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  s <- sd(ordinato)
  beta <- (s*sqrt(6))/pi
  mu <- m - 0.5772*beta
  adim <- (ordinato-mu)/beta
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- -log(-log(plotpos))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- -log(-log(pp))

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
   if (line==TRUE) lines(ordinato,adim,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
   if (line==TRUE) lines(adim, ordinato, lty=2)
  }
  else stop("gumbelplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

gumbelpoints <- function(x, a=0, orient="xF", ...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica di Gumbel della serie

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- -log(-log(plotpos))

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("gumbelpoints(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

frechetplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  ordinato <- sort(x)
  logordinato <- log(ordinato)
  m <- mean(logordinato)
  s <- sd(logordinato)
  beta <- (s*sqrt(6))/pi
  mu <- m - 0.5772*beta
  adim <- (logordinato-mu)/beta
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- -log(-log(plotpos))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .975, .99, .995, .999)
  qq <- -log(-log(pp))

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,log="x",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
   if (line==TRUE) lines(ordinato,adim,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
   if (line==TRUE) lines(adim, ordinato, lty=2)
  }
  else stop("frechetplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

weibullplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- log(-log(1 - plotpos))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .975, .99, .995, .999)
  qq <- log(-log(1 - pp))

  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,log="x",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
   #if (line==TRUE) lines(ordinato,adim,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
   #if (line==TRUE) lines(adim, ordinato, lty=2)
  }
  else stop("weibullplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

weibullpoints <- function(x, a=0, orient="xF", ...) {

  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- log(-log(1 - plotpos))

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("weibullpoints(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

unifplot <- function (x, a=0, orient="xF", line=FALSE, ...) {
    ordinato <- sort(x)
    u = (ordinato - min(ordinato))/(max(ordinato) - min(ordinato))
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- qunif(plotpos)
    pp <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
    qq <- qunif(pp)
    x = ordinato
    F = uteorica
    if (orient=="xF") {
     plot(x, F, axes = FALSE, ...)
     axis(1)
     axis(2, at = qq, labels = pp)
     abline(h = qq, col = "gray", lty = "dotted")
     grid(nx = NULL, ny = NA)
     if (line == TRUE) lines(ordinato, u, lty = 2)
    }
    else if (orient=="Fx") {
     plot(F, x, axes = FALSE, ...)
     axis(2)
     axis(1, at = qq, labels = pp)
     abline(v = qq, col = "gray", lty = "dotted")
     grid(ny = NULL, nx = NA)
     if (line == TRUE) lines(u, ordinato, lty = 2)
    }
    else stop("unifplot(x, a, orient, line): orient unknown")
    box()
}


# ---------------------------------------------------------------- #

unifpoints <- function (x, a=0, orient="xF", ...) {
 pointspos(x, a=0, orient="xF", ...)
}


# ---------------------------------------------------------------- #

expplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  # INPUT
  # x = colonna
  # line = linea rappresentante la distribuzione
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta probabilistica esponenziale della serie

  ordinato <- sort(x)
  m <- mean(ordinato)
  lambda <- 1/m
  adim <- lambda*ordinato
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- -log(1 - plotpos)
  pp <- c(.1, .25, .50, .75, .9, .95, .975, .99, .995, .999)
  qq <- -log(1 - pp)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
   if (line==TRUE) lines(ordinato,adim,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
   if (line==TRUE) lines(adim, ordinato, lty=2)
  }
  else stop("expplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

exppoints <- function(x, a=0, orient="xF", ...) {
  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- -log(1 - plotpos)

  x=ordinato
  F=uteorica
  if (orient=="xF") points(x,F,...)
  else if (orient=="Fx") points(F,x,...)
  else stop("gumbelpoints(x, a, orient): orient unknown")
}


# ---------------------------------------------------------------- #

paretoplot <- function(x, a=0, orient="xF", line=FALSE,...) {

  ordinato <- sort(x)
  logordinato <- log(ordinato)
  m <- mean(logordinato)
  lambda <- 1/m
  adim <- lambda*logordinato
  n <- length(logordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)
  uteorica <- -log(1 - plotpos)
  pp <- c(.1, .25, .50, .75, .9, .95, .975, .99, .995, .999)
  qq <- -log(1 - pp)

  #x11(width = 7, height = 7, pointsize = 12)
  x=ordinato
  F=uteorica
  if (orient=="xF") {
   plot(x,F,log="x",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
   if (line==TRUE) lines(ordinato,adim,lty=2)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
   if (line==TRUE) lines(adim, ordinato, lty=2)
  }
  else stop("expplot(x, a, orient, line): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

regionalplotpos <- function(x, cod, a=0, orient="xF", ...) {

  # INPUT
  # x = colonna
  # cod = .
  # ... = graphical parameters as xlab, ylab, main...
  # OUTPUT
  # rappresentazione in carta cartesiana della serie

  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)

  x=c(min(x),max(x))
  F=c(0,1)
  if (orient=="xF") plot(x,F,type="n",...)
  else if (orient=="Fx") plot(F,x,type="n",...)
  else stop("regionalplotpos(x, cod, a, orient): orient unknown")
  for(j in 1:k) {
   ordinato <- sort(X[cod==liv[j]])
   n <- length(ordinato)
   i <- 1:n
   plotpos <- (i - a)/(n + 1 - 2*a)
   if (orient=="xF") points(ordinato, plotpos, pch=j, col=j)
   else if (orient=="Fx") points(plotpos, ordinato, pch=j, col=j)
  }
  grid()
}


# --------------------------------------------------------------------------- #

regionalnormplot <- function(x, cod, a=0, orient="xF", ...) {
  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)
  x=c(min(X),max(X))
  F=c(qnorm((1-a)/(max(ni)+1-2*a)), qnorm((max(ni)-a)/(max(ni)+1-2*a)))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)
  if (orient=="xF") {
   plot(x,F,type="n",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- qnorm(plotpos)
    points(ordinato,uteorica,pch=j,col=j)
   }
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
  }
  else if (orient=="Fx") {
   plot(F,x,type="n",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- qnorm(plotpos)
    points(uteorica,ordinato,pch=j,col=j)
   }
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
  }
  else stop("regionalnormplot(x, cod, a, orient): orient unknown")
  box()
}


# --------------------------------------------------------------------------- #

regionallognormplot <- function(x, cod, a=0, orient="xF", ...) {
  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(qnorm((1-a)/(max(ni)+1-2*a)), qnorm((max(ni)-a)/(max(ni)+1-2*a)))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- qnorm(pp)
  if (orient=="xF") {
   plot(x,F,log="x",type="n",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- qnorm(plotpos)
    points(ordinato,uteorica,pch=j,col=j)
   }
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",type="n",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- qnorm(plotpos)
    points(uteorica,ordinato, pch=j,col=j)
   }
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
  }
  else stop("regionallognormplot(x, cod, a, orient): orient unknown")
  box()
}


# ---------------------------------------------------------------------- #

regionalgumbelplot <- function(x, cod, a=0, orient="xF", ...) {
  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(-log(-log((1-a)/(max(ni)+1-2*a))), -log(-log((max(ni)-a)/(max(ni)+1-2*a))))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- -log(-log(pp))
  if (orient=="xF") {
   plot(x,F,type="n",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(-log(plotpos))
    points(ordinato,uteorica,pch=j,col=j)
   }
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
  }
  else if (orient=="Fx") {
   plot(F,x,type="n",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(-log(plotpos))
    points(uteorica,ordinato,pch=j,col=j)
   }
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
  }
  else stop("regionalgumbelplot(x, cod, a, orient): orient unknown")
  box()
}


# ---------------------------------------------------------------------- #

regionalfrechetplot <- function(x, cod, a=0, orient="xF", ...) {
  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(-log(-log((1-a)/(max(ni)+1-2*a))), -log(-log((max(ni)-a)/(max(ni)+1-2*a))))
  pp <- c(.001, .01, .05, .1, .25, .50, .75, .9, .95, .99, .999)
  qq <- -log(-log(pp))
  if (orient=="xF") {
   plot(x,F,log="x",type="n",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(-log(plotpos))
    points(ordinato,uteorica,pch=j,col=j)
   }
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",type="n",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(-log(plotpos))
    points(uteorica,ordinato,pch=j,col=j)
   }
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
  }
  else stop("regionalfrechetplot(x, cod, a, orient): orient unknown")
  box()
}


# ---------------------------------------------------------------------- #

regionalexpplot <- function(x, cod, a=0, orient="xF", ...) {
  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(-log(1 - ((1-a)/(max(ni)+1-2*a))), -log(1 - ((max(ni)-a)/(max(ni)+1-2*a))))
  pp <- c(.1, .25, .50, .75, .9, .95, .975, .99, .995, .999)
  qq <- -log(1 - pp)
  if (orient=="xF") {
   plot(x,F,type="n",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(1 - plotpos)
    points(ordinato,uteorica,pch=j,col=j)
   }
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
  }
  else if (orient=="Fx") {
   plot(F,x,type="n",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(1 - plotpos)
    points(uteorica,ordinato,pch=j,col=j)
   }
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
  }
  else stop("regionalexpplot(x, cod, a, orient): orient unknown")
  box()
}


# ---------------------------------------------------------------------- #

regionalparetoplot <- function(x, cod, a=0, orient="xF", ...) {
  X=x
  cod <- factor(cod)
  k <- nlevels(cod)
  liv <- levels(cod)
  ni <- tapply(X,cod,length)

  x=c(min(X),max(X))
  F=c(-log(1 - ((1-a)/(max(ni)+1-2*a))), -log(1 - ((max(ni)-a)/(max(ni)+1-2*a))))
  pp <- c(.1, .25, .50, .75, .9, .95, .975, .99, .995, .999)
  qq <- -log(1 - pp)
  if (orient=="xF") {
   plot(x,F,log="x",type="n",axes=FALSE,...)
   axis(1)
   axis(2,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(1 - plotpos)
    points(ordinato,uteorica,pch=j,col=j)
   }
   abline(h=qq, col="gray", lty="dotted")
   grid(nx=NULL, ny=NA)
  }
  else if (orient=="Fx") {
   plot(F,x,log="y",type="n",axes=FALSE,...)
   axis(2)
   axis(1,at=qq,labels=pp)
   for(j in 1:k) {
    ordinato <- sort(X[cod==liv[j]])
    n <- length(ordinato)
    i <- 1:n
    plotpos <- (i - a)/(n + 1 - 2*a)
    uteorica <- -log(1 - plotpos)
    points(uteorica,ordinato,pch=j,col=j)
   }
   abline(v=qq, col="gray", lty="dotted")
   grid(ny=NULL, nx=NA)
  }
  else stop("regionalparetoplot(x, cod, a, orient): orient unknown")
  box()
}


# ---------------------------------------------------------------- #

plotposRP <- function(x, a=0, orient="xF", ...) {
  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)

  x=ordinato
  T=1/(1 - plotpos)
  if (orient=="xF") plot(x,T, log="y", ...)
  else if (orient=="Fx") plot(T,x, log="x", ...)
  else stop("plotposRP(x, a, orient): orient unknown") 
  grid(equilogs=FALSE)
}


# ---------------------------------------------------------------- #

pointsposRP <- function(x, a=0, orient="xF", ...) {
  ordinato <- sort(x)
  n <- length(ordinato)
  i <- 1:n
  plotpos <- (i - a)/(n + 1 - 2*a)

  x=ordinato
  T=1/(1 - plotpos)
  if (orient=="xF") points(x,T,...)
  else if (orient=="Fx") points(T,x,...)
  else stop("pointsposRP(x, a, orient): orient unknown")
}




# ---------------------------------------------------------------- #

plotposRPhist <- function(xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA, col12=c(1,1), a=0, orient="xF", ...) {
 # Method in Maidment pag 18.41-42:
 #  Jery R. Stedinger, Richard M. Vogel, and Eﬁ Foufoula-Georgiou. Frequency analysis of extreme events. In David R. Maidment, editor, Hand-Book of Hydrology, chapter 18. McGraw-Hill Companies, international edition, 1993.

  # see BayesianMCMC for the notation
  # col12 = colors for systematic and historical points
 if (all(is.na(c(xcont, xhist, infhist, suphist, seuil)))) {
  stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, ...): no input data")
 }
 else if (all(is.na(c(xhist, infhist, suphist, seuil)))) {
  # systematic data only
  n <- length(xcont)
  i <- n + 1 - rank(xcont, ties.method="first")
  q <- ((i - a)/(n + 1 - 2*a))
  T <- 1/q
  x <- xcont
  if (orient=="xF") plot(x,T, log="y", col=col12[1], ...)
  else if (orient=="Fx") plot(T,x, log="x", col=col12[1], ...)
  else stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, col12, a, orient=, ...): orient unknown")
  grid(equilogs=FALSE)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
  # historical data known
  if (!all(is.na(xcont))) {
   xcont_u <- xcont[xcont >= seuil]
   xcont_l <- xcont[xcont < seuil]
  } else {
   xcont_u <- NULL
   xcont_l <- NULL
  }
  xhist_u <- xhist[xhist >= seuil]
  xhist_l <- xhist[xhist < seuil]  
  x_u <- c(xcont_u, xhist_u)
   col_x_u <- rep(col12, c(length(xcont_u), length(xhist_u)))
  x_l <- c(xcont_l, xhist_l)
   col_x_l <- rep(col12, c(length(xcont_l), length(xhist_l)))

  r <- length(x_u)                                    # number of times the threshold is exceeded
  se <- length(x_l)                               # number of times the threshold is not exceeded by station-record data
  n <- length(xcont) + nbans                          # length of the all period (historic + recent)
  qe <- r/n                                           # exceedence probability of the perception threshold

  i <- r + 1 - rank(x_u, ties.method="first")
  qi <- qe*((i - a)/(r + 1 - 2*a))                    # eq. 18.6.11 in Maidment
  T_u <- 1/qi                                          # return periods above threshold

  j <- se + 1 - rank(x_l, ties.method="first")
  qj <- qe + (1 - qe)*((j - a)/(se + 1 - 2*a))        # eq. 18.6.12 in Maidment
  T_l <- 1/qj                                     # return periods below threshold

  x <- c(x_l, x_u)
  T <- c(T_l, T_u)
  colori <- c(col_x_l, col_x_u)
  if (orient=="xF") {
   plot(x,T, log="y", col=colori, ...)
    abline(v=seuil, lty=3)
  }
  else if (orient=="Fx") {
   plot(T,x, log="x", col=colori, ...) 
   abline(h=seuil, lty=3)
  }
  else stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, col12, a, orient=, ...): orient unknown")
  grid(equilogs=FALSE)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
  # just ranges for the historicals
  xhist <- (infhist + suphist)/2
  if (!all(is.na(xcont))) {
   xcont_u <- xcont[xcont >= seuil]
   xcont_l <- xcont[xcont < seuil]
  } else {
   xcont_u <- NULL
   xcont_l <- NULL
  }
  xhist_u <- xhist[xhist >= seuil]
  xhist_l <- xhist[xhist < seuil]
  x_u <- c(xcont_u, xhist_u)
   col_x_u <- rep(c(col12[1],NA), c(length(xcont_u), length(xhist_u)))
  x_l <- c(xcont_l, xhist_l)
   col_x_l <- rep(c(col12[1],NA), c(length(xcont_l), length(xhist_l)))

  r <- length(x_u)                                    # number of times the threshold is exceeded
  se <- length(x_l)                                   # number of times the threshold is not exceeded
  n <- length(xcont) + nbans                          # length of the all period (historic + recent)
  qe <- r/n                                           # exceedence probability of the perception threshold

  i <- r + 1 - rank(x_u, ties.method="first")
  qi <- qe*((i - a)/(r + 1 - 2*a))                    # eq. 18.6.11 in Maidment
  T_u <- 1/qi                                          # return periods above threshold

  j <- se + 1 - rank(x_l, ties.method="first")
  qj <- qe + (1 - qe)*((j - a)/(se + 1 - 2*a))        # eq. 18.6.12 in Maidment
  T_l <- 1/qj                                          # return periods below threshold

  x <- c(x_l, x_u)
  T <- c(T_l, T_u)

  Thist <- sort(T[x %in% xhist])[order(xhist)]
  colori <- c(col_x_l, col_x_u)
  if (orient=="xF") {
   plot(x,T, log="y", col=colori, ...)
    abline(v=seuil, lty=3)
    points(infhist, Thist, pch=24, bg=col12[2], col=col12[2])
    points(suphist, Thist, pch=25, bg=col12[2], col=col12[2])
    segments(infhist, Thist, suphist, Thist, lty=3, col=col12[2])
  }
  else if (orient=="Fx") {
   plot(T,x, log="x", col=colori, ...)
    abline(h=seuil, lty=3)
    points(Thist, infhist, pch=24, bg=col12[2], col=col12[2])
    points(Thist, suphist, pch=25, bg=col12[2], col=col12[2])
    segments(Thist, infhist, Thist, suphist, lty=3, col=col12[2])
  }
  else stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, col12, a, orient=, ...): orient unknown")
  grid(equilogs=FALSE)
  invisible(cbind(T,x))
 }
}


# ---------------------------------------------------------------- #

pointsposRPhist <- function(xcont, xhist=NA, infhist=NA, suphist=NA, nbans=NA, seuil=NA, col12=c(1,1), a=0, orient="xF", ...) {
 # Method in Maidment pag 18.41-42:
 #  Jery R. Stedinger, Richard M. Vogel, and Eﬁ Foufoula-Georgiou. Frequency analysis of extreme events. In David R. Maidment, editor, Hand-Book of Hydrology, chapter 18. McGraw-Hill Companies, international edition, 1993.

  # see BayesianMCMC for the notation
  # col12 = colors for systematic and historical points
 if (all(is.na(c(xcont, xhist, infhist, suphist, seuil)))) {
  stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, ...): no input data")
 }
 else if (all(is.na(c(xhist, infhist, suphist, seuil)))) {
  # systematic data only
  n <- length(xcont)
  i <- n + 1 - rank(xcont, ties.method="first")
  q <- ((i - a)/(n + 1 - 2*a))
  T <- 1/q
  x <- xcont
  if (orient=="xF") points(x,T, col=col12[1], ...)
  else if (orient=="Fx") points(T,x, col=col12[1], ...)
  else stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, col12, a, orient=, ...): orient unknown")
  grid(equilogs=FALSE)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(infhist, suphist))) & all(!is.na(c(xhist, seuil, nbans)))) {
  # historical data known
  if (!all(is.na(xcont))) {
   xcont_u <- xcont[xcont >= seuil]
   xcont_l <- xcont[xcont < seuil]
  } else {
   xcont_u <- NULL
   xcont_l <- NULL
  }
  xhist_u <- xhist[xhist >= seuil]
  xhist_l <- xhist[xhist < seuil]  
  x_u <- c(xcont_u, xhist_u)
   col_x_u <- rep(col12, c(length(xcont_u), length(xhist_u)))
  x_l <- c(xcont_l, xhist_l)
   col_x_l <- rep(col12, c(length(xcont_l), length(xhist_l)))

  r <- length(x_u)                                    # number of times the threshold is exceeded
  se <- length(x_l)                                   # number of times the threshold is not exceeded
  n <- length(xcont) + nbans                          # length of the all period (historic + recent)
  qe <- r/n                                           # exceedence probability of the perception threshold

  i <- r + 1 - rank(x_u, ties.method="first")
  qi <- qe*((i - a)/(r + 1 - 2*a))                    # eq. 18.6.11 in Maidment
  T_u <- 1/qi                                          # return periods above threshold

  j <- se + 1 - rank(x_l, ties.method="first")
  qj <- qe + (1 - qe)*((j - a)/(se + 1 - 2*a))        # eq. 18.6.12 in Maidment
  T_l <- 1/qj                                          # return periods below threshold

  x <- c(x_l, x_u)
  T <- c(T_l, T_u)
  colori <- c(col_x_l, col_x_u)
  if (orient=="xF") {
   points(x,T, col=colori, ...)
    abline(v=seuil, lty=3)
  }
  else if (orient=="Fx") {
   points(T,x, col=colori, ...) 
   abline(h=seuil, lty=3)
  }
  else stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, col12, a, orient=, ...): orient unknown")
  grid(equilogs=FALSE)
  invisible(cbind(T,x))
 }
 else if (all(is.na(c(xhist))) & all(!is.na(c(infhist, suphist, seuil, nbans)))) {
  # just ranges for the historicals
  xhist <- (infhist + suphist)/2
  if (!all(is.na(xcont))) {
   xcont_u <- xcont[xcont >= seuil]
   xcont_l <- xcont[xcont < seuil]
  } else {
   xcont_u <- NULL
   xcont_l <- NULL
  }
  xhist_u <- xhist[xhist >= seuil]
  xhist_l <- xhist[xhist < seuil]
  x_u <- c(xcont_u, xhist_u)
   col_x_u <- rep(c(col12[1],NA), c(length(xcont_u), length(xhist_u)))
  x_l <- c(xcont_l, xhist_l)
   col_x_l <- rep(c(col12[1],NA), c(length(xcont_l), length(xhist_l)))

  r <- length(x_u)                                    # number of times the threshold is exceeded
  se <- length(x_l)                                   # number of times the threshold is not exceeded
  n <- length(xcont) + nbans                          # length of the all period (historic + recent)
  qe <- r/n                                           # exceedence probability of the perception threshold

  i <- r + 1 - rank(x_u, ties.method="first")
  qi <- qe*((i - a)/(r + 1 - 2*a))                    # eq. 18.6.11 in Maidment
  T_u <- 1/qi                                          # return periods above threshold

  j <- se + 1 - rank(x_l, ties.method="first")
  qj <- qe + (1 - qe)*((j - a)/(se + 1 - 2*a))        # eq. 18.6.12 in Maidment
  T_l <- 1/qj                                          # return periods below threshold

  x <- c(x_l, x_u)
  T <- c(T_l, T_u)

  Thist <- T[x %in% xhist]
  colori <- c(col_x_l, col_x_u)
  if (orient=="xF") {
   points(x,T, col=colori, ...)
    abline(v=seuil, lty=3)
    points(infhist, Thist, pch=24, bg=col12[2], col=col12[2])
    points(suphist, Thist, pch=25, bg=col12[2], col=col12[2])
    segments(infhist, Thist, suphist, Thist, lty=3, col=col12[2])
  }
  else if (orient=="Fx") {
   points(T,x, col=colori, ...)
    abline(h=seuil, lty=3)
    points(Thist, infhist, pch=24, bg=col12[2], col=col12[2])
    points(Thist, suphist, pch=25, bg=col12[2], col=col12[2])
    segments(Thist, infhist, Thist, suphist, lty=3, col=col12[2])
  }
  else stop("plotposRPhist(xcont, xhist, infhist, suphist, nbans, seuil, col12, a, orient=, ...): orient unknown")
  grid(equilogs=FALSE)
  invisible(cbind(T,x))
 }
}

