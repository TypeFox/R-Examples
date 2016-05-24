
# Score test according to Koopman, as described by Gart and Nam (1988)

# ML estimate under restriction ratio=psi
# Gart&Nam,1988, eq.(3.3)

MLEScore<-function(psi, x0, x1, n0, n1){
va<-(n0+n1)*psi
vb<-(-((x0+n1)*psi + x1 + n0))
vc<-x0+x1
p0<-(-vb + c(-1,1)*sqrt(vb^2 - 4*va*vc))/(2*va)
p1<-p0*psi
return(c(p0=p0[1], p1=p1[1]))
}


varpiScore <- function(p0, p1, n0, n1){(1-p0)/(n0*p0) + (1-p1)/(n1*p1)}

# chi^2_K, page 326, Gart&Nam (1988)

QEscore <- function(psi, x0, x1, n0, n1, quantile, digits){
MLE <- MLEScore(psi=psi, x0=x0, x1=x1, n0=n0, n1=n1)
pt0 <- round(MLE["p0"], digits=digits)
pt1 <- pt0*psi
vpsipt0 <- 1/varpiScore(p0=pt0, p1=pt1, n0=n0, n1=n1)
return( ((x1-n1*pt1)^2 ) / ( ((1-pt1)^2) * vpsipt0 ) - quantile^2)
}


# Score test and intervals acc. to
# Gart&Nam, 1988, Section 3.3

RRscoreci <- function(x0, x1, n0, n1, conf.level=0.95, alternative="two.sided", loglow=-20, logupp=20, digits=12)
{


METHOD<-"Gart-Nam Score interval"

limitscore <- function(x0, x1, n0, n1, quantile, interval, digits)
{lwr<-uniroot(QEscore, interval=interval, x0=x0, x1=x1, n0=n0, n1=n1, quantile=quantile, digits=digits, tol=.Machine$double.eps^0.5); return(lwr$root)}

switch(alternative,
"two.sided"={QUANT<-qnorm(p=1-(1-conf.level)/2)},
"less"={QUANT<-qnorm(p=1-(1-conf.level))},
"greater"={QUANT<-qnorm(p=1-(1-conf.level))})

# distinguish the following special cases:

# c)
if(x0==0 & x1==0){estimate<-NaN; lower<-0; upper<-Inf}

# b1)
if(x0==n0 & x1==0){estimate<-0; lower<-0
switch(alternative,
"two.sided"={upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"less"={upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"greater"={upper<-Inf})
}

# b2)
if(0<x0 & x0<n0 & x1==0){estimate<-0; lower<-0
switch(alternative,
"two.sided"={upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"less"={upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"greater"={upper<-Inf})
}


# a1)
if(x0==0 & x1==n1){estimate<-Inf; upper<-Inf
switch(alternative,
"two.sided"={upperstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar},
"less"={lower<-0},
"greater"={upperstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar})
}

# a2)
if(x0==0 & 0<x1 & x1<n1){estimate<-Inf; upper<-Inf
switch(alternative,
"two.sided"={upperstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar},
"less"={lower<-0},
"greater"={upperstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar})
}

# e)
if(x0==n0 & x1==n1){estimate<-1; fold<-1+1/(n0+n1)
switch(alternative,
"two.sided"={
lower<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), 1/fold), digits=digits)
#upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(1*fold, exp(logupp)), digits=digits)
lowerstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), 1/fold), digits=digits)
upper<-1/lowerstar
},
"less"={
lower<-0
#upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(1*fold, exp(logupp)), digits=digits)
lowerstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), 1/fold), digits=digits)
upper<-1/lowerstar
},
"greater"={
lower<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), 1/fold), digits=digits)
upper<-Inf
})
}

# d1)
if(x0==n0 & 0<x1 & x1<n1){ estimate <- (x1/n1)/(x0/n0)
switch(alternative,
"two.sided"={
lower<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"less"={
lower<-0
upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"greater"={
lower<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-Inf
})
}


# d2)
if(0<x0 & x0<n0 & x1==n1){ estimate <- (x1/n1)/(x0/n0)
switch(alternative,
"two.sided"={
lowerstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), 1/estimate), digits=digits)
upperstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(1/estimate, exp(logupp)), digits=digits)
lower<-1/upperstar
upper<-1/lowerstar
},
"less"={
lowerstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), 1/estimate), digits=digits)
lower<-0
upper<-1/lowerstar
},
"greater"={
upperstar<-limitscore(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(1/estimate, exp(logupp)), digits=digits)
lower<-1/upperstar
upper<-Inf
}) 
}

if(0<x0 & x0<n0 & 0<x1 & x1<n1){ estimate <- (x1/n1)/(x0/n0)
switch(alternative,
"two.sided"={
lower<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"less"={
lower<-0
upper<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"greater"={
lower<-limitscore(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-Inf
})
}


conf.int<-c(lower=lower, upper=upper)
attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate,
conf.level=conf.level,
quantile=QUANT
)  
) 
}



#########################################################################

#  Miettinen Nurminen
# chi^2_MN, page 326, Gart&Nam (1988)

QEmn <- function(psi, x0, x1, n0, n1, quantile, digits){
p0h<-x0/n0
p1h<-x1/n1
n<-n0+n1
MLE <- MLEScore(psi=psi, x0=x0, x1=x1, n0=n0, n1=n1)
p0mle <- round(MLE["p0"], digits=digits)
p1mle <- p0mle*psi
return(((p1h-psi*p0h)^2)/(p1mle*(1-p1mle)/n1 + ((psi^2)*p0mle*(1-p0mle)/n0))*((n-1)/n) - quantile^2)
}


RRMNscoreci <- function(x0, x1, n0, n1, conf.level=0.95, alternative="two.sided", loglow=-20, logupp=20, digits=12)
{
METHOD<-"Miettinen-Nurminen Score interval"

limitmn <- function(x0, x1, n0, n1, quantile, interval, digits)
{lwr<-uniroot(QEmn, interval=interval, x0=x0, x1=x1, n0=n0, n1=n1, quantile=quantile, digits=digits, tol=.Machine$double.eps^0.5); return(lwr$root)}

switch(alternative,
"two.sided"={QUANT<-qnorm(p=1-(1-conf.level)/2)},
"less"={QUANT<-qnorm(p=1-(1-conf.level))},
"greater"={QUANT<-qnorm(p=1-(1-conf.level))})

# distinguish the following special cases:

# c)
if(x0==0 & x1==0){estimate<-NaN; lower<-0; upper<-Inf}

# b1)
if(x0==n0 & x1==0){estimate<-0; lower<-0
switch(alternative,
"two.sided"={upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"less"={upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"greater"={upper<-Inf})
}

# b2)
if(0<x0 & x0<n0 & x1==0){estimate<-0; lower<-0
switch(alternative,
"two.sided"={upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"less"={upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)},
"greater"={upper<-Inf})
}


# a1)
if(x0==0 & x1==n1){estimate<-Inf; upper<-Inf
switch(alternative,
"two.sided"={upperstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar},
"less"={lower<-0},
"greater"={upperstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar})
}

# a2)
if(x0==0 & 0<x1 & x1<n1){estimate<-Inf; upper<-Inf
switch(alternative,
"two.sided"={upperstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar},
"less"={lower<-0},
"greater"={upperstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), exp(logupp)), digits=digits)
 lower<-1/upperstar})
}

# e)
if(x0==n0 & x1==n1){estimate<-1; fold<-1+1/(n0+n1)
switch(alternative,
"two.sided"={
lower<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), 1/fold), digits=digits)
upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(1*fold, exp(logupp)), digits=digits)
},
"less"={
lower<-0
upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(1*fold, exp(logupp)), digits=digits)
},
"greater"={
lower<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), 1/fold), digits=digits)
upper<-Inf
})
}

# d1)
if(x0==n0 & 0<x1 & x1<n1){ estimate <- (x1/n1)/(x0/n0)
switch(alternative,
"two.sided"={
lower<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"less"={
lower<-0
upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"greater"={
lower<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-Inf
})
}


# d2)
if(0<x0 & x0<n0 & x1==n1){ estimate <- (x1/n1)/(x0/n0)
switch(alternative,
"two.sided"={
lowerstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), 1/estimate), digits=digits)
upperstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(1/estimate, exp(logupp)), digits=digits)
lower<-1/upperstar
upper<-1/lowerstar
},
"less"={
lowerstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(exp(loglow), 1/estimate), digits=digits)
lower<-0
upper<-1/lowerstar
},
"greater"={
upperstar<-limitmn(x0=x1, x1=x0, n0=n1, n1=n0, quantile=QUANT, interval=c(1/estimate, exp(logupp)), digits=digits)
lower<-1/upperstar
upper<-Inf
}) 
}

if(0<x0 & x0<n0 & 0<x1 & x1<n1){ estimate <- (x1/n1)/(x0/n0)
switch(alternative,
"two.sided"={
lower<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"less"={
lower<-0
upper<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(estimate, exp(logupp)), digits=digits)
},
"greater"={
lower<-limitmn(x0=x0, x1=x1, n0=n0, n1=n1, quantile=QUANT, interval=c(exp(loglow), estimate), digits=digits)
upper<-Inf
})
}


conf.int<-c(lower=lower, upper=upper)
attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate,
conf.level=conf.level,
quantile=QUANT
)  
) 


}



RRcrudenormalci<-function(x0, x1, n0, n1, conf.level=0.95, alternative="two.sided")
{
METHOD<-"Crude log interval"

x0I<-x0+0.5; x1I<-x1+0.5
n0I<-n0+0.5; n1I<-n1+0.5

estI <- log( (x1I/n1I)/(x0I/n0I) )
stderrlog <- sqrt( 1/x1I + 1/x0I - 1/n1I - 1/n0I )

estimate <- (x1/n1)/(x0/n0)

switch(alternative,
"two.sided"=={
   QUANT <- qnorm(p = 1-(1-conf.level)/2 ) 
   lower <- estI - QUANT * stderrlog 
   upper <- estI + QUANT * stderrlog
  },
"less"={
   QUANT <- qnorm(p = conf.level ) 
   lower <- (-Inf)
   upper <- estI + QUANT * stderrlog
  },
"greater"={ 
   QUANT <- qnorm(p = conf.level )
   lower <- estI - QUANT * stderrlog
   upper <- Inf
  })
if(is.na(lower)){lower <- -Inf}
if(is.na(upper)){upper <- Inf}

conf.int<-exp(c(lower=lower, upper=upper))
attr(conf.int, which="methodname")<-METHOD

return(list(conf.int=conf.int,
estimate=estimate,
conf.level=conf.level,
quantile=QUANT
))
}





RRmoverrci<-function(x0, x1, n0, n1, conf.level=0.95, alternative="two.sided")
{
METHOD<-"Method of variance estimates recovery (Donner, Zou, 2012)"
# Wilson(1927) CI for one binomial proportion (Agresti, Coull, 1998)

Wilsonci<-function(x, n, quant){pe<-x/n; (pe + (quant^2)/(2*n) + c(-1,1)*quant*sqrt((pe*(1-pe) + (quant^2)/(4*n))/n))/(1+(quant^2)/n)}

p1<-(x1/n1); p0<-(x0/n0)

estimate<-(x1/n1)/(x0/n0)

# Eq. (9), Donner and Zou, 2012, Stat Methods Med Res 2012 21: 347-359.

switch(alternative,
"two.sided"={
QUANT <- qnorm(1-(1-conf.level)/2)
wilci1 <- Wilsonci(x=x1, n=n1, quant=QUANT)
wilci0 <- Wilsonci(x=x0, n=n0, quant=QUANT)
lower <- (p1*p0 - sqrt( (p1*p0)^2 - wilci1[1]*wilci0[2]*(2*p1 - wilci1[1])*(2*p0-wilci0[2])))/(wilci0[2]*(2*p0 - wilci0[2]))
if(x0==0){ upper<-Inf}else{
upper <- (p1*p0 + sqrt( (p1*p0)^2 - wilci1[2]*wilci0[1]*(2*p1 - wilci1[2])*(2*p0-wilci0[1])))/(wilci0[1]*(2*p0 - wilci0[1]))}
},

"greater"={
QUANT <- qnorm(1-(1-conf.level))
wilci1 <- Wilsonci(x=x1, n=n1, quant=QUANT)
wilci0 <- Wilsonci(x=x0, n=n0, quant=QUANT)
if(x1==0){ lower <- 0}else{
lower <- (p1*p0 - sqrt( (p1*p0)^2 - wilci1[1]*wilci0[2]*(2*p1 - wilci1[1])*(2*p0-wilci0[2])))/(wilci0[2]*(2*p0 - wilci0[2]))}
upper <- Inf
},

"less"={
QUANT <- qnorm(1-(1-conf.level))
wilci1 <- Wilsonci(x=x1, n=n1, quant=QUANT)
wilci0 <- Wilsonci(x=x0, n=n0, quant=QUANT)
lower <- 0
if(x0==0){ upper <- Inf}else{
upper <- (p1*p0 + sqrt( (p1*p0)^2 - wilci1[2]*wilci0[1]*(2*p1 - wilci1[2])*(2*p0-wilci0[1])))/(wilci0[1]*(2*p0 - wilci0[1]))}
})

conf.int<-c(lower=lower, upper=upper)
attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate,
conf.level=conf.level,
quantile=QUANT
)  
) 
}







####################################


"Prop.ratio" <-
function(x, y, conf.level=0.95, alternative="two.sided", CImethod=c("Score","MNScore","MOVER","GNC"))
{

CImethod<-match.arg(CImethod)
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

 if( is.data.frame(x) & is.data.frame(y) )
  {
   colsx<-colSums(x)
   colsy<-colSums(y)
   n1<-sum(colsx)
   n0<-sum(colsy)
   x1<-colsx[1]
   x0<-colsy[1]
  }
  else
   {
    if(!((is.numeric(x)|is.integer(x)) & (is.numeric(y)|is.integer(y))) ){stop("x, y must be integer (or numeric) vectors!")}  
    if(length(x)!=2 | length(y)!=2 ){stop("x, y should be vectors of length 2!")}
      n1<-sum(x)
      n0<-sum(y)
      x1<-x[1]
      x0<-y[1]
   }


if(any(c(x0<0, x1<0, n0<0, n1<0 ))){stop("x, y must be positive!")}

# CHECK for non-integer input
if(any( abs(c( x0-as.integer(x0), x1-as.integer(x1), n0-as.integer(n0), n1-as.integer(n1))) > sqrt(.Machine$double.eps) ) ){warning("At least one observed count in x, y is no integer value.")}

# Warning for small sample case:
if(any(c(x0<5, x1<5,(n0-x0)<5, (n1-x1)<5))){warning("Chi-square (or normal) approximation might be incorrect (at least one cell count is below 5).")}

switch(CImethod,
"GNC"={
OUT<-RRcrudenormalci(x0=x0, x1=x1, n0=n0, n1=n1, conf.level=conf.level, alternative=alternative)
},


"Score"={
OUT<-RRscoreci(x0=x0, x1=x1, n0=n0, n1=n1, conf.level=conf.level, alternative=alternative, loglow=-20, logupp=20)
},


"MNScore"={
OUT<-RRMNscoreci(x0=x0, x1=x1, n0=n0, n1=n1, conf.level=conf.level, alternative=alternative, loglow=-20, logupp=20)
},


"MOVER"={
OUT<-RRmoverrci(x0=x0, x1=x1, n0=n0, n1=n1, conf.level=conf.level, alternative=alternative)
})

names(OUT$estimate)<-paste("ratio of proportions ((",x1,"/",n1,")/(",x0,"/",n0,"))",sep="")

return(OUT) 
}


