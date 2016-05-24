CIlnpv <-
function(x0, x1, p, conf.level=0.95,
 alternative=c("two.sided", "less", "greater"))
{
alternative<-match.arg(alternative)

expit<-function(p){exp(p)/(1+exp(p))}

switch(alternative,
two.sided={
z<-qnorm(p=1-(1-conf.level)/2)
seest<-setil(x1=x1, k=z)
spest<-sptil(x0=x0, k=z)
varestlnpv<-varlnpv(x0=x0, x1=x1, k=z)
estlnpv <- logitnpv(p=p, se=seest[2], sp=spest[2])
llwr<-estlnpv - z*sqrt(varestlnpv[2])
lupr<-estlnpv + z*sqrt(varestlnpv[2])
},
less={
z<-qnorm(p=conf.level)
seest <- setil(x1=x1, k=z)
spest <- sptil(x0=x0, k=z)
varestlnpv <- varlnpv(x0=x0, x1=x1, k=z)
estlnpv <- logitnpv(p=p, se=seest[2], sp=spest[2])
llwr <- (-Inf)
lupr <- estlnpv + z*sqrt(varestlnpv[2])
},
greater={
z<-qnorm(p=conf.level)
seest <- setil(x1=x1, k=z)
spest <- sptil(x0=x0, k=z)
varestlnpv <- varlnpv(x0=x0, x1=x1, k=z)
estlnpv <- logitnpv(p=p, se=seest[2], sp=spest[2])
llwr <- estlnpv - z*sqrt(varestlnpv[2])
lupr <- Inf
}
)

conf.int <- c(expit(llwr), expit(lupr))
names(conf.int)<-c("lower","upper")
estimate <- expit(estlnpv)
names(estimate)<-NULL
return(list(conf.int=conf.int, estimate=estimate))
}

