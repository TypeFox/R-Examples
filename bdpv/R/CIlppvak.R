CIlppvak <-
function(x0, x1, p, conf.level=0.95,
 alternative=c("two.sided", "less", "greater"))
{
alternative<-match.arg(alternative)

expit<-function(p){exp(p)/(1+exp(p))}

switch(alternative,
two.sided={
z <- qnorm(p=1-(1-conf.level)/2)
seest <- setil(x1=x1, k=z)
spest <- sptil(x0=x0, k=z)
varestlppv <- varlppv(x0=x0, x1=x1, k=z)
estlppv <- logitppv(p=p, se=seest[1], sp=spest[1])
llwr<-estlppv - z*sqrt(varestlppv[1])
lupr<-estlppv + z*sqrt(varestlppv[1])
},
less={
z<-qnorm(p=conf.level)
seest <- setil(x1=x1, k=z)
spest <- sptil(x0=x0, k=z)
varestlppv <- varlppv(x0=x0, x1=x1, k=z)
estlppv <- logitppv(p=p, se=seest[1], sp=spest[1])
llwr <- (-Inf)
lupr <- estlppv + z*sqrt(varestlppv[1])
},
greater={
z<-qnorm(p=conf.level)
seest <- setil(x1=x1, k=z)
spest <- sptil(x0=x0, k=z)
varestlppv <- varlppv(x0=x0, x1=x1, k=z)
estlppv <- logitppv(p=p, se=seest[1], sp=spest[1])
llwr <- estlppv - z*sqrt(varestlppv[1])
lupr <- Inf
}
)

conf.int <- c(expit(llwr), expit(lupr))
names(conf.int)<-c("lower","upper")
estimate <- expit(estlppv)
names(estimate)<-NULL
return(list(conf.int=conf.int, estimate=estimate))

}

