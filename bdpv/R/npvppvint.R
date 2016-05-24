npv <-
function(p, se, sp){sp*(1-p)/((1-se)*p + (sp)*(1-p))}

logitnpv <-
function(p, se, sp){log(sp*(1-p)/((1-se)*p))}

ppv <-
function(p, se, sp){se*p/(se*p + (1-sp)*(1-p))}


logitppv <-
function(p, se, sp){log(se*p/((1-sp)*(1-p)))}

setil <-
function(x1, k){
n1 <- sum(x1)
sehat <- x1[1]/n1
setil <- ( n1*sehat + k^(2)/2 )/(n1+k^(2))
return(c(setil=setil, sehat=sehat))
}

sptil <-
function(x0, k){
n0 <- sum(x0)
sphat <- x0[2]/n0
sptil <- ( n0*sphat + k^(2)/2 )/(n0+k^(2))
return(c(sptil=sptil, sphat=sphat))
}

varlnpv <-
function(x0, x1, k){
n1<-sum(x1)
n0<-sum(x0)
seest<-setil(x1=x1, k=k)
spest<-sptil(x0=x0, k=k)
vartil <- (seest[1]/(1-seest[1]))*(1/(n1 + k^(2))) + ((1-spest[1])/spest[1])*(1/(n0 + k^(2)))
varhat <- (seest[2]/(1-seest[2]))*(1/n1) + ((1-spest[2])/spest[2])*(1/n0)
return(c(vartil=vartil, varhat=varhat))
}

varlppv <-
function(x0, x1, k, p){
n1<-sum(x1)
n0<-sum(x0)
seest<-setil(x1=x1, k=k)
spest<-sptil(x0=x0, k=k)
vartil <- ((1-seest[1])/seest[1])*(1/(n1+k^(2))) + (spest[1]/(1-spest[1]))*(1/(n0+k^(2)))
varhat <- ((1-seest[2])/seest[2])*(1/n1) + (spest[2]/(1-spest[2]))*(1/n0)
return(c(vartil=vartil, varhat=varhat))
}


