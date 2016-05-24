distfun=function(x){
#if(!is.loaded("distfun")){
#  dyn.load("distfun.so")
#}
n<-nrow(x)
p <- ncol(x)
x[is.na(x)]=0
mode(x)="single"
n2=n*(n-1)/2
junk=.Fortran("distfun",
         x,
        as.integer(n),
       as.integer(p),
       as.integer(n2),
       d=single(n2*p), PACKAGE="sparcl"
)
return(junk$d)
}

