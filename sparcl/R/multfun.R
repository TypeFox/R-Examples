multfun=function(x){
#if(!is.loaded("multfun")){
#  dyn.load("multfun.so")
#}
n<-nrow(x)
p <- ncol(x)
x[is.na(x)]=0
mode(x)="single"
n2=n*(n-1)/2
junk=.Fortran("multfun",
         x,
        as.integer(n),
       as.integer(p),
       as.integer(n2),
       d=single(n2*p), PACKAGE="sparcl"
)
return(junk$d)
}

