argmin=function(x){
    vx=as.vector(x)
    imax=order(vx)[1]
if(!is.matrix(x))imax
    else{
        d=dim(x)
        c1=as.vector(outer(seq(d[1]),rep(1,d[2])))[imax]
        c2=as.vector(outer(rep(1,d[1]),seq(d[2])))[imax]
        c(c1,c2)

    }
}
