int.Simp<-function(xs,fxs){
    n<-length(xs)
    h<-xs[2]-xs[1]
    int<-3*(fxs[1]+fxs[n])/8+7*(fxs[2]+fxs[n-1])/6+23*(fxs[3]+fxs[n-2])/24+sum(fxs[4:(n-3)])
    return(int*h)
}