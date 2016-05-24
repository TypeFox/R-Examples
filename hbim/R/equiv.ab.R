`equiv.ab` <-
function(effab1,ab1,effab2,ab2,npts=100){
    abpts<- min(ab1)+ (max(ab1)-min(ab1))*(1:npts)/(npts+1)
    equiv.eff2<-approx(ab2,effab2,xout=abpts)$y
    equiv.ab1<-approx(effab1,ab1,xout=equiv.eff2)$y
    equiv.eff1<-approx(ab1,effab1,xout=abpts)$y
    out<-list(abpts=abpts,abpts10=10^abpts,equiv.eff2=equiv.eff2,equiv.ab1=equiv.ab1,
       equiv.ab110=10^equiv.ab1,equiv.eff1=equiv.eff1,x=equiv.ab1-abpts,y=equiv.eff1)
    out
}

