eva.prod<-function(x,marginal="student",g=1)
{
if (marginal=="student"){
    d<-1
    vakio<-gamma((g+d)/2)/((g*pi)^(d/2)*gamma(g/2))
    y<-vakio*(1+x^2/g)^(-(g+d)/2)
    val<-prod(y)
}
if (marginal=="student.old"){
    d<-1
    vakio<-gamma((g+d)/2)/((g*pi)^(d/2)*gamma(g/2))
    x1<-vakio*(1+x[1]^2/g)^(-(g+d)/2)
    x2<-vakio*(1+x[2]^2/g)^(-(g+d)/2)
    val<-x1*x2
}
if (marginal=="studentR"){
    #x1<-dt(x[1],df=g)
    #x2<-dt(x[2],df=g)
    y<-dt(x,df=g)
    val<-prod(y)  
}
if (marginal=="polyno.old"){
    vakio<-2*(1-1/(g+1))
    y<-vakio*abs(1-x)^g
    val<-prod(y)
}
if (marginal=="polyno"){ 
    vakio<-1/(2*(1-1/(g+1)))
    y<-vakio*abs(1-abs(x)^g)
    val<-prod(y)
}
if (marginal=="double"){
    vakio<-1/2
    y<-exp(-abs(x))
    val<-prod(y)
}
if (marginal=="gauss"){
    vakio<-(2*pi)^(-1/2)
    y<-exp(-x^2/2)
    val<-prod(y)
}

  
return(val)
}

