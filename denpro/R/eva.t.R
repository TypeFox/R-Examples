eva.t<-function(x,df,mu,Sigma)
{
norma<-gamma((df+1)/2)/((df*pi)^(1/2)*gamma(df/2))
fu<-(1+x^2/df)^(-(df+1)/2)
val<-norma*fu

return(val)
}


