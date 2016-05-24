
####Swiss Re curves####

swissRe<-function(c)
{
  out<-numeric(2)
  b <- exp(3.1 - 0.15*c*(1+c))
  g <-exp(c*(0.78 + 0.12*c))
  out<-c(b,g)
  names(out)<-c("b","g")
  return(out)
}
