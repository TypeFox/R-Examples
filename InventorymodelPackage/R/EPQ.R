EPQ <-
function(n=NA,a=NA,d=NA,h=NA,m=NA,r=NA,s=NA){
if (sum(is.na(m)==T)!=length(m)|sum(is.na(d)==T)==length(d)){ 
pedido<-2*a*m/h
d=m*pedido
} else {
pedido<-sqrt(2*a*d*(h+s)/(h*(1-d/r)*s))
m=d/pedido
}
faltantes<-sqrt(2*a*d*h*(1-d/r)/(s*(h+s)))
coste_pedido<-2*a*m
sol<-list(pedido,faltantes,coste_pedido)
names(sol)<-c("Optimal order","Optimal shortages","Order costs")
return(sol)}
