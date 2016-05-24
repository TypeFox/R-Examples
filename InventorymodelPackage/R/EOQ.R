EOQ <-
function(n=NA,a=NA,d=NA,h=NA,m=NA){
if (sum(is.na(m)!=T)==length(m)|sum(is.na(d)==T)==length(d)){ 
#caso demanda "d" desconocida
pedido<-2*a*m/h
} else {
pedido<-sqrt(2*a*d/h)
m=d/pedido
}
coste_pedido<-2*a*m
sol<-list(pedido,coste_pedido)
names(sol)<-c("Optimal order","Order costs")
return(sol)}
