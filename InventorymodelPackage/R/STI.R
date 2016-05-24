STI <-
function(n=NA,a=NA,av=NA,d=NA,h=NA,m=NA){
if (sum(is.na(m)!=T)==length(m)|sum(is.na(d)==T)==length(d)){ 
pedido<-2*(a+av)*m/h
} else {
pedido<-sqrt(2*(a+av)*d/h)
m=d/pedido
}
coste_pedido<-2*(a+av)*m
sol<-list(pedido,coste_pedido)
names(sol)<-c("Optimal order","Order cost")
return(sol)}
