#==============================================================
# agosto 05
# se actualiza en marzo 7/07
# ordenamientos de inercia 
# modificacion septiembre 8
# entra parti = la lista de salida de la funcion partial.wwm
#    coro="col" "row" escoje inercia para filas o columnas
#    pato = "par" "tot" escoje inercia para puntos parciales
#          o la inercia de la nube de un punto
#   can = número de filas mayores y menores en la salida
#   dec = numero de decimales
#         ax = 0, para el subespacio de nf dimension 
# la inercia intra se multiplica x10000
#==============================================================
wwinertia.ord <- function(parti,ax=0,coro="row",pato="tot",can=5,dec=1)
{
# control de entrada
  if (!inherits(parti, "parwwm")) 
        stop("non convenient data")
names(parti$cw) <- rownames(parti$col.wit)
cat("\n Within inertias are multiplicate by 10000 \n")
if(ax != 0) 
  {
 	if (coro=="row" & pato=="tot"){
   		veciner <- parti$row.wit[,ax]
   		w <- parti$lw
   	}
 	if (coro=="col" & pato=="tot" ){
   		veciner <- parti$col.wit[,ax]
   		w <- parti$cw
   	}
 	if (coro=="row" & pato=="par" ){
   		veciner <- parti$row.cwit[,ax]
   		w <- NULL
   	}
 	if (coro=="col" & pato=="par" ){
   		veciner <- parti$col.cwit[,ax]
   		w <- NULL
   	}
  }
if(ax==0)
   {   
 	cat("\n Subspace dimension ", parti$nf,"\n")
 	if (coro=="row" & pato=="tot" ){
   		veciner <- parti$row.witS
   		w <- parti$lw
   	}
 	if (coro=="col" & pato=="tot"){
   		veciner <- parti$col.witS
   		w <- parti$cw
   	}
	if (coro=="row" & pato=="par"){
   		veciner <- parti$row.cwitS
		w <- NULL
   		w[1:length(parti$lw)] <- 0
   	}
	if (coro=="col" & pato=="par"){
   		veciner <- parti$col.cwitS
   		cat("\n Subspace dimension ", parti$nf,"\n")
   		w <-NULL
		w[1:length(parti$cw)] <- 0
   	}
   }
first <- sort(veciner,decreasing=TRUE)[1:can]
wf <-w[names(first)]
last <- sort(veciner,decreasing=FALSE)[1:can]
cumlast <- cumsum(last)
total <- sum(veciner)
last <- sort(last,decreasing=TRUE)
wl <- w[names(last)]
cumlast <- sort(cumlast,decreasing=TRUE)
Iord <- c(first,last)
cumfirst <- cumsum(first)
Iord <- cbind(c(1:can,can:1),Iord*10000,Iord*100/total,c(cumfirst,cumlast)*100/total)
colnames(Iord) <- c("No","WitIner","PercWit","PercCum")
if(pato=="tot") {
    Iord <- cbind(Iord,c(wf,wl)*100)
    colnames(Iord) <- c("No","WitIner","PercWit","PercCum","Weight")
}
Iord <- round(Iord,dec)
return(Iord)
}
#==============================================================

