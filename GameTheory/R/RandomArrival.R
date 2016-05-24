RandomArrival <-
function(E,C,Names=NULL){

z<-length(C)
P<-permn(z)	
elementos<-length(P)
pagos<-matrix(0,nrow=elementos,ncol=z)
for (i in 1:elementos){
	orden<-P[[i]]
	for (j in 1:z){
		pagos[i, orden[j] ]<- min(C[orden[j]],max(0,E-sum(pagos[i,])))
		}
	}
	
	
	res<-apply(pagos,2,mean)
	R<-as.matrix(res)
 	rownames(R)<-Names
 	colnames(R)<-"RA"
 	Output<-list(Results=R,Claims=C,Method="Random Arrival",Short="RA",E=E,Names=Names)
    class(Output)<-"ClaimsRule"
    return(Output)  
	
}
