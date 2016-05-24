mwhc2c<-
function(n=NA,a=NA,b=NA,d=NA,K=NA,c1=NA,c2=NA,cooperation=c(0,1),allocation=c(0,1)){

if (sum(sort(d/K)==(d/K))!=n) {
print("Warning: agents, d and K are not in the order indicated by the ratios d/K.")
d0K0<-d/K
d<-d[order(d0K0)]
K<-K[order(d0K0)]
}
if (cooperation==0){
	pedido<-sqrt(b*d/(2*a+2*(c1-c2)*K+b*K^2/d))
	costes<-sqrt(b*d*(2*a+2*(c1-c2)*K+b*K^2/d))-b*K+c2*d
	faltantes<-d/pedido-K
	sol<-list(pedido,costes,faltantes)
	names(sol)<-c("Number of orders per time unit","Costs","Optimal shortages")
}
if (cooperation==1){
matriz0<-coalitions(n)
matriz<-cbind(as.vector(matriz0[[2]]),rep(0,2^n),rep(0,2^n))
matriz0<-coalitions(n)[[1]]
matriz1<-matrix(0,ncol=n+1,nrow=2^n)
for (i in 2:nrow(matriz)){
	matriz0[i,]->coa
	aux<-which(coa==1)
	coaux<-aux
	dcoa<-d[aux];Kcoa<-K[aux];bcoa<-b[aux]
	s=length(aux)
	T<-rep(0,length(coa))
	coste<-Inf;pedido<-0

	for (l in s:1){
		k=l;T[coaux[k]]=1;T
		aux1<-which(T==1)
		xT<-sqrt(sum(b[aux1]*d[aux1])/(2*a+2*(c1-c2)*sum(K[aux1])+sum(b[aux1]*K[aux1]^2/d[aux1])))
		SxT<-which(xT<dcoa/Kcoa);SxT2<-which(xT>=dcoa/Kcoa)
		xTaux<-xT

		if (length(SxT)>length(aux1)){xTaux<-dcoa[SxT[1]]/Kcoa[SxT[1]]}
		if (length(SxT)<length(aux1)){xTaux<-dcoa[which(SxT2>=l)[1]]/Kcoa[which(SxT2>=l)[1]]}
		
		coste2<-xTaux*(a+1/xTaux*c1*sum(dcoa[SxT2])+c1*sum(Kcoa[SxT]))+
		1/xTaux*sum(bcoa[SxT]/(2*dcoa[SxT])*(-Kcoa[SxT]*xTaux+dcoa[SxT])^2)+
		sum(c2*(dcoa[SxT]-Kcoa[SxT]*xTaux))
		if (coste>coste2){pedido<-xTaux;coste<-coste2}
	}
	matriz[i,2]<-pedido
	matriz[i,3]<-coste
	if (allocation==1){
		SxT<-which(pedido<dcoa/Kcoa)
		suma1<-0
		for (k in 1:length(SxT)){suma1<-suma1+aux[SxT[k]]*10^(length(SxT)-k)}
		matriz1[i,aux[SxT]+1]<-(bcoa[SxT]*dcoa[SxT])/2*pedido*(2*a+2*(c1-c2)*sum(Kcoa[SxT])+sum(bcoa[SxT]*Kcoa[SxT]^2/dcoa[SxT]))/sum(bcoa[SxT]*dcoa[SxT])+1/pedido*(bcoa[SxT]*dcoa[SxT])/2+c2*dcoa[SxT]-bcoa[SxT]*Kcoa[SxT]
		matriz1[i,aux[-SxT]+1]<-c1*dcoa[-SxT]
		matriz1[i,1]<-suma1
		
	}
}
colnames(matriz)<-c("Coalitions","Optimal orders","Costs")
rownames(matriz)<-rep(" ",2^n)
sol<-matriz
if (allocation==1){
	colnames(matriz1)<-c("Coalition_SxT",1:n)
	rownames(matriz1)<-rep(" ",2^n)	
	sol<-list(matriz,matriz1);names(sol)<-c("Optimal policies","GR-rule")
}
}
return(sol)}
