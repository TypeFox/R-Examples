likEvid<-function(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq)
{	
	sortieR<-0
	# code to change freq in case of rare alleles
	Repliste2<-Repliste[which(Repliste!=0)]#remove zeros from Repliste
	rare<-Repliste2[which(!Repliste2 %in% names(freq),arr.ind=TRUE)]
	if(length(rare)!=0){ 
	freq[as.character(rare)]<-1/(2*2085)
	# rare<-which(Repliste2[[j]] %in% in freq)}
	# warning('allele',rare,'has been added to the database frequency')
	print(paste('WARNING: allele',rare,'has been added to the database with frequency',1/(2*2085),sep=' '))

	}

	appelC<-function(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq,sortieR)
	{

		lenRepliste<-length(Repliste)
		lenT<-length(T)
		lenV<-length(V)
		lenHom<-length(prDHom)
		lenHet<-length(prDHet)
		allele<-names(freq)
		freq0<-as.numeric(freq)
		lenFreq<-length(freq)
		sortieR<-0
		.C('evidenceC3', as.double(Repliste), as.integer(lenRepliste),as.double(T), as.integer(lenT), as.double(V), as.integer(lenV),
		as.integer(x), as.double(theta), as.double(prDHet),as.integer(lenHet),as.double(prDHom),as.integer(lenHom),as.double(prC),as.character(allele),as.double(freq0),as.integer(lenFreq), sortieR=as.double(sortieR),PACKAGE='forensim')
		
		
	}
	
	tmp=(appelC(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq,sortieR))
	return(tmp$sortieR)
}


LR<-function(Repliste,Tp,Td,Vp,Vd,xp,xd,theta,prDHet,prDHom,prC,freq){

num<-likEvid(Repliste,T=Tp,V=Vp,x=xp,theta,prDHet,prDHom,prC,freq)
deno<-likEvid(Repliste,T=Td,V=Vd,x=xd,theta,prDHet,prDHom,prC,freq)


list('num'=num,'deno'=deno,'LR'=num/deno)}