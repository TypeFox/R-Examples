data.generation <-
function(gmt.file.raw,dataset,phenotype,gene.set.size)
{
	#====  UnitTests for the function bic.generation
	num.gene.data=sum(apply(is.finite(dataset),1,sum))
	if(num.gene.data!=dim(dataset)[1]*dim(dataset)[2])
		stop("Some entries in the dataset matrix are not numerical")	
	num.feno=phenotype%in%c(1:dim(dataset)[2])
	if(sum(num.feno)!=length(phenotype))
		stop("Some entries of the phenotype identifiers do not correspond to the dataset matrix")
	if(gene.set.size<2)
		stop("Gene set most have at least 2 genes")
	nom.genes=rownames(dataset)
	if(is.character(nom.genes)==F)
		stop("Dataset most have gene symbols as row identifiers")	
	#=======================================================
	
	datos=dataset[,phenotype]
	genes.diff.uncs.bns=rownames(datos)
	db.proc=readLines(gmt.file.raw)
	db.proc=lapply(db.proc,function(x){unlist(strsplit(x,"\t"))})
	db.proc.gns.diff=list()
	nom.proc=0
	for(i in 1:length(db.proc))
   	{
   		aux=match(db.proc[[i]][-(1:2)],genes.diff.uncs.bns,nomatch=0)
   		db.proc.gns.diff[[i]]=db.proc[[i]][-(1:2)][which(aux!=0)]
   		nom.proc=c(nom.proc,db.proc[[i]][1])
   	}
	nom.proc=nom.proc[-1]
	names(db.proc.gns.diff)=nom.proc

	orden.nom.proc=sort(nom.proc,index.return=T)
	db.proc.gns.diff=db.proc.gns.diff[orden.nom.proc$ix]
	nom.proc=orden.nom.proc$x
	tam.procs=unlist(lapply(db.proc.gns.diff,length))
	names(tam.procs)=NULL
	cls.proc=which(tam.procs>=gene.set.size)
	db.proc.gns.diff.cls=db.proc.gns.diff[cls.proc]
	gns.en.proc=sort(unique(unlist(db.proc.gns.diff.cls)))
	gns.dnd.proc=list()
	for(i in 1:length(gns.en.proc))
   	{
   		aux=lapply(db.proc.gns.diff.cls,function(x){ match(gns.en.proc[i],x,nomatch=0)})
    	gns.dnd.proc[[i]]=which(aux!=0)
    }
	names(gns.dnd.proc)=gns.en.proc
	donde.genes=match(gns.en.proc,genes.diff.uncs.bns)
	donde.affy=donde.genes
	db.proc.ids.diff.cls=lapply(db.proc.gns.diff.cls,function(x){ match(x,gns.en.proc) })
	cls.gns=sort(unique(unlist(db.proc.ids.diff.cls)))
	datos.fc=lapply(db.proc.ids.diff.cls,function(x){
                               colMeans(datos[donde.affy[x],],na.rm=T)/apply(datos[donde.affy[x],],2,function(x){
                              	      sd(x,na.rm=T)
                              })})                                                
	datos.fc=do.call(rbind,datos.fc)
	if(length(datos.fc[,1])==length(unique(datos.fc[,1])))
	{
		db.proc.ids.diff.cls.2=db.proc.ids.diff.cls
		tam.ids.proc=lapply(db.proc.ids.diff.cls.2,length)
		tam.ids.proc=unlist(tam.ids.proc)
		names(tam.ids.proc)=NULL
		nrow.proc.ids=length(db.proc.ids.diff.cls.2)
		ncol.proc.ids=max(tam.ids.proc)
		db.proc.ids.diff.cls.mat.2=matrix(0,nrow=nrow.proc.ids,ncol=ncol.proc.ids)
		for(i in 1:nrow.proc.ids)
   			db.proc.ids.diff.cls.mat.2[i,1:tam.ids.proc[i]]=db.proc.ids.diff.cls.2[[i]]
		db.proc.ids.diff.cls.mat.2=cbind(db.proc.ids.diff.cls.mat.2,tam.ids.proc)
		
		gns.dnd.proc.2=gns.dnd.proc	           
		tam.ids.gns=lapply(gns.dnd.proc.2,length)
		tam.ids.gns=unlist(tam.ids.gns)
		names(tam.ids.gns)=NULL
		nrow.gns.ids=length(gns.dnd.proc)
		ncol.gns.ids=max(tam.ids.gns)
		gns.dnd.proc.mat.2=matrix(0,nrow=nrow.gns.ids,ncol=ncol.gns.ids)
		for(i in 1:nrow.gns.ids)
   			gns.dnd.proc.mat.2[i,1:tam.ids.gns[i]]=gns.dnd.proc.2[[i]]
		gns.dnd.proc.mat.2=cbind(gns.dnd.proc.mat.2,tam.ids.gns)
	}else{
		no.estan=match(unique(datos.fc[,1]),datos.fc[,1])
		aux=match(seq(1,length(datos.fc[,1])),sort(no.estan),nomatch=0)
		no.estan=which(aux==0)
		datos.fc=datos.fc[-no.estan,]

		db.proc.ids.diff.cls.2=db.proc.ids.diff.cls[-no.estan]
		tam.ids.proc=lapply(db.proc.ids.diff.cls.2,length)
		tam.ids.proc=unlist(tam.ids.proc)
		names(tam.ids.proc)=NULL
		nrow.proc.ids=length(db.proc.ids.diff.cls.2)
		ncol.proc.ids=max(tam.ids.proc)
		db.proc.ids.diff.cls.mat.2=matrix(0,nrow=nrow.proc.ids,ncol=ncol.proc.ids)
		for(i in 1:nrow.proc.ids)
   			db.proc.ids.diff.cls.mat.2[i,1:tam.ids.proc[i]]=db.proc.ids.diff.cls.2[[i]]
		db.proc.ids.diff.cls.mat.2=cbind(db.proc.ids.diff.cls.mat.2,tam.ids.proc)
		
		gns.dnd.proc.2=lapply(gns.dnd.proc,function(x){
	           						aux.1=which(aux[x]!=0)
	           						x=aux[x][aux.1]
	           						x})
		tam.ids.gns=lapply(gns.dnd.proc.2,length)
		tam.ids.gns=unlist(tam.ids.gns)
		names(tam.ids.gns)=NULL
		nrow.gns.ids=length(gns.dnd.proc)
		ncol.gns.ids=max(tam.ids.gns)
		gns.dnd.proc.mat.2=matrix(0,nrow=nrow.gns.ids,ncol=ncol.gns.ids)
		for(i in 1:nrow.gns.ids)
   			gns.dnd.proc.mat.2[i,1:tam.ids.gns[i]]=gns.dnd.proc.2[[i]]
		gns.dnd.proc.mat.2=cbind(gns.dnd.proc.mat.2,tam.ids.gns)
	}
	FC.ntwrk.obj.0=1:length(db.proc.ids.diff.cls.2)
	FC.ntwrk.src.0=rep(0,length(db.proc.ids.diff.cls.2))
	nu.ceros=1
	while(nu.ceros!=0)
  	{
	    FC.ntwrk.src.0=sample(FC.ntwrk.obj.0,length(FC.ntwrk.obj.0))
   		nu.ceros=which((FC.ntwrk.obj.0-FC.ntwrk.src.0)==0)
   		nu.ceros=length(nu.ceros)
	}

	G.ntwrk.obj.0=1:length(gns.dnd.proc.2)
	G.ntwrk.src.0=rep(0,length(gns.dnd.proc.2))
	nu.ceros=1
	while(nu.ceros!=0)
  	{
   		G.ntwrk.src.0=sample(G.ntwrk.obj.0,length(G.ntwrk.obj.0))
   		nu.ceros=which((G.ntwrk.obj.0-G.ntwrk.src.0)==0)
	    nu.ceros=length(nu.ceros)
  	}
  	
  	result=list(datos,datos.fc,donde.affy,gns.dnd.proc.mat.2,db.proc.ids.diff.cls.mat.2,FC.ntwrk.obj.0,FC.ntwrk.src.0,G.ntwrk.obj.0,G.ntwrk.src.0)
  	nombres=c("gene.data","set.data","affy.loc","gene2set.mat","set2gene.mat","Set.obj","Set.src","G.obj","G.src")
  	names(result)=nombres
  	return(result)
}
