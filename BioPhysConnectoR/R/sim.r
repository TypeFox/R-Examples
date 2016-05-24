#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



sim<-function(pdbs,mj1=NULL,mj2=NULL,mj.avg=FALSE,alpha=82,cuts=169,path=getwd(),cluster=NULL){
	if(is.null(mj1)){
		mj1<-mat.read(system.file("mj1.txt",package="BioPhysConnectoR"))
		}
	if(is.null(mj2)){
		mj2<-mat.read(system.file("mj2.txt",package="BioPhysConnectoR"))
		}
	if(mj.avg){
		mj1<-matrix(data=mean(mj1),nrow=20,ncol=20)
		mj2<-matrix(data=mean(mj2),nrow=20,ncol=20)
		}
	cores<-1
        if(!is.null(cluster)){
		clusterEvalQ(cluster,library(BioPhysConnectoR))
		cores<-length(cluster)              
	}

	sim.func<-function(arg){
		library.dynam("BioPhysConnectoR",package="BioPhysConnectoR");

		p<-extractPDB(arg)
		fn<-unlist(strsplit(gsub("\\\\","/",arg),split="/"))
		fn<-fn[length(fn)]
		file.name<-unlist(strsplit(fn,split=".pdb"))[1]
		seq<-p$caseq
		n<-p$lca
		nr3<-3*n
		b<-p$b
		d<-p$chains

		s<-aa2num(seq,0,verbose=FALSE) #0 da es an C-Programm weitergereicht wird

		interaction.mat<-build.interact(s,mj1,mj2,d,alpha)

		out<-build.contacts(n,cuts,p$coords)

		contact.mat<-out$cm
		deltas<-out$deltas

		cov.mat<-get.cov(contact.mat,interaction.mat,deltas);
		beta<-get.bfacs(cov.mat)
		write.table(cov.mat,col.names=FALSE,row.names=FALSE,quote=FALSE,file=paste(path,"/cov_matrix_",file.name,".out",sep=""))
		write(beta,file=paste(path,"/b_factors_",file.name,".out",sep=""),sep="\n")

	}

	fun<-function(pdbs){
		return(lapply(pdbs,sim.func))
	}

	pdbs<-as.list(pdbs)

	lpdb<-length(pdbs)
  if(cores > lpdb){
    cores<-lpdb
  }
	ll<-floor(lpdb/cores)
	overlap<-lpdb%%cores
	npdb<-list()
	it<-0
	for(i in 1:cores){
		it<-it+1
		if(overlap>=i){
		itn<-it+ll
		}else{
			itn<-it+ll-1
		}
		npdb[[i]]<-pdbs[it:itn]
		it<-itn
	}
	pdbs<-npdb
	if(!is.null(cluster)){
        	parLapply(cluster,pdbs,fun)
	}else{
		lapply(pdbs,fun)
	}
	return(1)
}
