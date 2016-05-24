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



simc<-function(pdb,mj1=NULL,mj2=NULL,mj.avg=FALSE,cl=NULL,alpha=82,cuts=169,path=getwd(),inv2file=FALSE,bfacs=TRUE,frob=TRUE, loc=NULL,norm=FALSE,file=NULL,cluster=NULL){
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
	p<-extractPDB(pdb)

	seq<-p$caseq
	n<-p$lca
	b<-p$b
	d<-p$chains

	s<-aa2num(seq,0)

	interaction.mat<-build.interact(s,mj1,mj2,d,alpha)
	out<-build.contacts(n,cuts,p$coords)

	contact.mat<-out$cm
	deltas<-out$deltas
	cc<-out$cnr
	#umformatieren und auf gueltige contacts ueberpruefen
	if(is.null(cl)){
		cl<-get.contact.list(contact.mat,cumsum(d),single=TRUE)
	}else{
		cons<-get.contact.list(contact.mat,cumsum(d),single=FALSE)
		if(is.data.frame(cl)){
			#bedeutet: ist data.frame -> ok: check, format, unwahrscheinlicher fall
			cl<-as.list(unname(cl))
		}else if(!is.null(dim(cl))){
			#bedeutet: ist matrix -> ok: check, format
			if(dim(cl)[1]==2){
				cl<-cl
			}else if(dim(cl)[2]==2){
				cl<-t(cl)
			}else{
				stop("ERROR: not a contact list\n")
			}
			cl<-as.list(unname(as.data.frame(cl)))
		}
		cl<-cl[(cl%in%cons)]
	}
  lcl<-length(cl)
	cat("Number of Contacts:",lcl,"\n")
  if(cores > lcl){
    cores<-lcl
  }
  fn<-unlist(strsplit(gsub("\\\\","/",pdb),split="/"))
	fn<-fn[length(fn)]
	file.name<-unlist(strsplit(fn,split=".p"))[1]

	hess<-build.hess(contact.mat,interaction.mat,deltas)
	out<-get.svd(hess)

	evf<-out$ev
	ihf<-build.invhess(out)

	if(inv2file){
		write.table(ihf,file=paste(path,"/cov_mat_",file.name,"_0_0.out",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
		}
	if(bfacs){
		bees<-get.bfacs(ihf)
		write(bees,file=paste(path,"/b_factors_",file.name,"_0_0.out",sep=""),sep="\n")
		}
	if(norm){
		nihf<-mat.norm(ihf)
	}
	simc.func<-function(arg){
		library.dynam("BioPhysConnectoR",package="BioPhysConnectoR");
		contact.mat.new<-contact.mat
		contact.mat.new[arg[1],arg[2]]<-contact.mat.new[arg[2],arg[1]]<-0
		hessian.mat<-build.hess(contact.mat.new,interaction.mat,deltas)
		out<-get.svd(hessian.mat)
		cov.mat<-build.invhess(out)
		if(bfacs){
			bees<-get.bfacs(cov.mat)
		write(bees,file=paste(path,"/b_factors_",file.name,"_",arg[1],"_",arg[2],".out",sep=""),sep="\n")
			}
		if(inv2file){
		write.table(cov.mat,file=paste(path,"/cov_matrix_",file.name,"_",arg[1],"_",arg[2],".out",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
			}

		if(frob){
			if(norm){
				mat.cov<-mat.norm(cov.mat)
				mat.ih<-nihf
			}else{
				mat.cov<-cov.mat
				mat.ih<-ihf
			}
			loc<-rbind(c(1,n*3,1,n*3),loc)
			fn<-fnorm(mat.ih,mat.cov)
			l<-length(loc)
			if(l>0 && (l%%4==0)){
				fn<-vector()
				for(i in 1:(l/4)){
					ff<-fnorm(mat.ih[loc[i,1]:loc[i,2],loc[i,3]:loc[i,4]],mat.cov[loc[i,1]:loc[i,2],loc[i,3]:loc[i,4]])
					fn<-c(fn,ff)
				}
			}else{
				stop("wrong length")
			}
		}

		if(frob){
			return(c(arg[1],arg[2],fn))
			}
		}

	##callingfunction for parLapply / lapply
	fun<-function(cl){
		return(lapply(cl,simc.func))
	}
	##list formatting to get a list of lists with length cores
	overlap<-lcl%%cores
	ll<-floor(lcl/cores)
	ncl<-list()
	it<-0
	for(i in 1:cores){
		it<-it+1
		if(overlap>=i){
			itn<-it+ll
		}else{
			itn<-it+ll-1
		}
		ncl[[i]]<-cl[it:itn]
		it<-itn
	}

	##calling lapply
	if(!is.null(cluster)){
		tqq<-parLapply(cluster,ncl,fun)
	}else{
		tqq<-lapply(ncl,fun)
	}

    if(frob){
	    ##formatting to matrix
	    colmat<-length(tqq[[1]][[1]])
	    qq<-matrix(unlist(tqq),ncol=colmat,byrow=TRUE)
    }

	if(is.null(file)){
		file<-file.name
	}

	if(frob){
		write.table(qq,quote=FALSE,row.names=FALSE,col.names=FALSE,file=paste(path,"/results_",file,".out",sep=""))
        
    	return(qq)
    }
}

