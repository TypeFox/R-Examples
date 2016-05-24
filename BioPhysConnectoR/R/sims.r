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

sims<-function(pdb,alignment,mj1=NULL,mj2=NULL,mj.avg=FALSE,alpha=82,cuts=169,path=getwd(),mimethod="ORMI",gapchar="NOGAPCHAR",inv2file=FALSE,bfacs=TRUE,frob=TRUE,loc=NULL,norm=FALSE,cluster=NULL){
	# ### Begin of the original bio3d functions "read.fasta" and "aa123" as provided in bio3d 1.0-6 under GPL version2 by Grant, Rodrigues, ElSawy, McCammon, Caves, (2006) {Bioinformatics} 22, 2695--2696.
	aa123<-function (aa) {

	# convert one-letter IUPAC amino-acid code into
	# three-letter PDB style, for instance "A" into "ALA".

	aa1 <- c("-","X",
			"A","C","D","E","F","G",
			"H","I","K","L","M","N","P","Q",
			"R","S","T","V","W","Y")
	aa3 <- c("---","UNK",
			"ALA", "CYS", "ASP", "GLU", "PHE", "GLY",
			"HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
			"ARG", "SER", "THR", "VAL", "TRP", "TYR")

	convert <- function(x) {
		if(is.na(x)) return(NA)
		if (all(x != aa1)) {
		warning("Unknown one letter code for aminoacid")
	#      return(NA)
		return("UNK")
		}
		else {
		return(aa3[which(x == aa1)])
		}
	}
	return(as.vector(unlist(sapply(aa, convert))))
	}
	read.fasta<-function(file, rm.dup=TRUE, to.upper=FALSE, to.dash=TRUE) {

	## Version   0.3 ... Thu Apr 26 19:17:09 PDT 2007
	##                    uses scan instead of read.table

	raw.fa <- scan(file, what=character(0), sep="\n", quiet = TRUE)
	ind <- grep(">", raw.fa) ## seq id lines
	if(length(ind) == 0) {
		stop("read.fasta: no '>' id lines found, check file format")
	}

	if (to.dash) { raw.fa[-ind] <- gsub("[/.]","-", raw.fa[-ind]) }
	if (to.upper) { raw.fa[-ind] <- toupper(raw.fa[-ind]) }

	ind.s <- ind+1           ## seq start and end lines
	ind.e <- c((ind-1)[-1], length(raw.fa))

	seq.dim <-  apply(cbind(ind.s, ind.e), 1,
						function(x) sum( nchar(raw.fa[ (x[1]:x[2])]) ))

	seq.format <- function(x, max.seq=max(seq.dim)) {
		fa <- rep("-",max.seq)
		fa[ c(1:x[3]) ] <- unlist(strsplit( raw.fa[ (x[1]:x[2]) ], split=""));
		return(fa)
	}

	##seq.format( cbind(ind.s[1], ind.e[1], seq.dim[1]) )
	store.fa <- t(apply(cbind(ind.s, ind.e, seq.dim), 1, seq.format))
	rownames(store.fa) <- gsub("^>| .*", "",raw.fa[ind], perl=TRUE)

	##  if (to.dash) { store.fa <- gsub("[/.]","-", store.fa ) }
	##  if (to.upper) { store.fa <- toupper(store.fa) }

	if (rm.dup) {        ## remove duplicated seq id's
		dups <- as.numeric(duplicated(row.names(store.fa)))
		if (sum(dups) > 0) {
		print(paste(" ** Duplicated sequence id's: ",
					row.names(store.fa)[dups]," **",sep=""))
		store.fa <- store.fa[!dups,]
		}
	}

	output <- list(id=rownames(store.fa), ali=store.fa)
	class(output) <- "fasta"
	return(output)
	}
	# ### End of bio3d functions
	cores<-1
	if(!is.null(cluster)){
		clusterEvalQ(cluster,library(BioPhysConnectoR))
		cores<-length(cluster)
	}

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

	p<-extractPDB(pdb)
	fn<-unlist(strsplit(gsub("\\\\","/",pdb),split="/"))
	fn<-fn[length(fn)]
	file.name<-unlist(strsplit(fn,split=".pdb"))[1]
	aln<-read.fasta(alignment)$ali

	if(mimethod=="ORMI"){
		bool<-FALSE
	}else{
		bool<-TRUE
	}

	H<-get.entropy(aln,bool=bool,gapchar=gapchar)
	MI<-get.mie(aln,method=mimethod,gapchar=gapchar)

	seq<-p$caseq
	n<-p$lca
	b<-p$b
	d<-p$chains

	s<-aa2num(seq,0)

	out<-build.contacts(n,cuts,p$coords)

	contact.mat<-out$cm
	deltas<-out$deltas
	cc<-out$cnr

	interaction.mat<-build.interact(s,mj1,mj2,d,alpha)

	hess<-build.hess(contact.mat,interaction.mat,deltas)
	out<-get.svd(hess)
	evf<-out$ev
	ih<-build.invhess(out)
	if(norm){
		nih<-mat.norm(ih)
	}
	if(inv2file){
		write.table(ih,file=paste(path,"/cov_matrix_",file.name,"_pdb.out",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
		}
	if(bfacs){
		beta<-get.bfacs(ih)
		write(beta,file=paste(path,"/b_factors_",file.name,"_pdb.out",sep=""),sep="\n")
		}

	sims.func<-function(arg){
		library.dynam("BioPhysConnectoR",package="BioPhysConnectoR");
		seq<-aa123(aln[arg,])
		s<-aa2num(seq,0,verbose=FALSE)

		interaction.mat<-build.interact(s,mj1,mj2,d,alpha)

		hessian.mat<-build.hess(contact.mat,interaction.mat,deltas)
		out<-get.svd(hessian.mat)
		covariance.mat<-build.invhess(out)

		if(bfacs){
			bees<-get.bfacs(covariance.mat)
			write(bees,file=paste(path,"/b_factors_",file.name,"_",arg,"_",rn[arg],".out",sep=""),sep="\n")
			}
		if(frob){
			if(norm){
				mat.cov<-mat.norm(covariance.mat)
				mat.ih<-nih
			}else{
				mat.cov<-covariance.mat
				mat.ih<-ih
			}
			loc<-rbind(c(1,n*3,1,n*3),loc)
			fno<-fnorm(mat.ih,mat.cov)
			l<-length(loc)
			if(l>0 && (l%%4==0)){
				fno<-vector()
				for(i in 1:(l/4)){
					ff<-fnorm(ih[loc[i,1]:loc[i,2],loc[i,3]:loc[i,4]],mat.cov[loc[i,1]:loc[i,2],loc[i,3]:loc[i,4]])
					fno<-c(fno,ff)
				}
			}else{
				stop("wrong length")
			}
			}
		if(inv2file){
			write.table(covariance.mat,file=paste(path,"/covariance_matrix_",file.name,"_",arg,"_",rn[arg],".out",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
			}
		if(frob){
			return(c(arg,rn[arg],fno))
			}
	}

	fun<-function(seq.list){
		return(lapply(seq.list,sims.func))
	}

	rn<-rownames(aln)
	seq.list<-as.list(1:nrow(aln))

	lsl<-length(seq.list)
  if(cores > lsl){
    cores<-lsl
  }
	overlap<-lsl%%cores
	ll<-floor(lsl/cores)
	nsl<-list()
	it<-0
	for(i in 1:cores){
		it<-it+1
		if(overlap>=i){
			itn<-it+ll
		}else{
			itn<-it+ll-1
		}
		nsl[[i]]<-seq.list[it:itn]
		it<-itn
	}
	seq.list<-nsl
	if(!is.null(cluster)){
        	qq<-parLapply(cluster,seq.list,fun)
	}else{
		qq<-lapply(seq.list,fun)
	}

	if(frob){
		colmat<-length(qq[[1]][[1]])
		mqq<-matrix(unlist(qq),ncol=colmat,byrow=TRUE)
			write.table(mqq,file=paste(path,"/results_",file.name,".out",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
	ret<-list(entropy=H,mi=MI,res=mqq)
	return(ret)
	}
