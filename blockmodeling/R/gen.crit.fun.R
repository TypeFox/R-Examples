"gen.crit.fun" <-
function(
	#function for generting a function for computing criteria function of a blockmodel and preparing data
	M,	#matrix or an array in case of multirelational networks with dimensions N x N x R, where N is the number of vertices and R the number of relations,
	k,	#number partitions for each mode
#	e1="default",	#weight of the error of binearized matrix
#	e2="default",	#weight of the error of valued conenctions
	approach,	#the approach used - can be one of "ss","ad","bin","val", "imp" (yet to be implementer) , "bv", "bi". A vector if different approaches should be used for different relations.
	cut = if(length(dim(M))==2) min(M[M>0]) else apply(M,3,function(x)min(x[x>0])),	#the cutting parameter used to binerize a valued network	m=NULL,	#suficient value for individual cells, can be specified as a number or a function of a block (or only one of its rows or coulms anslyzed) (max for implicit appraoch)
	m= NULL,
#	s="default",	#suficient value for colum and row statistics
	FUN="max",	#function to calculate row and colum statistics
	blocks=NULL,	#permissible block types and their ordering, can be also on of 'structural', 'regular', 'regular.ext' or 'all'. A list if different permissible block types and their ordering should be applied to different relations.
	BLOCKS=NULL,	#prespecified model, a list if different blockmodels should be applied to different relations
	diag = TRUE,	#should the diagonal blocks be treated as diagonal blocks
	relWeights = NULL,	#the weights for different relations in the array M - can be also used for combining criterion functions of several approaches to blockmodeling
	normMto2rel=FALSE, #create two-realation netowrk from one relational network through row and column normalization.
	normMto2relRegToRreCre = normMto2rel,	#when normMto2rel is TRUE, should regular blocks be converted to row-regular in row-normalized relation and to column-regular in column-normalized relation.
	sameModel=normMto2rel, #should we damand the same blockmodel for all relations. If set to TRUE, it demands that accros all relations the ideal block on the same position in the matrix BLOCKS should be chosen. Usually, these positions are occupied by the same blocks. If not, use with caution.
	mindim = 2,	#minimal dimension for rew/column regular, dominant and functional blocks
	mindimreg = FALSE, #should mindim be also used for regular blocks
	regDiagSep = FALSE, #should the diagonal be treated sepeartly in regular blocks -> not used in evalution of funtion f by rows or columns but as the diagonal in null and complete blocks? Currently only supported for homogeneity blockmodeling
	blockWeights=NULL,	#weights for all block types - if only some of them are specified, the other remain 1
	positionWeights=NULL,	#weigths for positions in the blockmodel (the dimensions must be the same as the error matrix)
	save.err.v=sameModel,	#save a vector of errors of all block tipes for all blocks
	changeT=FALSE,	#do we also want a function for computing only the chages in the blockmodel/errors
	dn="default",	#the density reshold for density block - default = mean(M>=cut)
	av="default",	#the treshold for average block - default = mean(M)
	norm= FALSE,	#should the block errors be normalised by the size of the blocks
	normbym = FALSE,#should the block errors be normalised by m
	allow.max0 = !"null" %in% unique(c(unlist(blocks),unlist(BLOCKS))),	#Should the maximum that is the basis for calculation of inconsistencies in implicit blockmodeling be allowed to be 0. If FALSE, the maximum is in such case set to the maximum of the network (if maximum of a block is 0) or to the maxsimum of the block (if row or column maxismum is 0)
			#Used only in implicit blockmodeling
	allow.dom0=FALSE,	#should the dominant row or column be allowed to be 0. Used only in implicit blockmodeling.
	domMax=TRUE,	#should it be allowed that the dominant row or column is not the one with the largest sum, max,... Not yet used.
	max.con.val=if(normbym && approach!="imp") "m" else "non" ,	#used only in valued and implicit blockmodeling - should the largest values be cencored, limited to (larger values set to) - resonoble values are:
																# "m" or "s" - the maximum value equals the value of the parameter m/s
																# "non" - no transformation is done
																# numeriacal values larger then parameters m and s and lower the the maximum
	BLOCKS.CV=NULL, #prespecified block central values for homogeneity blockmodeling - must have the same dimensions as BLOCKS
	CV.use=NULL, #how should the prespecified block central values for homogeneity blockmodeling be used - possible values are "fixed","min" and "max".
	return.CV=TRUE,	#should a matrix of block central values (for homogeneity blockmodeling only) be returned.
	use.for=TRUE, #should fotran subrutines be used when possible
  ...
){
#	  cat("Start: approach = ",approach,", normbym = ", normbym,", norm = ", norm, "\n", sep="")
	
#	norm
#	normbym
#	max.con.val

	if(!is.null(positionWeights))assign(x="pWeights",value=positionWeights,envir=parent.frame())
	
	if(is.null(blocks)&&is.null(BLOCKS))stop("Either blocks or BLOCKS must be specified!\n")
	bin.all.ideal.blocks<-c("null","com","rdo","cdo","reg","rre","cre","rfn","cfn","den","dnc")
	val.all.ideal.blocks<-c("null","com","rdo","cdo","reg","rre","cre","rfn","cfn","avg","dnc")
	hom.all.ideal.blocks<-c("null","com","rdo3","cdo3","rdo1","cdo1","rdo2","cdo2","rfn","cfn","reg","rre","cre","dnc") #blocks "rdo1","cdo1","rdo2","cdo2","rfn","cfn" are in testing phase and are also not implemented for pre-specified blocks central values.
	tmpfun<-"crit.fun.tmp"
 	M<-as.array(M)
	if(length(dim(M))==3){
 	  nr<-dim(M)[3]
  } else nr<-1

	if(sameModel){
		if(!is.null(BLOCKS)){
			if(is.list(blocks)){
				for(i in 2:length(blocks)){
					if(!all(blocks[[1]]==blocks[[i]])){
						warning("sameModel used with diferent allowed blocks for different relations\n")
						break
					}
				}
			}
		}else{
			if(is.list(BLOCKS)){
				for(i in 2:length(BLOCKS)){
					if(!all(BLOCKS[[1]]==BLOCKS[[i]])){
						warning("sameModel used with diferent BLOCKS for different relations\n")
						break
					}
				}
			}
		}
	}
	
	nmode<-length(k)
  	
  
	if(nmode==1){
		k<-c(k,k)
	} else diag<-FALSE
	
	if(normMto2rel){
		if(is.matrix(M)){
			LM<-array(NA,dim=c(dim(M),2))
			LM[,,1]<-diag(1/rowSums(M))%*%M
			LM[,,2]<- M%*%diag(1/colSums(M))
			LM[is.nan(LM)]<-0
			M<-LM
			if(is.null(BLOCKS)){
				BLOCKS.org<-array(blocks,dim=c(length(blocks),k[1],k[2])) #,dimnames=list(NULL,nclu,nclu)
			} else if(!is.list(BLOCKS)) BLOCKS.org<-BLOCKS
			
			if(!is.list(BLOCKS)){
				BLOCKS2<-BLOCKS1<-BLOCKS.org
				if(normMto2relRegToRreCre){
					BLOCKS1[BLOCKS1=="reg"]<-"rre"
					BLOCKS2[BLOCKS2=="reg"]<-"cre"
				}
				BLOCKS<-list(BLOCKS1,BLOCKS2)
			} else if(length(BLOCKS)!=2) stop("BLOCKS must be a list of legnth 2 for the use with normMto2rel!\n")
			nr<-2
		}else stop("normMto2rel = TRUE can be only used on one relational networks.\n")
	}
	if(is.matrix(M)) class(M)<-"mat"

	if(nr==1){
		directed<-nmode==2||!all(M==t(M))
	} else{
		directed<-nmode==2||!all(apply(M,3,function(x)all(x==t(x))))
	}
	
	if(nr!=length(approach)){
		if(nr==1) M<-array(M,dim=c(dim(M),length(approach)))
	}


	if(nmode>2){
		kmode<-k
		k<-c(sum(k),sum(k))
	}


	if(nr==1){
		sameModel<-FALSE
		if(is.null(BLOCKS)){
			BLOCKS<-array(blocks,dim=c(length(blocks),k[1],k[2])) #,dimnames=list(NULL,nclu,nclu)
		}else {
			blocks<-unique(na.omit(as.vector(BLOCKS)))
		}
		n.types<-length(blocks)
	}else{
		if(is.null(relWeights))relWeights<-rep(1,nr)
		if(length(relWeights)!=nr) stop("Length of relWeights does not match the number of relations (the third dimmension of M).")
		if(is.null(BLOCKS)){
				if(!is.list(blocks)){
					blocks<-rep(list(blocks),times=nr)
				}
				BLOCKS<-lapply(blocks,function(x)array(x,dim=c(length(x),k[1],k[2])))
			}else {
				if(!is.list(BLOCKS)){
				BLOCKS<-rep(list(BLOCKS),times=nr)
				}
				blocks<-lapply(BLOCKS,function(x)unique(as.vector(x)))
		}
		n.types<-sapply(blocks,length)
	}

	#	dimnames(BLOCKS)<-list(NULL,nclu,nclu)
	#}
	
	 assign(x="expBLOCKS",value=BLOCKS,envir=parent.frame())

	
	if(sameModel&&nr>1&&(sum(apply(sapply(BLOCKS,dim),1,ss))!=0)){
		stop("sameModel used with BLOCKS with diferent dimensions\n")
	}


	normal<-TRUE

	
	if(nr==1 && all(approach=="ss") && all(blocks=="com")&& all(BLOCKS=="com") && nmode==1 && use.for){
		assign(x="useM",value=M,envir=parent.frame())
		fun1<-(parse(text=paste(c(tmpfun,"<-function(M,clu,...){
		n<-dim(M)[1]
		dn<-dimnames(M)
		k<-as.integer(",k[1],")
		clu<-as.integer(factor(clu))
		res<- .Fortran('critfunsscom',M=matrix(as.double(M),ncol=n),n=n,clu=clu, k=k,diag=",diag,",err=0,E=diag(k), BM=diag(k))
		class(res)<-'crit.fun'
		dn->dimnames(res$M)
    IM<-res$E
    IM[,] <- 'com'
    res$IM<-IM
		return(res[c('err','clu','E','IM','BM')])}
	"),sep="",collapse="")))

		if(changeT){
			return(list(fun1 = fun1, fun2 = fun1))
		} else {
			return(fun1)
		}

		normal<-FALSE
	}

	if(normal){
		if(nr==1){
			if(!is.null(blockWeights)){
				use.weights<-TRUE
				w<-rep(1, times=length(blocks))
				names(w)<-blocks
				w[names(blockWeights)]<-blockWeights
				blockWeights<-w
				assign(x="bweights",value=blockWeights,envir=parent.frame())
			} else use.weights<-FALSE
		}else{
			if(!is.null(blockWeights)){
				use.weights<-TRUE
				ww<-NULL
				for(i in 1:nr){
					w<-rep(1, times=length(blocks[[i]]))
					names(w)<-blocks[[i]]
					if(is.list(blockWeights)){
						w[names(blockWeights[[i]])]<-blockWeights[[i]]
					}else{
						w[names(blockWeights)]<-blockWeights
					}
					ww<-c(ww,list(w))
				}
				blockWeights<-ww
				assign(x="bweights",value=blockWeights,envir=parent.frame())
			} else use.weights<-FALSE
		}




		fun<-c(tmpfun,"<-function(M,clu,mindim=",mindim,if(use.weights)",blockWeights=bweights" else NULL, if(sameModel)",BLOCKS=expBLOCKS" else NULL,if(!is.null(positionWeights))",positionWeights=pWeights",if(any(approach%in%c("val","imp")))",m=exm","){\n")
		if(changeT) {fun2<-c(tmpfun,"<-function(M,clu,mindim=",mindim,if(use.weights)",blockWeights=bweights" else NULL, if(sameModel)",BLOCKS=expBLOCKS" else NULL,if(!is.null(positionWeights))",positionWeights=pWeights",if(any(approach%in%c("val","imp")))",m=exm",",change",if(nmode==2)", modechange" else NULL,",res.old){\n")}
		fun<-c(fun,"E<-array(NA,dim=c(",k[1],",",k[2],if(nr>1)c(",",nr)else NULL,"))\n","IM<-array(NA,dim=c(",k[1],",",k[2],if(nr>1)paste(",",nr,sep="")else NULL,"))\n")
		if(changeT) {
			fun2<-c(fun2,"E<-res.old$E\n",if(normMto2rel) "IM = res.old$completeIM\n" else "IM<-res.old$IM\n")	#image matrix - matrix of types of blocks #,dimnames=list(nclu,nclu)
			length.fun<-length(fun)
		}


		if(save.err.v||sameModel) {
			if(nr==1){
				fun<-c(fun,"ERR.V<-array(NA,dim=c(",n.types,",",k[1],",",k[2],"),dimnames = list(",paste("c('",paste(blocks,collapse="','"),"')",sep=""),",NULL,NULL))\n") #,dimnames=list(nclu,nclu,all.block.types)
			}else{
				if(sameModel){
					fun<-c(fun,"ERR.V<-array(NA,dim=c(",n.types[1],",",k[1],",",k[2],",",nr,"))\n") #,dimnames=list(nclu,nclu,all.block.types)
				}else{
					fun<-c(fun,"ERR.V<-list(array(NA,dim=c(",n.types[1],",",k[1],",",k[2],"),dimnames = list(",paste("c('",paste(blocks[[1]],collapse="','"),"')",sep=""),",NULL,NULL))\n") #,dimnames=list(nclu,nclu,all.block.types)
					for(i in 2:nr){
						fun<-c(fun,"ERR.V<-c(ERR.V,list(array(NA,dim=c(",n.types[i],",",k[1],",",k[2],"),dimnames = list(",paste("c('",paste(blocks[[i]],collapse="','"),"')",sep=""),",NULL,NULL))\n") #,dimnames=list(nclu,nclu,all.block.types)
					}
				}
			}
			if(changeT&&save.err.v) {
				fun2<-c(fun2,"ERR.V<-res.old$ERR.V\n")
				length.fun<-length(fun)
			}
		}

		if(nr>1){
			if(length(FUN)==1) FUN<-rep(list(FUN),times=nr)			
			if(length(approach)==1) approach<-rep(approach,times=nr)
			if(length(m)==1) m<-rep(m,times=nr)
			if(length(normbym)==1) normbym<-rep(normbym,times=nr)
			if(length(norm)==1) norm<-rep(norm,times=nr)
			if(length(cut)==1) cut<-rep(cut,times=nr)
			if(length(max.con.val)==1) max.con.val<-rep(max.con.val,times=nr)
		}
		

		if(any(approach=="imp")){
			if(is.null(m)) m<-rep(1,times=length(approach))
			m[approach=="imp"]<-"max"
			approach[approach=="imp"]<-"val"
			print(m)
		}

		if(any(approach %in% c("ad","ss"))){
			if(any(approach == "ss")){
				fun<-c(fun,"dev<-function(x)sum((x-mean(x))^2)\n")
				fun<-c(fun,"idev<-function(x)(x-mean(x))^2\n")
				if(is.null(CV.use)){
					fun<-c(fun,"dev.cv<-function(x,cv)sum((x-cv)^2)\n")
					fun<-c(fun,"idev.cv<-function(x,cv)(x-cv)^2\n")			
				}else{
					fun<-c(fun,"dev.cv<-function(x,cv,use='fixed'){cv<-switch(use,min=max(cv,mean(x)),max=min(cv,mean(x)),fixed=cv,free=mean(x));sum((x-cv)^2)}\n")
					fun<-c(fun,"idev.cv<-function(x,cv,use='fixed'){cv<-switch(use,min=max(cv,mean(x)),max=min(cv,mean(x)),fixed=cv,free=mean(x));(x-cv)^2}\n")				
				}
				fun<-c(fun,"cv<-mean\n")
			} else{
				fun<-c(fun,"dev<-function(x)sum(abs(x-median(x)))\n")
				fun<-c(fun,"idev<-function(x)abs(x-median(x))\n")
				if(is.null(CV.use)){
					fun<-c(fun,"idev.cv<-function(x,cv)abs(x-cv)\n")				
					fun<-c(fun,"dev.cv<-function(x,cv)sum(abs(x-cv))\n")
				}else{
					fun<-c(fun,"dev.cv<-function(x,cv,use='fixed'){cv<-switch(use,min=max(cv,median(x)),max=min(cv,median(x)),fixed=cv);sum((x-cv)^2)}\n")
					fun<-c(fun,"idev.cv<-function(x,cv,use='fixed'){cv<-switch(use,min=max(cv,median(x)),max=min(cv,median(x)),fixed=cv);(x-cv)^2}\n")				
				}
				fun<-c(fun,"cv<-median\n")
			}
			approach[approach %in% c("ad","ss")]<-"hom"
			if(any(normbym[approach %in% c("ad","ss")])) {
				normbym[approach %in% c("ad","ss")]<-FALSE
				warning("'normbym' can be only used for valued ('val') and implicit ('imp') apporach")
			}
		}

		#### for binary approach
		if(any(approach=="bin")){
			if(nr==1){
				if(!all(blocks%in%bin.all.ideal.blocks)) stop("The block(s)",blocks[!blocks%in%bin.all.ideal.blocks] ," is (are) not defined!")
				bM <- (M>=cut)+0
			}else {
				for(iblocks in blocks[approach=="bin"]){
					if(!all(iblocks%in%bin.all.ideal.blocks)) stop("The block(s)",blocks[!blocks%in%bin.all.ideal.blocks] ," is (are) not defined!")
				}
				bM <- M
				for(i in 1:nr) bM[,,i]<-(M[,,i]>=cut[i])+0
			}

      assign(x="bM",value=bM,envir=parent.frame())
			if("den"%in% unlist(blocks)) fun<-c(fun,"dn<-",dn,"\n")
			if(any(normbym[approach=="bin"])) {
				normbym[approach=="bin"]<-FALSE
				warning("'normbym' can be only used for valued ('val') and implicit ('imp') apporach")
			}
		}


		#### for valued approach
		if(any(approach=="val")){
			if(is.null(m))stop("m must be specified for valued blockmodeling")
			if(nr==1){
				if(!all(blocks%in%val.all.ideal.blocks)) stop("The block(s)",blocks[!blocks%in%val.all.ideal.blocks] ," is (are) not defined!")
			}else {
				for(iblocks in blocks[approach=="val"])
				if(!all(iblocks%in%val.all.ideal.blocks)) stop("The block(s)",blocks[!blocks%in%val.all.ideal.blocks] ," is (are) not defined!")
			}


			if(any(max.con.val!="non")){
				if(nr==1){
					if(max.con.val=="m")max.con.val<-m
					if(max.con.val=="s")max.con.val<-s
					if(!is.numeric(max.con.val))stop('"max.con.val" must me numeric or a strig with a value "m", "s" or "non"!')
	
					M[M>max.con.val]<-max.con.val
					M[M<(-max.con.val)]<- (-max.con.val)
				}else{
					for(i in 1:nr){
						tmax.con.val<-max.con.val[i]
						if(tmax.con.val=="m")tmax.con.val<-m[i]
						if(tmax.con.val=="s")tmax.con.val<-s[i]						
						tmax.con.val<-as.numeric(tmax.con.val)
						if(is.na(tmax.con.val))stop('"max.con.val" must me numeric or a strig with a value "m", "s" or "non"!')
						Mi<-M[,,i]
						Mi[Mi>tmax.con.val]<-tmax.con.val
						Mi[Mi<(-tmax.con.val)]<- (-tmax.con.val)
						Mi->M[,,i]
					}
				}
			}
			
			misfun<-!grepl(pattern="^[1234567890]+$",x=m)
			mfun<-m

	#		fun<-c(fun,"m <- ",m,"\n")
			if("avg"%in% unlist(blocks)) fun<-c(fun,"av<-",av,"\n")
		}else{
			misfun<-NULL
			mfun<-NULL
		}

		assign(x="useM",value=M,envir=parent.frame())



		if(nr>1){
#			Lm<-m
			LFUN<-FUN
			Lcut<-cut
			LBLOCKS<-BLOCKS
			Lblocks<-blocks
			Lmfun<-mfun
			Lapproach<-approach
			Lnormbym<-normbym
			Lnorm<-norm
			LBLOCKS.CV<-BLOCKS.CV
			LCV.use<-CV.use
			Lmisfun<-misfun
		}


		if(!is.null(m)){
			assign(x="exm",value=m,envir=parent.frame())
			if(nr>1) fun<-c(fun,"Lm<-m\n")
		}
		
		
		for(i1 in 1:k[1]){
			for(i2 in if(directed){1:k[2]}else{1:i1}){
				for(inr in seq(length=nr)){
					if(nmode<=2||{cl<-cut(c(i1,i2),breaks=c(0,cumsum(kmode))+0.5,labels=FALSE);cl[1]!=cl[2]}){
						if(nr>1){
#							Lm[inr]->m
							if(!is.null(m)) fun<-c(fun,"m<-Lm[",inr,"]\n")
							LFUN[[inr]]->FUN
							Lcut[inr]->cut
							LBLOCKS[[inr]]->BLOCKS
							LBLOCKS.CV[[inr]]->BLOCKS.CV
							LCV.use[[inr]]->CV.use
							Lblocks[[inr]]->blocks
							Lmisfun[inr]->misfun
							Lmfun[inr]->mfun
							Lapproach[inr]->approach
							Lnormbym[inr]->normbym
							Lnorm[inr]->norm
						}
						if(changeT) {
							fun2<-c(fun2,fun[-(1:length.fun)],"if(any(c(",if(nmode==2)c("if(modechange==1) ",i1," else ",i2) else c(i1,",",i2),") %in% change)){\n")
							length.fun<-length(fun)
						}


						if(approach=="bin"){
							fun<-c(fun,"B<-bM[clu",if(nmode==2) "[[1]]" else NULL,"==",i1,",clu",if(nmode==2) "[[2]]" else NULL,"==",i2,if(nr>1) c(",",inr) else NULL,",drop=FALSE]\n")
							fun<-c(fun,"nr<-dim(B)[1]\n")	#numer of rows
							fun<-c(fun,"nc<-dim(B)[2]\n")	#numer of colums
							#if(any(c("null","com","rdo","cdo","reg","rre","cre","rfn","cfn","den","dnc") %in% BLOCKS[,i1,i2]))
							#if(any(c("reg","rre","cre") %in% BLOCKS[,i1,i2])&i1==i2&diag&regDiagSep)fun<-c(fun,"Bd<-B\ndiag(Bd)<-0\n")							
							if(any(c("null","com","rfn","cfn","den") %in% BLOCKS[,i1,i2])) fun<-c(fun,"st<-sum(B)\n")
							if(any(c("reg","rdo","rre","rfn") %in% BLOCKS[,i1,i2])) fun<-c(fun,"sr<-rowSums(B)\n")
							if(any(c("reg","cdo","cre","rfn") %in% BLOCKS[,i1,i2])) fun<-c(fun,"sc<-colSums(B)\n")
							if(any(c("reg","rre","rfn") %in% BLOCKS[,i1,i2])) fun<-c(fun,"pr<-sum(sr>0)\n")
							if(any(c("reg","cre","rfn") %in% BLOCKS[,i1,i2])) fun<-c(fun,"pc<-sum(sc>0)\n")
						} else 	{
							fun<-c(fun,"B<-M[clu",if(nmode==2) "[[1]]" else NULL,"==",i1,",clu",if(nmode==2) "[[2]]" else NULL,"==",i2,if(nr>1) c(",",inr) else NULL,",drop=FALSE]\n")
							#if(any(c("reg","rre","cre") %in% BLOCKS[,i1,i2])&i1==i2&diag&regDiagSep)fun<-c(fun,"Bd<-B\ndiag(Bd)<-0\n")
						}
						fun<-c(fun,"dim(B)<-dim(B)[1:2]\n")
						
						fun<-c(fun,"nr<-dim(B)[1]\n")	#numer of rows
						fun<-c(fun,"nc<-dim(B)[2]\n")	#numer of colums
						
						if(approach=="val"){
							if(normbym) fun<-c(fun,"if(max(B)!=0){\n")
							#if(any(c("null","com","rdo","cdo","reg","rre","cre","rfn","cfn","den","dnc") %in% BLOCKS[,i1,i2]))
							if(misfun){
								fun<-c(fun,"m<-do.call(",mfun,",list(B))\n")
								if(!allow.max0) fun<-c(fun,"if(m==0)m<-do.call(",mfun,",list(M",if(nr>1) c("[,,",inr,"]") else NULL,"))\n")
							}

							#if(any(c("com","reg","rre","cre","rdo","cdo","rfn","cfn","avg") %in% BLOCKS[,i1,i2]) || norm || (i1==i2 && diag)) fun<-c(fun,"nr<-dim(B)[1]\n")	#numer of rows
							#if(any(c("com","reg","rre","cre","rdo","cdo","rfn","cfn","avg") %in% BLOCKS[,i1,i2]) || norm) fun<-c(fun,"nc<-dim(B)[2]\n")	#numer of colums
							if(any(c("null","rfn","cfn","avg") %in% BLOCKS[,i1,i2])) fun<-c(fun,"sumB<-sum(B)\n")


							if(any(c("reg","rre") %in% BLOCKS[,i1,i2])) fun<-c(fun,"sr<-apply(B,1,",FUN,")\n",
													"er<-m-sr",if(!misfun) c("[sr<m]") else NULL,"\n", #row errors for regular blocks
													"pr <- sum(sr>0)\n")

							if(any(c("reg","cre") %in% BLOCKS[,i1,i2])) fun<-c(fun,"sc<-apply(B,2,",FUN,")\n",
													"ec<-m-sc",if(!misfun) c("[sc<m]") else NULL,"\n", #column errors for regular blocks
													"pc <- sum(sc>0)\n")

							if("rfn" %in% BLOCKS[,i1,i2]) {
								if(any(c("reg","rre") %in% BLOCKS[,i1,i2]) && FUN=="max"){
									fun<-c(fun,"srm<-sr\n")
								}else fun<-c(fun,"srm<-apply(B,1,max)\n")
							}
							if((c("cfn") %in% BLOCKS[,i1,i2])) {
								if(any(c("reg","cre") %in% BLOCKS[,i1,i2]) && FUN=="max"){
									fun<-c(fun,"scm<-sc\n")
								}else fun<-c(fun,"scm<-apply(B,2,max)\n")
							}
						}

						if(approach == "hom"){
							if(any(c("null","com") %in% BLOCKS[,i1,i2]) && diag && i1==i2) fun<-c(fun,"B2<-B\ndiag(B2)<-NA\n")
							if(any(c("reg","rre","cre") %in% BLOCKS[,i1,i2]) && i1==i2 && diag && regDiagSep)fun<-c(fun,"Bd<-B\ndiag(Bd)<-0\n")
							#fun<-c(fun,"nr<-dim(B)[1]\n")	#numer of rows
							#fun<-c(fun,"nc<-dim(B)[2]\n")	#numer of colums
						}

						fun<-c(fun,"err<-numeric()\n")


						
						for(iEq in 1:dim(BLOCKS)[1]){
							eq=BLOCKS[iEq,i1,i2]
							if(!is.na(eq)){
								if(approach=="bin") fun<-c(fun, switch(EXPR=eq,
									"null" = c("err['",eq,"'] <-", "st","\n"),
									"com" = c("err['",eq,"'] <-", "nc*nr-st","\n"),
									"rdo" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , "(nc - max(sr))*nr","\n"),
									"cdo" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , "(nr - max(sc))*nc","\n"),
									"reg" = c("err['",eq,"'] <-",if(mindimreg&mindim>1)"if(nc<mindim|nr<mindim) Inf else " else NULL , "((nc-pc)*nr+(nr-pr)*pc)","\n"),
									"rre" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , "(nr-pr)*nc","\n"),
									"cre" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , "(nc-pc)*nr","\n"),
									"rfn" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , "st - pr + (nr - pr)*nc","\n"),
									"cfn" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , "st - pc + (nc - pc)*nr","\n"),
									"den" = c("err['",eq,"'] <-", "max(0,dn*nr*nc-st)","\n"),
									"dnc" = c("err['",eq,"'] <-", "0","\n")
								))

								if(approach=="val") fun<-c(fun, switch(EXPR=eq,
									"null" = c("err['",eq,"'] <-", "sumB",if(normbym) "/m","\n"),
									"com" = c("err['",eq,"'] <-", "sum(m-B[B<m])",if(normbym) "/m","\n"),
									"rdo" = c(if(mindim>1) c("if(nr<mindim) err['",eq,"'] <-Inf else"),"{",if(misfun) c("mr<-apply(B,1,",mfun,")\n", if(!allow.dom0) c("mr[mr==0]<-m\n")),"err['",eq,"'] <- min(",if(normbym)"(" else NULL,"apply(",if(misfun) c("rep(mr,times=nc) - B,1,sum)*nr") else c("m - B,1,sumpos)*nr"),if(diag && i1==i2){c(" + useneg(sum(diag(B))*nr - ",if(misfun) "(diag(mr - B))*nr)" else c("usepos(diag(m - B))*nr)"))} else NULL,if(normbym) {if(misfun) {if(allow.dom0)")/{mr[mr==0]<-1;mr})" else ")/mr)"} else c(")/m)")} else ")","}\n"),
									"cdo" = c(if(mindim>1) c("if(nc<mindim) err['",eq,"'] <-Inf else"),"{",if(misfun) c("mc<-apply(B,2,",mfun,")\n", if(!allow.dom0) c("mc[mc==0]<-m\n")),"err['",eq,"'] <- min(",if(normbym)"(" else NULL,"apply(",if(misfun) c("rep(mc,each =nr) - B,2,sum)*nc") else c("m - B,2,sumpos)*nc"),if(diag && i1==i2){c(" + useneg(sum(diag(B))*nc - ",if(misfun) "(diag(mc - B))*nc)" else c("usepos(diag(m - B))*nc)"))} else NULL,if(normbym) {if(misfun) {if(allow.dom0)")/{mc[mc==0]<-1;mc})" else ")/mc)"} else c(")/m)")} else ")","}\n"),
									"reg" = c("err['",eq,"'] <-",if(mindimreg&mindim>1)"if(nc<mindim|nr<mindim) Inf else " else NULL , if(normbym)"("else NULL, "sum(er)*nc+ sum(ec)*nr - sum(pmin(rep(er,times=length(ec)),rep(ec,each=length(er))))",if(normbym) ")/m","\n"),
									"rre" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , "sum(er)*nc",if(normbym) "/m","\n"),
									"cre" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , "sum(ec)*nr",if(normbym) "/m","\n"),
									"rfn" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , if(normbym)"("else NULL, "sumB- sum(srm) + sumpos(m - srm)*nc",if(normbym) "/m","\n"),
									"cfn" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , if(normbym)"("else NULL, "sumB- sum(scm) + sumpos(m - scm)*nr",if(normbym) "/m","\n"),
									"avg" = c("err['",eq,"'] <-", "max(0,av*nr*nc-sumB)", if(normbym) "/m", "\n"),
									"dnc" = c("err['",eq,"'] <-", "0","\n")
								))

								if(approach=="hom") {
									if(is.null(BLOCKS.CV)) {
										#cat("i1 = ", i1,", i2 = ",i2,",diag = ",diag,", regDiagSep = ", regDiagSep,", diag && i1==i2 && regDiagSep = ", diag && i1==12 && regDiagSep, "\n",sep="") 
										fun<-c(fun, switch(EXPR=eq,
											"null" = c("err['",eq,"'] <-",if(diag && i1==i2)"if(nr>=2) dev.cv(na.omit(as.vector(B2)),cv=0) + dev(diag(B)) else dev.cv(B,cv=0)" else "dev.cv(B,cv=0)","\n"),
											"com" = c("err['",eq,"'] <-",if(diag && i1==i2)"if(nr>=2) dev(na.omit(as.vector(B2))) + dev(diag(B)) else 0" else "dev(B)","\n"),
											"rdo1" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL ,if(diag && i1==i2)"min(","{tmp<-apply(B,1,cv);tmp<-idev(tmp);min(apply(B,1,dev)[which(tmp==max(tmp))])}*nc",if(diag && i1==i2)",min((apply(B2,1,function(x)dev(na.omit(x)))+idev(diag(B)))[which(tmp==max(tmp))])*nc)","\n"),
											"cdo1" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL ,if(diag && i1==i2)"min(","{tmp<-apply(B,2,cv);tmp<-idev(tmp);min(apply(B,2,dev)[which(tmp==max(tmp))])}*nr",if(diag && i1==i2)",min((apply(B2,2,function(x)dev(na.omit(x)))+idev(diag(B)))[which(tmp==max(tmp))])*nr)","\n"),
											"rdo2" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL ,if(diag && i1==i2)"min(","min(apply(B,1,dev))*nc",if(diag && i1==i2)",min((apply(B2,1,function(x)dev(na.omit(x)))+idev(diag(B))))*nc)","\n"),
											"cdo2" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL ,if(diag && i1==i2)"min(","min(apply(B,2,dev))*nr",if(diag && i1==i2)",min((apply(B2,2,function(x)dev(na.omit(x)))+idev(diag(B))))*nr)","\n"),
											"rdo3" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL ,if(diag && i1==i2)"min(","{tmp<-apply(B,1,cv);min(apply(B,1,dev)[which(tmp==max(tmp))])}*nc",if(diag && i1==i2)",min((apply(B2,1,function(x)dev(na.omit(x)))+idev(diag(B)))[which(tmp==max(tmp))])*nc)","\n"),
											"cdo3" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL ,if(diag && i1==i2)"min(","{tmp<-apply(B,2,cv);min(apply(B,2,dev)[which(tmp==max(tmp))])}*nr",if(diag && i1==i2)",min((apply(B2,2,function(x)dev(na.omit(x)))+idev(diag(B)))[which(tmp==max(tmp))])*nr)","\n"),
											"rfn" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL ," {tmp<-apply(B,1,max);dev(tmp)+dev.cv(B,cv=0)-dev.cv(tmp,cv=0)}","\n"),
											"cfn" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL ," {tmp<-apply(B,2,max);dev(tmp)+dev.cv(B,cv=0)-dev.cv(tmp,cv=0)}","\n"),
											"reg" = c("err['",eq,"'] <-",if(mindimreg&mindim>1)"if(nc<mindim|nr<mindim) Inf else " else NULL,if(!(diag && i1==i2 && regDiagSep)){
													 c("max(dev(apply(B,1,",FUN,"))*nc,dev(apply(B,2,",FUN,"))*nr)" )
												} else c("max(dev(apply(Bd,1,",FUN,"))*(nc-1),dev(apply(Bd,2,",FUN,"))*(nr-1))+ss(dev(diag(B)))" ),"\n"),
											"rre" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , if(!(diag && i1==i2 && regDiagSep)){
													c("dev(apply(B,1,",FUN,"))*nc")
												}else c("dev(apply(Bd,1,",FUN,"))*(nc-1)+dev(diag(B))"),"\n") ,
											"cre" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , if(!(diag && i1==i2 && regDiagSep)){
													c("dev(apply(B,2,",FUN,"))*nr")
												} else c("dev(apply(Bd,2,",FUN,"))*(nr-1) + dev(diag(B))"),"\n"),
											"dnc" = c("err['",eq,"'] <-", "0","\n")
										))
									}else{
										UseNotNull<-!is.null(CV.use)
										fun<-c(fun, switch(EXPR=eq,
											"null" = c("err['",eq,"'] <-",if(diag && i1==i2)"if(nr>=2) dev.cv(na.omit(as.vector(B2)),cv=0) + dev(diag(B)) else dev.cv(B,cv=0)" else "dev.cv(B,cv=0)","\n"),
											"com" = c("err['",eq,"'] <-",if(diag && i1==i2) c("if(nr>=2) dev.cv(na.omit(as.vector(B2)),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull)c(", use=CV.use[",iEq,",",i1,",",i2,"]"),") + dev(diag(B)) else dev.cv(B,cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull)c(", use=CV.use[",iEq,",",i1,",",i2,"]"),")") else c("dev.cv(B,cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull)c(", use=CV.use[",iEq,",",i1,",",i2,"]"),")"),"\n"),					
											"reg" = c("err['",eq,"'] <-",if(mindimreg&mindim>1)"if(nc<mindim|nr<mindim) Inf else " else NULL,if(!(diag && i1==i2 && regDiagSep)){
													c("max(dev.cv(apply(B,1,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull) c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*nc,dev.cv(apply(B,2,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull)c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*nr)" ,"\n")
												} else c("max(dev.cv(apply(Bd,1,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull) c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*(nc-1),dev.cv(apply(Bd,2,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull)c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*(nr-1))+ss(dev(diag(B)))" ,"\n")),
											"rre" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , if(!(diag && i1==i2 && regDiagSep)){
													c("dev.cv(apply(B,1,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull) c(", use=CV.use[",iEq,",",i1,",",i2,"]") else NULL,")*nc","\n")
												} else c("dev.cv(apply(Bd,1,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull) c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*(nc-1)+dev(diag(B))","\n")),
											"cre" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , if(!(diag && i1==i2 && regDiagSep)){
													c("dev.cv(apply(B,2,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull) c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*nr","\n")
												} else c("dev.cv(apply(Bd,2,",FUN,"),cv=BLOCKS.CV[",iEq,",",i1,",",i2,"]",if(UseNotNull) c(", use=CV.use[",iEq,",",i1,",",i2,"]")else NULL,")*(nr-1)+dev(diag(B))","\n")) ,
											"dnc" = c("err['",eq,"'] <-", "0","\n")
										))
									}
								}
							}
						}

						if(diag && i1==i2 && any(c("null","com","rdo","cdo") %in% BLOCKS[,i1,i2])){
							if(approach=="bin"){
								fun<-c(fun,"if(nr>=2){\n")
								fun<-c(fun,"sd<-sum(diag(B))\n")
								fun<-c(fun,"d.err.null<-(nr - 2*sd)\n")
								for(eq in BLOCKS[,i1,i2]){
									if(!is.na(eq)){
										fun<-c(fun,switch(EXPR=eq,
											"null" = c("err['",eq,"'] <-", "err['",eq,"']"," + min(0,d.err.null)","\n"),
											"com" = c("err['",eq,"'] <-", "err['",eq,"']"," + min(0,-d.err.null)","\n"),
											"rdo" = c("err['",eq,"'] <-",if(mindim>1)"if(nr<mindim) Inf else " else NULL , "err['",eq,"']","  - ifelse(sd==0,nc,0)","\n"),
											"cdo" = c("err['",eq,"'] <-",if(mindim>1)"if(nc<mindim) Inf else " else NULL , "err['",eq,"']","  - ifelse(sd==0,nr,0)","\n")
										))
									}
								}
								fun<-c(fun,"}\n")
							}

							if(any(c("null","com") %in% BLOCKS[,i1,i2])){
								if(approach=="val"){
									fun<-c(fun,"if(nr>=2){\n")
									fun<-c(fun,"d.err.com<- (sum(diag(B))-sumpos(m-diag(B)))\n")
									for(eq in BLOCKS[,i1,i2]){
										if(!is.na(eq)){
											fun<-c(fun,switch(EXPR=eq,
												"null" = c("err['",eq,"'] <-", "err['",eq,"']"," + min(0,-d.err.com)", if(normbym)"/m","\n"),
												"com" = c("err['",eq,"'] <-", "err['",eq,"']","  + min(0,d.err.com)",  if(normbym)"/m", "\n")
											))
										}
									}
									fun<-c(fun,"}\n")
								}
							}
						}




						if(use.weights)	fun<-c(fun,"err<-err*blockWeights",if(nr>1)c("[[",inr,"]]") else NULL,"[names(err)]\n")
						if(norm) fun<-c(fun,"err<-err/nr/nc\n")
						if(approach=="val"&& normbym) {
							fun<-c(fun,"} else {\nerr<-numeric()\nerr[",paste("c('",paste(BLOCKS[,i1,i2],sep="",collapse="','"),"')",sep=""),"]<-",if(allow.max0) "0" else "1" ,"\nif('null' %in% ",paste("c('",paste(BLOCKS[,i1,i2],sep="",collapse="','"),"')",sep=""),") err['null']<-0\n",if(!norm)"err<-err*nr*nc\n" else NULL,"}\n")
						}


						if(save.err.v && !sameModel) {
							fun<-c(fun,"ERR.V",if(nr>1)c("[[",inr,"]]") else NULL,"[,",i1,",",i2,"]<- err[BLOCKS",if(nr>1)c("[[",inr,"]]") else NULL,"[,",i1,",",i2,"]]\n")
						}
						
						fun<-c(fun,"if(length(err)==0) err['non']<-Inf \n")

						if(sameModel){
							fun<-c(fun,"ERR.V[,",i1,",",i2,",",inr,"]<- err[BLOCKS",if(nr>1)c("[[",inr,"]]") else NULL,"[,",i1,",",i2,"]]\n")
							if(!directed) fun<-c(fun,"ERR.V[,",i2,",",i1,",",inr,"]<- err[BLOCKS",if(nr>1)c("[[",inr,"]]") else NULL,"[,",i2,",",i1,"]]\n")
						}else{
							fun<-c(fun,"E[",i1,",",i2,if(nr>1)c(",",inr)else NULL,"]<-min(err)\n")
							fun<-c(fun,"IM[",i1,",",i2,if(nr>1)c(",",inr)else NULL,"]<-names(which(err==min(err))[1])\n")
						}


						if(changeT){
							fun2<-c(fun2,fun[-(1:length.fun)],"}\n")
							length.fun<-length(fun)
						}


					}

				}
			}
		}


		if(!directed && !sameModel){
			if(nr==1){
				fun<-c(fun,"E[upper.tri(E)]<-t(E)[upper.tri(E)]\nIM[upper.tri(IM)]<-t(IM)[upper.tri(IM)]\n")
			}else{
				for(inr in 1:nr){
				fun<-c(fun,"E[,,",inr,"][upper.tri(E[,,",inr,"])]<-t(E[,,",inr,"])[upper.tri(E[,,",inr,"])]\nIM[,,",inr,"][upper.tri(IM[,,",inr,"])]<-t(IM[,,",inr,"])[upper.tri(IM[,,",inr,"])]\n")
				}
			}
		}

		if(sameModel){
			fun<-c(fun,"allE<-apply(ERR.V,c(1,2,3),function(x)sum(x*c(",paste(relWeights,collapse=","),")))
			oneE<-apply(allE,c(2,3),min,na.rm=TRUE)
			oneEarray<-array(NA,dim=c(",n.types[1],",",k[1],",",k[2],"))\n")
			for(idim in 1:n.types[1]) fun<-c(fun,"oneEarray[",idim,",,]<-oneE\n")
			fun<-c(fun,"oneIM<-apply(allE==oneEarray,c(2,3),function(x)which(x)[1])\n")
			fun<-c(fun,"for(inr in 1:",nr,"){
				for(i1 in 1:",k[1],"){
					for(i2 in 1:",k[2],"){
						E[i1,i2,inr]<-ERR.V[oneIM[i1,i2],i1,i2,inr]
						IM[i1,i2,inr]<-BLOCKS[[inr]][oneIM[i1,i2],i1,i2]
					}
				}
			}\n")
			if(!is.null(positionWeights))fun<-c(fun,"oneE<-oneE*positionWeights\n")
			fun<-c(fun,"totalErr<-sum(oneE)\n")
		} else 	{
			if(nr>1){
				fun<-c(fun,"oneE<-apply(E,c(1,2),function(x)sum(x*c(",paste(relWeights,collapse=","),")))\n")
				if(!is.null(positionWeights))fun<-c(fun,"oneE<-oneE*positionWeights\n")
				fun<-c(fun,"totalErr<-sum(oneE)\n")
			}else {
				if(!is.null(positionWeights))fun<-c(fun,"E<-E*positionWeights\n")
				fun<-c(fun,"totalErr<-sum(E)\n")
			}
		}
		
		if(normMto2rel&sameModel){
			fun<-c(fun,"oneIM<-IM[,,1]\n")
			if(normMto2relRegToRreCre)fun<-c(fun,"oneIM[oneIM=='rre']<-'reg'\n")
		}
		
		fun<-c(fun,"return(list(err=totalErr,clu=clu,E=E,IM=",if(normMto2rel&sameModel)"oneIM" else "IM",",approach=approach,normMto2rel=",normMto2rel,if(sameModel){if(normMto2rel) ",completeIM = IM" else ",oneIM=oneIM"},if(nr>1)",oneE=oneE" else NULL, if(save.err.v)",ERR.V=ERR.V" else NULL ,"))\n")
		if(changeT){
			fun2<-c(fun2,fun[-(1:length.fun)],"}\n")
			fun<-c(fun,"}\n")
#			cat(file="fun1.R",fun,sep="")
#			cat(file="fun2.R",fun2,sep="")
			return(list(fun1= parse(text=paste(fun,sep="",collapse="")),fun2 = parse(text=paste(fun2,sep="",collapse=""))))
		} else {
			fun<-c(fun,"}\n")
#			cat(file="fun1.R",fun,sep="")
			return(parse(text=paste(fun,sep="",collapse="")))
		}
	}
}
