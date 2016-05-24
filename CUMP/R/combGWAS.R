combGWAS<-function(project="mv",traitlist,traitfile,comb_method=c("z"),betasign=rep(1,length(traitlist)),snpid,beta=NULL,SE=NULL,Z=NULL,coded_all,AF_coded_all,n_total=NULL,pvalue=NULL){
	trait<-paste(traitlist,collapse="_",sep="")
	print(traitlist)
	
	if (is.null(snpid) | is.null(AF_coded_all) | is.null(coded_all) | ((is.null(beta) | is.null(SE)) & is.null(Z))){
		stop("Missing necessary columns!")
	}else if (!is.null(beta) & !is.null(SE) & !is.null(n_total)){
		headername<-c(snpid,beta,SE,AF_coded_all,coded_all,n_total)
	}else if (!is.null(beta) & !is.null(SE) & is.null(n_total)){
		headername<-c(snpid,beta,SE,AF_coded_all,coded_all)
	}else if ((is.null(beta) | is.null(SE)) & !is.null(n_total)){
		headername<-c(snpid,Z,AF_coded_all,coded_all,n_total)
	}else if ((is.null(beta) | is.null(SE)) & is.null(n_total)){
		headername<-c(snpid,Z,AF_coded_all,coded_all)
	}	 

	if ("beta"%in%comb_method & (is.null(beta)|is.null(SE))){
		stop("beta or SE can not be missing for the beta method!")
	}


	n<-length(traitfile)
	data<-list()
	for (i in 1:n){
		if (length(scan(file=traitfile[i], what="character", nlines=1, sep=""))>1) which.sep <- "" else
		if (length(scan(file=traitfile[i], what="character", nlines=1, sep="\t"))>1) which.sep <- "" else 
   		if (length(scan(file=traitfile[i], what="character", nlines=1, sep=" "))>1) which.sep <- " " else 
   		if (length(scan(file=traitfile[i], what="character", nlines=1, sep=","))>1) which.sep <- "," else
		if (length(scan(file=traitfile[i], what="character", nlines=1, sep=";"))>1) which.sep <- ";" else
		stop(paste("Separator field not recognized for ",traitfile[i],sep=""))
		data[[i]]<-read.csv(traitfile[i],as.is=TRUE,sep=which.sep,header=TRUE)[,headername]  
		if (is.null(beta) | is.null(SE)){
			data[[i]]$BETA<-data[[i]][,Z]
			data[[i]]$se<-1
			data[[i]]<-data[[i]][,c(snpid,"BETA","se",AF_coded_all,coded_all,n_total)]

		}
	}	


	if (length(scan(file=traitfile[1], what="character", nlines=1, sep=""))>1) which.sep <- "" else
	if (length(scan(file=traitfile[1], what="character", nlines=1, sep="\t"))>1) which.sep <- "" else 
   	if (length(scan(file=traitfile[1], what="character", nlines=1, sep=" "))>1) which.sep <- " " else 
   	if (length(scan(file=traitfile[1], what="character", nlines=1, sep=","))>1) which.sep <- "," else
	if (length(scan(file=traitfile[1], what="character", nlines=1, sep=";"))>1) which.sep <- ";" else
	stop(paste("Separator field not recognized for ",traitfile[1],sep=""))
	data0<-read.csv(traitfile[1],as.is=TRUE,sep=which.sep,header=TRUE)


	if (is.null(n_total)){
		for (i in 1:n){
			names(data[[i]])<-c("SNPID",paste("beta",i,sep=""),paste("SE",i,sep=""),paste("AF_coded_all",i,sep=""),paste("coded_all",i,sep=""))
		}
	}else{
		for (i in 1:n){
			names(data[[i]])<-c("SNPID",paste("beta",i,sep=""),paste("SE",i,sep=""),paste("AF_coded_all",i,sep=""),paste("coded_all",i,sep=""),paste("N",i,sep=""))
		}
	}

	dat<-data[[1]]
	for (i in 2:n){
		if (identical(data[[1]]$SNPID,data[[i]]$SNPID)) dat <- cbind(dat,data[[i]][,-1]) else dat=merge(dat,data[[i]],by="SNPID",all=TRUE)
	}
	
	
	## flip sign of beta if coded alleles different between two datasets
        for(i in 1:n){
        	tmp<-na.omit(dat[,paste("coded_all",c(1,i),sep="")])
                if(!identical(tmp[,1],tmp[,2])) {
                   warning(paste("coded alleles differ between datasets 1 and",i,",sign of beta in",i," is reversed for those SNPs"))
		   switch<-dat[,"coded_all1"]!=dat[,paste("coded_all",i,sep="")] &!is.na(dat[,"coded_all1"]!=dat[,paste("coded_all",i,sep="")])
		   dat[switch,paste("beta",i,sep="")]<- - dat[switch,paste("beta",i,sep="")]
		   dat[switch,paste("coded_all",i,sep="")]<- dat[switch,"coded_all1"]
                }
        }



	##reverse the beta sign if needed
	for (i in 1:n){
		dat[,paste("beta",i,sep="")]<-dat[,paste("beta",i,sep="")]*betasign[i]
	}
	


	
	CombN_zchisq<-function(){
		z<-matrix(0,nrow(dat),n)
		maf<-matrix(0,nrow(dat),n)
		for (i in 1:n){
			z[,i]<-ifelse(dat[,paste("SE",i,sep="")]>0,dat[,paste("beta",i,sep="")]/dat[,paste("SE",i,sep="")],NA)
			maf[,i]<-pmin(dat[,paste("AF_coded_all",i,sep="")],1-dat[,paste("AF_coded_all",i,sep="")])
		}
	
		MAF<-rowMeans(maf)

		outcol<-names(data0)[!names(data0)%in%c(beta,SE,pvalue,Z)]
       	col.keep<-data0[,outcol]

		if (!is.null(n_total)){
			meanN<-apply(dat[,c(paste("N",1:n,sep=""))],1,mean,na.rm=TRUE)
	    		minN<-apply(dat[,c(paste("N",1:n,sep=""))],1,min,na.rm=TRUE)
	     		maxN<-apply(dat[,c(paste("N",1:n,sep=""))],1,max,na.rm=TRUE)
			Z=data.frame(dat$SNP,MAF,z,meanN,minN,maxN)
		}else{
			Z=data.frame(dat$SNP,MAF,z)
		}

		for (i in 1:n){
			names(Z)[names(Z)==paste("X",i,sep="")]<-paste("z",i,sep="")
		}

		cov.all=matrix(0,n,n)
		corr=matrix(0,n,n)
		for (i in 1:n){
              	for (j in 1:n){
                     	cov.all[i,j]=cov(Z[,i+2],Z[,j+2],use="complete.obs")
                    }
              }
		for (i in 1:n){
              	cov.all[i,i]=var(Z[,i+2],na.rm=TRUE)
	       }
		corr<-cor(Z[,3:(n+2)],use="complete.obs")
       	if(min(eigen(corr)$values)<0.01){
			stop ("The correlation matrix is nearly singular.")
		}
         
       	list(col.keep=col.keep,Z=Z,cov.all=cov.all,corr=corr)
	}






	#### For OBrien's method

	doN_z<-function(){
		dat<-CombN_zchisq()
		outfile<-paste(project,"_",trait,"_z.csv",sep="")
		print(paste("output is in ", outfile))
		if (!is.null(n_total)){
			out<-dat$Z[,c("dat.SNP",paste("z",1:n,sep=""),"meanN","minN","maxN")]
		}else{
			out<-dat$Z[,c("dat.SNP",paste("z",1:n,sep=""))]
		}
		out$beta <- t(matrix(rep(1,n),1,n)%*%solve(dat$cov.all)%*%t(dat$Z[,c(paste("z",1:n,sep=""))]))
		var <-matrix(rep(1,n),1,n)%*%solve(dat$cov.all)%*%matrix(rep(1,n),n,1)
		out$SE<-ifelse(!is.na(out$beta),sqrt(var),NA)
		out$Z.comb <- out$beta/out$SE
		out$pval<-2*pnorm(abs(out$Z.comb),lower.tail=FALSE)
		for (i in 1:n){
			out[,paste("p",i,sep="")]<-2*pnorm(abs(out[,paste("z",i,sep="")]),lower.tail=FALSE)
		}
		out<-merge(dat$col.keep,out,by.x=snpid, by.y="dat.SNP",all=TRUE)

		out1<-out[order(out$Z.comb^2, decreasing=TRUE),]
		##names(out1)[1]<-"SNPID"
		write.table(out1, outfile, sep=",",quote=FALSE,na="",row.names=FALSE)
		write.table(dat$corr,"correlation_df.out",sep=" ",quote=FALSE,col.names=traitlist,row.names=traitlist)
	}


	#### For chi-square method

	doN_chisq<-function(){
		dat<-CombN_zchisq()
		outfile<-paste(project,"_",trait,"_chisq.csv",sep="")
		print(paste("output is in ", outfile))
		if (!is.null(n_total)){
			out<-dat$Z[,c("dat.SNP",paste("z",1:n,sep=""),"meanN","minN","maxN")]
		}else{
			out<-dat$Z[,c("dat.SNP",paste("z",1:n,sep=""))]
		}
	
		comp_cs<-function(x){
			cs<-matrix(unlist(x)[1:n],nrow=1)%*%solve(dat$cov.all)%*%matrix(unlist(x)[1:n])
			cs
		}

		out$chisq.comb <- apply(out[,c(paste("z",1:n,sep=""))],1,comp_cs)
		out$pval<-pchisq(out$chisq.comb,df=n,lower.tail=FALSE)
		for (i in 1:n){
			out[,paste("p",i,sep="")]<-2*pnorm(abs(out[,paste("z",i,sep="")]),lower.tail=FALSE)
		}
		out<-merge(dat$col.keep,out,by.x=snpid, by.y="dat.SNP",all=TRUE)

		out1<-out[order(out$chisq.comb, decreasing=TRUE),]
		##names(out1)[1]<-"SNPID"
		write.table(out1, outfile, sep=",",quote=FALSE,na="",row.names=FALSE)
		write.table(dat$corr,"correlation_df.out",sep=" ",quote=FALSE,col.names=traitlist,row.names=traitlist)
	}




	#### For sum of square method

	doN_sumsq<-function(){
		dat<-CombN_zchisq()
		outfile<-paste(project,"_",trait,"_sumsq.csv",sep="")
		print(paste("output is in ", outfile))
		if (!is.null(n_total)){
			out<-dat$Z[,c("dat.SNP",paste("z",1:n,sep=""),"meanN","minN","maxN")]
		}else{
			out<-dat$Z[,c("dat.SNP",paste("z",1:n,sep=""))]
		}

	
		comp_cs<-function(x){
			cs<-matrix(unlist(x)[1:n],nrow=1)%*%matrix(unlist(x)[1:n])
			cs
		}

		out$chisq.comb <- apply(out[,c(paste("z",1:n,sep=""))],1,comp_cs)
		c<-eigen(dat$cov.all)$values
		a<-sum(c^3)/sum(c^2)
		b<-sum(c)-((sum(c^2))^2)/sum(c^3)
		d<-((sum(c^2))^3)/((sum(c^3))^2)
		out$pval<-pchisq((out$chisq.comb-b)/a,df=d,lower.tail=FALSE)
		for (i in 1:n){
			out[,paste("p",i,sep="")]<-2*pnorm(abs(out[,paste("z",i,sep="")]),lower.tail=FALSE)
		}
		out<-merge(dat$col.keep,out,by.x=snpid, by.y="dat.SNP",all=TRUE)

		out1<-out[order(out$chisq.comb, decreasing=TRUE),]
		##names(out1)[1]<-"SNPID"
		write.table(out1, outfile, sep=",",quote=FALSE,na="",row.names=FALSE)
		write.table(dat$corr,"correlation_df.out",sep=" ",quote=FALSE,col.names=traitlist,row.names=traitlist)
		write(paste("sumsq df=",d,sep=""),"correlation_df.out",append=TRUE)
	}




	#### For beta method

	CombN_beta<-function(){
		beta<-matrix(0,nrow(dat),n)
		se<-matrix(0,nrow(dat),n)
		maf<-matrix(0,nrow(dat),n)
		for (i in 1:n){
			beta[,i]<-dat[,paste("beta",i,sep="")]
			se[,i]<-dat[,paste("SE",i,sep="")]
			maf[,i]<-pmin(dat[,paste("AF_coded_all",i,sep="")],1-dat[,paste("AF_coded_all",i,sep="")])
		}
		
		MAF<-rowMeans(maf)

		outcol<-names(data0)[!names(data0)%in%c(beta,SE,pvalue)]
       	col.keep<-data0[,outcol]
	
		if (!is.null(n_total)){
			meanN<-apply(dat[,c(paste("N",1:n,sep=""))],1,mean,na.rm=TRUE)
       		minN<-apply(dat[,c(paste("N",1:n,sep=""))],1,min,na.rm=TRUE)
       		maxN<-apply(dat[,c(paste("N",1:n,sep=""))],1,max,na.rm=TRUE)
			Z=data.frame(dat$SNP,MAF,beta,se,meanN,minN,maxN)
		}else{
			Z=data.frame(dat$SNP,MAF,beta,se)
		}
		

		for (i in 1:n){
			names(Z)[names(Z)==paste("X",i,sep="")]<-paste("beta",i,sep="")
			names(Z)[names(Z)==paste("X",i,".1",sep="")]<-paste("SE",i,sep="")
		}

		cov.all=matrix(0,n,n)
		corr=matrix(0,n,n)
              for (i in 1:n){
              	for (j in 1:n){
                     	cov.all[i,j]=cor(Z[,i+2],Z[,j+2],use="complete.obs")
              	}
              }
             	for (i in 1:n){
              	cov.all[i,i]=var(Z[,i+2],na.rm=TRUE)
	       }
		corr<-cor(Z[,3:(n+2)],use="complete.obs")
		if(min(eigen(corr)$values)<0.01){
			stop ("The correlation matrix is nearly singular.")
		}

      		list(col.keep=col.keep,Z=Z,cov.all=cov.all,corr=corr)
	}


	doN_beta<-function(){
		n<-length(traitfile)
		dat<-CombN_beta()
		outfile<-paste(project,"_",trait,"_beta.csv",sep="")
		print(paste("output is in ", outfile))
		if (!is.null(n_total)){
			out<-dat$Z[,c("dat.SNP",paste("beta",1:n,sep=""),paste("SE",1:n,sep=""),"meanN","minN","maxN")]
		}else{
			out<-dat$Z[,c("dat.SNP",paste("beta",1:n,sep=""),paste("SE",1:n,sep=""))]
		}

		newse<-function(x){
			cov.mat<-diag(unlist(x)[(n+1):(2*n)])%*%dat$cov.all%*%diag(unlist(x)[(n+1):(2*n)])
			if (any(is.na(unlist(x))) | any(unlist(x)[(n+1):(2*n)]==0)|det(cov.mat)==0){
				beta<-NA
				SE<-NA
				Z.comb <- NA
				pval<-NA
			}else{			
				beta<-matrix(rep(1,n),1,n)%*%solve(cov.mat)%*%matrix(unlist(x)[1:n])
				var<-matrix(rep(1,n),1,n)%*%solve(cov.mat)%*%matrix(rep(1,n),n,1)
				SE<-ifelse(!is.na(beta),sqrt(var),NA)
				Z.comb <- beta/SE
				pval<-2*pnorm(abs(Z.comb),lower.tail=FALSE)
			}
			list(beta=beta,SE=SE,Z.comb=Z.comb,pval=pval)
		}

		out2<-apply(out[,c(paste("beta",1:n,sep=""),paste("SE",1:n,sep=""))],1,newse)
		out3<-unlist(out2)
		out4<-matrix(out3,nrow(out),4,byrow=TRUE,dimnames=list(NULL,c("beta","SE","Z.comb","pval")))
		out5<-as.data.frame(out4)
		out<-cbind(out,out5)
		for (i in 1:n){
			out[,paste("p",i,sep="")]<-2*pnorm(abs(ifelse(dat$Z[,paste("SE",i,sep="")]>0,dat$Z[,paste("beta",i,sep="")]/dat$Z[,paste("SE",i,sep="")],NA)),lower.tail=FALSE)
		}

		out<-merge(dat$col.keep,out,by.x=snpid, by.y="dat.SNP",all=TRUE)
		out<-out[,!names(out)%in%paste("SE",1:n,sep="")]
		
		out1<-out[order(out$Z.comb^2, decreasing=TRUE),]
		##names(out1)[1]<-"SNPID"
		write.table(out1, outfile, sep=",",quote=FALSE,na="",row.names=FALSE)
		write.table(dat$corr,"correlation_df.out",sep=" ",quote=FALSE,col.names=traitlist,row.names=traitlist)
	}


	if (FALSE%in%(comb_method%in%c("z","beta","chisq","sumsq"))){
		stop("Please check the method name!")
	}
	
	if ("beta"%in%comb_method){
		doN_beta()
	}

	if ("z"%in%comb_method){
		doN_z()
	}

	if ("chisq"%in%comb_method){
		doN_chisq()
	}

	if ("sumsq"%in%comb_method){
		doN_sumsq()
	}


}









