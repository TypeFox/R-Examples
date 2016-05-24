getSampleCCF <- function(id,dd,VCFdir,thr.coverage=10,upper.cov=0.99,nc=10,tc=11,AD=3,filter=TRUE,TCGA=TRUE){
    getID <- function(sID,num.res=3){
        mm=c()
        for(id in sID){
            tmp=unlist(strsplit(id,'-'))
            if(length(tmp)==1){
                tmp=unlist(strsplit(id,'\\.'))
            }
            ll='TCGA'
            for(j in 2:num.res){
                ll=paste(ll,tmp[j],sep='-')
            }
            mm=c(mm,ll)
        }
        return(mm)
    }
	vcffiles<-dir(VCFdir,full.names=TRUE)
	cat('loading vcf file',id,'\n')
	if(TCGA)vv.id<-grep(getID(id),vcffiles) else vv.id<-grep(id,vcffiles)
	if(length(vv.id)==0){
		cat('No VCF file!')
		return(NA)
	}
	vcf<-read.table(vcffiles[vv.id],header=FALSE,stringsAsFactors=FALSE)
	vv.a<-which(vcf[,tc]!='./.'&vcf[,nc]!='./.')
	vcf<-vcf[vv.a,]
	y<-unlist(lapply(strsplit(vcf[,tc],':'),function(x)x[[AD]]))
	y<-strsplit(y,',')
	y0<-unlist(lapply(strsplit(vcf[,nc],':'),function(x)x[[AD]]))
	y0<-strsplit(y0,',')
	vv.bi<-which(unlist(lapply(y,length))==2&unlist(lapply(y0,length))==2)
	y<-y[vv.bi]
	y0<-y0[vv.bi]
	vcf<-vcf[vv.bi,]
	a<-as.numeric(unlist(lapply(y,function(x)x[[1]])))
	b<-as.numeric(unlist(lapply(y,function(x)x[[2]])))
	s<-a+b
	a0<-as.numeric(unlist(lapply(y0,function(x)x[[1]])))
	b0<-as.numeric(unlist(lapply(y0,function(x)x[[2]])))
	s0<-a0+b0
	QUAL<-vcf[,6]
	thr.upper<-quantile(s,upper.cov)
	gt.n<-unlist(lapply(strsplit(vcf[,nc],':'),function(x)x[[1]]))
	gt.t<-unlist(lapply(strsplit(vcf[,tc],':'),function(x)x[[1]]))
	vv.somatic<-which(gt.n=='0/0'&gt.t!='0/0'|gt.n=='1/1'&gt.t!='1/1')	## somatic mutations
	vv.cov<-which(s>=thr.coverage&s0>=thr.coverage)	## passed coverage requirement
	if(filter)vv.pass<-which(vcf[,7]=='PASS') else vv.pass<-1:nrow(vcf)	## passed filtered mutations
	#vv.pass=which(b>=s.coverage&QUAL>=10)
	sam.dd<-dd[rownames(dd)==id,]
	MutInfo<-c()
	n<-nrow(sam.dd)
	vv.int<-intersect(vv.somatic,vv.cov)
	vv.int<-intersect(vv.int,vv.pass)
    cc<-seq(0.02,1,by=0.01)
    computeSD <- function(N,S,f){
        M1list<-c()
        M2list<-c()
        MLElist<-c()
        for(ii in 1:length(N)){
            PF<-sum(dbinom(S[ii],N[ii],f),na.rm=TRUE)
            M1<-sum(dbinom(S[ii],N[ii],f)*cc,na.rm=TRUE)/PF
            M2<-sum(dbinom(S[ii],N[ii],f)*cc^2,na.rm=TRUE)/PF
            M1list<-c(M1list,M1)
            M2list<-c(M2list,M2)
            MLElist<-c(MLElist,cc[which.max(dbinom(S[ii],N[ii],f))])
        }
        return(list(M1=MLElist,SD=sqrt(M2list-M1list^2)))
    }
	tmp.vcf<-vcf[vv.int,]
	for(k in 1:n){
		vv<-which(tmp.vcf[,1]==sam.dd[k,1]&tmp.vcf[,2]>=sam.dd[k,2]&tmp.vcf[,2]<=sam.dd[k,3])
		if(length(vv)==0)next
		nb<-sam.dd[k,9]
		nt<-sam.dd[k,10]
		tmp.b<-b[vv.int[vv]]
		tmp.s<-s[vv.int[vv]]
		tmp.a<-tmp.s-tmp.b
		tmp.b0<-b0[vv.int[vv]]
		tmp.s0<-s0[vv.int[vv]]
		nv<-length(vv)
        if(is.na(sam.dd[k,8])){
            INFO<-paste(rep(id,nv),tmp.b,tmp.s,tmp.b0,tmp.s0,NA,NA,rep(NA,nv),rep(NA,nv),rep(NA,nv),rep(NA,nv),sep=';')
            tmp.dd<-tmp.vcf[vv,]
            tmp.dd[,8]<-INFO
            MutInfo<-rbind(MutInfo,tmp.dd)
            next
        }
        p<-sam.dd[k,8]
		if(p>1)p<- 1
		if(nb==1 & nt==2){
			## unchanged genome
            ff=cc/2
            Ms=computeSD(tmp.s,tmp.b,ff)
			CCF<-Ms$M1
            SDs<-Ms$SD
			IsEarly<-rep('C',nv)
		}
		else if(nt==1){
			nc<-nt*p+2*(1-p)
            ff<-cc/nc
            Ms<-computeSD(tmp.s,tmp.b,ff)
			CCF<-Ms$M1
            SDs<-Ms$SD
            fh.ea<-(p*(nt-nb)+1-p)/nc        ## early, major allele, high limit
            fl.ea<-(p*(nt-nb))/nc    ## early, major allele, low limit
            fh.t<-p/nc       ## late, in tumor, high limit
            fh.e<-(1-p)/nc   ## in euploid cell, high limit
            pEarly.a<-pbeta(rep(fh.ea,nv),tmp.b+1,tmp.a+1)-pbeta(rep(fl.ea,nv),tmp.b+1,tmp.a+1)
            pLate<-pbeta(rep(fh.t,nv),tmp.b+1,tmp.a+1)
            pEuploid<-pbeta(rep(fh.e,nv),tmp.b+1,tmp.a+1)
            Ptot<-pEarly.a+pLate+pEuploid
            cp.A<-pEarly.a/Ptot
            cp.CD<- 1-cp.A
            cp.C<-pLate/Ptot
            cp.D<-pEuploid/Ptot
			cp.AC<- 1-cp.D
			cp.AD<- 1-cp.C
            vv.A<-which(cp.A>=0.95)
            vv.CD<-which(cp.CD>=0.95&cp.C<0.95&cp.D<0.95)
            vv.C<-which(cp.C>=0.95)
            vv.D<-which(cp.D>=0.95)
			vv.AC<-which(cp.AC>=0.95&cp.A<0.95&cp.C<0.95)
			vv.AD<-which(cp.AD>=0.95&cp.A<0.95&cp.D<0.95)
			IsEarly<-rep('A1/B/C',nv)
			IsEarly[vv.A]<-'A1'
			IsEarly[vv.C]<-'B'
			IsEarly[vv.D]<-'C'
			IsEarly[vv.AC]<-'A1/B'
			IsEarly[vv.AD]<-'A1/C'
			IsEarly[vv.CD]<-'B/C'
		}
		else if(nb==0|nt==2*nb){
			nc<-nt*p+2*(1-p)
			fh.ea<-(p*(nt-nb)+1-p)/nc	## early, major allele, high limit
			fl.ea<-(p*(nt-nb))/nc	## early, major allele, low limit
			fh.t<-p/nc	## late, in tumor, high limit
			fh.e<-(1-p)/nc	## in euploid cell, high limit
			pEarly.a<-pbeta(rep(fh.ea,nv),tmp.b+1,tmp.a+1)-pbeta(rep(fl.ea,nv),tmp.b+1,tmp.a+1)
			pLate<-pbeta(rep(fh.t,nv),tmp.b+1,tmp.a+1)
			pEuploid<-pbeta(rep(fh.e,nv),tmp.b+1,tmp.a+1)
			Ptot<-pEarly.a+pLate+pEuploid
			cpEarly.a<-pEarly.a/Ptot
			cpLate.eup<-1-cpEarly.a
			cpLate<-pLate/Ptot
			cpEup<-pEuploid/Ptot
			vv.early.a<-which(cpEarly.a>=0.95)
			vv.late.eup<-which(cpLate.eup>=0.95)
			vv.late<-which(cpLate>=0.95)
			vv.Eup<-which(cpEup>=0.95)
			z<-tmp.b/tmp.s
			CCF<-rep(NA,nv)
            SDs<-rep(NA,nv)
            ff.early=(cc-p+(nt-nb)*p)/nc
            Ms.early=computeSD(tmp.s[vv.early.a],tmp.b[vv.early.a],ff.early)
			CCF[vv.early.a]<-Ms.early$M1
            SDs[vv.early.a]<-Ms.early$SD
            ff.late=cc/nc
            Ms.late=computeSD(tmp.s[c(vv.late,vv.late.eup)],tmp.b[c(vv.late,vv.late.eup)],ff.late)
			CCF[c(vv.late,vv.late.eup)]<-Ms.late$M1
            SDs[c(vv.late,vv.late.eup)]<-Ms.late$SD
            CCF<-round(CCF,3)
			IsEarly<-rep('A1/B/C',nv)
			IsEarly[vv.early.a]<-'A1'
			IsEarly[vv.late.eup]<-'B/C'
			IsEarly[vv.late]<-'B'
			IsEarly[vv.Eup]<-'C'
		}
		else if(nb>=1 & nt>2){
			nc<-nt*p+2*(1-p)
			fh.ea<-(p*(nt-nb)+1-p)/nc	## early, major allele, high limit
			fl.ea<-(p*(nt-nb))/nc	## early, major allele, low limit
			fh.eb<-(nb*p+1-p)/nc	## early, minor allele, high limit
			fl.eb<-nb*p/nc	## early, minor allele, low limit
			fh.t<-p/nc	## late, in tumor, high limit
			fh.e<-(1-p)/nc	## euploid, high limit
			pEarly.a<-pbeta(rep(fh.ea,nv),tmp.b+1,tmp.a+1)-pbeta(rep(fl.ea,nv),tmp.b+1,tmp.a+1)
			pEarly.b<-pbeta(rep(fh.eb,nv),tmp.b+1,tmp.a+1)-pbeta(rep(fl.eb,nv),tmp.b+1,tmp.a+1)
			pLate<-pbeta(rep(fh.t,nv),tmp.b+1,tmp.a+1)
			pEuploid<-pbeta(rep(fh.e,nv),tmp.b+1,tmp.a+1)
			Ptot<-pEarly.a+pEarly.b+pLate+pEuploid
			cp.A<-pEarly.a/Ptot
			cp.B<-pEarly.b/Ptot
			cp.C<-pLate/Ptot
			cp.D<-pEuploid/Ptot
			cp.AB<- 1-cp.C-cp.D
			cp.AC<- 1-cp.B-cp.D
			cp.AD<- 1-cp.B-cp.D
			cp.BC<- 1-cp.A-cp.D
			cp.BD<- 1-cp.A-cp.C
			cp.CD<- 1-cp.A-cp.B
			cp.ABC<- 1-cp.D
			cp.ABD<- 1-cp.C
			cp.ACD<- 1-cp.B
			cp.BCD<- 1-cp.A
			vv.A<-which(cp.A>=0.95)
			vv.B<-which(cp.B>=0.95)
			vv.C<-which(cp.C>=0.95)
			vv.D<-which(cp.D>=0.95)
			vv.CD<-which(cp.CD>=0.95&cp.C<0.95&cp.D<0.95)
			vv.AB<-which(cp.AB>=0.95&cp.A<0.95&cp.B<0.95)
			vv.AC<-which(cp.AC>=0.95&cp.A<0.95&cp.C<0.95)
			vv.AD<-which(cp.AD>=0.95&cp.A<0.95&cp.D<0.95)
			vv.BC<-which(cp.BC>=0.95&cp.B<0.95&cp.C<0.95)
			vv.BD<-which(cp.BD>=0.95&cp.B<0.95&cp.D<0.95)
            vv.BCD<-which(cp.BCD>=0.95&cp.BC<0.95&cp.BD<0.95&cp.CD<0.95&cp.B<0.95&cp.C<0.95&cp.D<0.95)
			vv.ABC<-which(cp.ABC>=0.95&cp.BC<0.95&cp.AB<0.95&cp.AC<0.95&cp.B<0.95&cp.C<0.95&cp.A<0.95)
			vv.ABD<-which(cp.ABD>=0.95&cp.AB<0.95&cp.AD<0.95&cp.BD<0.95&cp.B<0.95&cp.D<0.95&cp.A<0.95)
			vv.ACD<-which(cp.ACD>=0.95&cp.AC<0.95&cp.AD<0.95&cp.CD<0.95&cp.A<0.95&cp.D<0.95&cp.C<0.95)
			z<-tmp.b/tmp.s
			CCF<-rep(NA,nv)
            SDs<-rep(NA,nv)
            ff.A<-(cc-p+(nt-nb)*p)/nc
            Ms.A<-computeSD(tmp.s[vv.A],tmp.b[vv.A],ff.A)
            CCF[vv.A]<-Ms.A$M1
            SDs[vv.A]<-Ms.A$SD
            ff.B<-(cc-p+nb*p)/nc
            Ms.B<-computeSD(tmp.s[vv.B],tmp.b[vv.B],ff.B)
            CCF[vv.B]<-Ms.B$M1
            SDs[vv.B]<-Ms.B$SD
            ff.C<-cc/nc
            Ms.C<-computeSD(tmp.s[c(vv.C,vv.D,vv.CD)],tmp.b[c(vv.C,vv.D,vv.CD)],ff.C)
            CCF[c(vv.C,vv.D,vv.CD)]<-Ms.C$M1
            SDs[c(vv.C,vv.D,vv.CD)]<-Ms.C$SD
			if(nb==1){
                ff.BCD<-cc/nc
                Ms.BCD<-computeSD(tmp.s[c(vv.BCD,vv.BC,vv.BD)],tmp.b[c(vv.BCD,vv.BC,vv.BD)],ff.BCD)
                CCF[c(vv.BCD,vv.BC,vv.BD)]<-Ms.BCD$M1
                SDs[c(vv.BCD,vv.BC,vv.BD)]<-Ms.BCD$SD
			}
			IsEarly<-rep('A1/A2/B/C',nv)
			IsEarly[vv.A]<-'A1'
			IsEarly[vv.BCD]<-'A2/B/C'
			IsEarly[vv.CD]<-'B/C'
			IsEarly[vv.B]<-'A2'
			IsEarly[vv.C]<-'B'
			IsEarly[vv.D]<-'C'
			IsEarly[vv.AB]<-'A1/A2'
			IsEarly[vv.AC]<-'A1/B'
			IsEarly[vv.AD]<-'A1/C'
			IsEarly[vv.BC]<-'A2/B'
			IsEarly[vv.BD]<-'A2/C'
			IsEarly[vv.ACD]<-'A1/B/C'
			IsEarly[vv.ABD]<-'A1/A2/C'
			IsEarly[vv.ABC]<-'A1/A2/B'
		}
		p<-round(p,3)
		CCF<-round(CCF,3)
        SDs<-round(SDs,3)
		INFO<-paste(rep(id,nv),tmp.b,tmp.s,tmp.b0,tmp.s0,CCF,SDs,rep(p,nv),rep(nb,nv),rep(nt,nv),IsEarly,sep=';')
		#tmp.mat=cbind(as.numeric(tmp.vcf[vv,1]),as.numeric(tmp.vcf[vv,2]),tmp.vcf[vv,4],tmp.vcf[vv,5],tmp.b,tmp.s,CCF,matrix(rep(c(p,nb,nt),nv),ncol=3,byrow=T),IsEarly)
		tmp.dd<-tmp.vcf[vv,]
		tmp.dd[,8]<-INFO
		MutInfo<-rbind(MutInfo,tmp.dd)
	}
	if(is.null(MutInfo))return()
	#MutInfo=as.matrix(MutInfo)
	#rownames(MutInfo)=rep(id,nrow(MutInfo))
	return(MutInfo)
}