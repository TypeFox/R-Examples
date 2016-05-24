
linkp<-function(coef)
{
	l<-length(coef)
	out<-matrix(NA,l,l)
	for(i in 1:l)
		for (j in 1:l)
			out[i,j]<-sum(rownames(coef[[i]])%in%rownames(coef[[j]]))
	return(out)
}


modIRT<-function(coef,var=NULL,names=NULL,ltparam=TRUE,lparam=TRUE,display=TRUE,digits=2)
{
	coef<-lapply(coef,FUN=function(x) x<-as.matrix(x[,colSums(x==1)!=nrow(x)]))
	coef<-lapply(coef,FUN=function(x) x<-as.matrix(x[,colSums(x==0)!=nrow(x)]))
	itmp<-sapply(coef,ncol)
	if (length(unique(itmp))!=1) stop("Mixed model types not allowed.")
	itmp<-itmp[1]
	if (itmp==1 | itmp==2) lparam<-FALSE
	lc<-length(coef)
	if (is.null(names)) names<-paste("T",1:lc,sep="")
	mods<-list()
	for (i in 1:lc) {
		c0<-c1<-c2<-NULL
		D1<-D2<-matrix(NA,0,0)
		coefi<-coef[[i]]
		if (!is.null(var)) vari<-var[[i]]
		else vari<-NULL
		if (!is.null(vari)) if (any(is.na(vari))) vari<-NULL
		n1<-nrow(coefi)
		if (itmp==1)
			c1<-coefi[,1]
		if (itmp==2) { 
			c1<-coefi[,1]
			c2<-coefi[,2]
		}
		if (itmp==3) { 
			c0<-coefi[,1]
			c1<-coefi[,2]
			c2<-coefi[,3]
		}
		if (ltparam & itmp==1) {
			c1<--c1
		}
		if (ltparam & itmp>1) {
			D1<-rbind(diag(-1/c2),diag(c1/c2^2))
			D1<-cbind(D1,rbind(matrix(0,n1,n1),diag(1,n1)))
			c1<--c1/c2
		}
		if (lparam) {
			D2<-diag(exp(c0)/(1+exp(c0))^2)
			c0<-exp(c0)/(1+exp(c0))
		}
		if (!lparam & itmp==3) {
			D2<-diag(1,n1)
		}
		if (!is.null(vari)) {
			if ((ltparam | lparam) & itmp>1) {
				D1<-blockdiag(D2,D1)
				vari<-t(D1)%*%vari%*%D1
			}
		}
		if (itmp==3) names(c0)<-paste("Gussng",names(c0),sep=".")
		names(c1)<-paste("Dffclt",names(c1),sep=".")
		if (itmp>1) names(c2)<-paste("Dscrmn",names(c2),sep=".")
		coefi<-c(c0,c1,c2)
		if (!is.null(vari)) rownames(vari)<-colnames(vari)<-names(coefi)
		mods[[i]]<-list(coefficients=coefi,var=vari,itmp=itmp)
		if (display) {
			if (!is.null(var)) out<-cbind(coefi,diag(vari)^0.5)
			else out<-cbind(coefi,NA)
			out<-round(out,digits)
			colnames(out)<-c("value","std.err")
			cat("Form:",names[i],"\n")
			print(out)
		}
	}
	names(mods)<-names
	class(mods)<-"modIRT"
	return(mods)
}


irtp1<-function(ab,diff,discr,guess,D)
{
	elp<-exp(D*discr*(ab-diff))
	guess+(1-guess)*elp/(1+elp)
}



obj<-function(eqc,P2,ab,a1,b1,c1,met,itmp,wt,D=D){
	ifelse(itmp==1,A<-1,A<-eqc[1])
	B<-eqc[2]

	b12<-A*b1+B
	a12<-a1/A

	ni<-ncol(P2)
	P1<-matrix(NA,length(ab),ni)
	for (i in 1:ni)
		P1[,i]<-irtp1(ab,diff=b12[i],discr=a12[i],guess=c1[i],D=D)

	if (met=="Haebara") f<-0.5*sum(rowSums((P2-P1)^2)*wt)
	if (met=="Stocking-Lord") f<-0.5*sum(((rowSums(P2)-rowSums(P1))^2)*wt)
	return(f)
}


direc<-function(mod1,mod2,method="mean-mean",suff1=".1",suff2=".2",D=1,quadrature=TRUE,nq=30)
{
	if (method!="mean-mean" & method!="mean-sigma" & method!="mean-gmean" & method!="Haebara" & method!="Stocking-Lord") warning("Method not implemented.")
	name1<-names(mod1)
	name2<-names(mod2)
	forms<-paste(name1,name2,sep=".")
	mod1<-mod1[[1]]
	mod2<-mod2[[1]]
	tab1<-data.frame(value1=mod1$coef)
	tab2<-data.frame(value2=mod2$coef)
	var1<-mod1$var
	var2<-mod2$var
	itmp1<-mod1$itmp
	itmp2<-mod2$itmp
	if (itmp1!=itmp2) stop("Mixed type models not allowed.")
	itmp<-itmp1
	tab<-merge(tab1,tab2,by=0)
	if (nrow(tab)==0) {
		warning("no common items")
		out<-list(ni=0)
	}
	if (nrow(tab)>0) {
		comuni<-tab$Row.names
		taball<-merge(tab1,tab2,by=0,all=T)
		niall<-nrow(taball)/itmp
		ni<-nrow(tab)/itmp
		if (ni==0) comuni<-"None"
		tabDff<-tab[substr(tab$Row.names,1,6)=="Dffclt",]
		tabDsc<-tab[substr(tab$Row.names,1,6)=="Dscrmn",]
		tabGss<-tab[substr(tab$Row.names,1,6)=="Gussng",]
		if (nrow(tabDff)==0 & ni>0) warning("missing difficulty parameters")
		if (nrow(tabGss)!=nrow(tabDff) & itmp==3) warning("missing guessing parameters but argument itmp=3")
		if (nrow(tabDsc)!=nrow(tabDff) & itmp==2) warning("missing discrimination parameters but argument itmp=2")
		ifelse (itmp>1,a1<-tabDsc$value1,a1<-rep(1,ni))
		ifelse (itmp>1,a2<-tabDsc$value2,a2<-rep(1,ni))
		b1<-tabDff$value1
		b2<-tabDff$value2
		ifelse (itmp==3,c1<-tabGss$value1,c1<-rep(0,ni))
		ifelse (itmp==3,c2<-tabGss$value2,c2<-rep(0,ni))

		if (method=="mean-sigma" & itmp>1) A<-sd(b2)/sd(b1)
		if (method=="mean-sigma" & itmp==1) A<-1
		if (method=="mean-mean")  A<-mean(a1)/mean(a2)
		if (method=="mean-gmean") A<-exp(sum(log(a1/a2)))^(1/ni)
		if (method=="mean-sigma" | method=="mean-mean" | method=="mean-gmean") 
			B<-mean(b2)-A*mean(b1)
		if (method=="Haebara" | method=="Stocking-Lord") {
			if (quadrature) {
				  gq<-gauss.quad.prob(nq,dist="normal")
				  ab<-gq$nodes
				  wt<-gq$weights
			}
			else {
				ab<-seq(-4,4,l=40)
				wt<-rep(1,40)
			}
			P2<-matrix(NA,length(ab),ni)
			for (i in 1:ni)
				P2[,i]<-irtp1(ab,diff=b2[i],discr=a2[i],guess=c2[i],D=D)
			par<-nlminb(start=c(1,0),objective=obj,P2=P2,ab=ab,a1=a1,b1=b1,c1=c1,met=method,itmp=itmp,wt=wt,D=D)$par

			A<-par[1]
			B<-par[2]
			b12<-A*b1+B
			a12<-a1/A
			P1<-matrix(NA,length(ab),ni)
			for (i in 1:ni)
				P1[,i]<-irtp1(ab,diff=b12[i],discr=a12[i],guess=c1[i],D=D)
		}
		if (!is.null(var1) & !is.null(var2)) {
			if (method=="mean-sigma") {
				partialA_b2<-A*sd(b2)^(-2)*(b2-mean(b2))/ni
				partialA_b1<--A*sd(b1)^(-2)*(b1-mean(b1))/ni
				if (itmp==1) {partialA_b2<-rep(0,ni)
					partialA_b1<-rep(0,ni)}
				partialB_b2<-1/ni-partialA_b2*mean(b1)
				partialB_b1<- -partialA_b1*mean(b1)-A/ni
				partialA_a2<-rep(0,ni)
				partialA_a1<-rep(0,ni)
				partialB_a2<-rep(0,ni)
				partialB_a1<-rep(0,ni)
				partialA_c2<-rep(0,ni)
				partialA_c1<-rep(0,ni)
				partialB_c2<-rep(0,ni)
				partialB_c1<-rep(0,ni)
			}
			if (method=="mean-mean") {
				partialA_b2<-rep(0,ni)
				partialA_b1<-rep(0,ni)
				partialB_b2<-rep(1/ni,ni)
				partialB_b1<-rep(-A/ni,ni)
				partialA_a2<-rep(-sum(a1)/(sum(a2)^2),ni)
				partialA_a1<-rep(1/sum(a2),ni)
				partialB_a2<--partialA_a2*mean(b1)
				partialB_a1<--partialA_a1*mean(b1)
				partialA_c2<-rep(0,ni)
				partialA_c1<-rep(0,ni)
				partialB_c2<-rep(0,ni)
				partialB_c1<-rep(0,ni)
			}
			if (method=="mean-gmean") {
				partialA_b2<-rep(0,ni)
				partialA_b1<-rep(0,ni)
				partialB_b2<-rep(1/ni,ni)
				partialB_b1<-rep(-A/ni,ni)
				partialA_a2<--1/ni*A/a2
				partialA_a1<-1/ni*A/a1
				partialB_a2<--partialA_a2*mean(b1)
				partialB_a1<--partialA_a1*mean(b1)
				partialA_c2<-rep(0,ni)
				partialA_c1<-rep(0,ni)
				partialB_c2<-rep(0,ni)
				partialB_c1<-rep(0,ni)
			}
			if (method=="Haebara") {
				P1<-t(P1)
				P2<-t(P2)
				tmp<-((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2)*(P1-c1)/(1-c1)^2*(1-P1)*a1^2*D
				tmp<-apply(tmp,2,sum)*wt
				partialSIR_AB<-matrix(0,2,2)
				for (i in 1:length(ab)) 
					partialSIR_AB<- partialSIR_AB-tmp[i]*c(ab[i],1)%*%t(c((ab[i]-B)/A^2,1/A))
				if (itmp==1) {partialSIR_AB[1,]<-0
					partialSIR_AB[,1]<-0}
				abMAT<-matx(ab,ni)
				wtMAT<-matx(wt,ni)
				tmp<-(P1-c1)/(1-c1)*(1-P1)*a1*(1-P2)/(1-c2)*wtMAT
				tmp_a2<-tmp*(P2-c2)*D*(abMAT-b2)
				partialSIR_a2<-rbind(rowSums(tmp_a2*abMAT),rowSums(tmp_a2))
				tmp_b2<-tmp*(P2-c2)*D*(-a2)
				partialSIR_b2<-rbind(rowSums(tmp_b2*abMAT),rowSums(tmp_b2))
				partialSIR_c2<-rbind(rowSums(tmp*abMAT),rowSums(tmp))
				tmp_a1<-(((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2)*a1/(1-c1)*
					D*(abMAT-A*b1-B)/A+P2-P1)*(P1-c1)/(1-c1)*(1-P1)*wtMAT
				partialSIR_a1<-rbind(rowSums(tmp_a1*abMAT),rowSums(tmp_a1))
				tmp_b1<-((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2)*
					(P1-c1)/(1-c1)^2*(1-P1)*(-D*a1^2)*wtMAT
				partialSIR_b1<-rbind(rowSums(tmp_b1*abMAT),rowSums(tmp_b1))
				tmp_c1<-((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2-(P2-P1)*(1-P1))*
					(1-P1)*a1/(1-c1)^2*wtMAT
				partialSIR_c1<-rbind(rowSums(tmp_c1*abMAT),rowSums(tmp_c1))
			}
			if (method=="Stocking-Lord") {
				P1<-t(P1)
				P2<-t(P2)
				tmp1<--(P1-c1)/(1-c1)*(1-P1)*a1
				tmp2<-P2-P1
				tmp3<-(1+c1-2*P1)/(1-c1)*a1
				tmp4<-(P1-c1)/(1-c1)*(1-P1)*a1*D
				tmp<-colSums(tmp1)*colSums(tmp4)+colSums(tmp2)*colSums(tmp3*tmp4)
				tmp<-tmp*wt
				partialSIR_AB<-matrix(0,2,2)
				for (i in 1:length(ab)) 
					partialSIR_AB<- partialSIR_AB-tmp[i]*c(ab[i],1)%*%t(c((ab[i]-B)/A^2,1/A))
				if (itmp==1) {partialSIR_AB[1,]<-0
					partialSIR_AB[,1]<-0}
				abMAT<-matrix(rep(ab,each=ni),nrow=ni)
				abMAT<-matx(ab,ni)
				wtMAT<-matx(wt,ni)
				tmp<-(P1-c1)/(1-c1)*(1-P1)*a1
				tmp<-matx(colSums(tmp),ni)*(1-P2)/(1-c2)*wtMAT
				tmp_a2<-tmp*(P2-c2)*D*(abMAT-b2)
				partialSIR_a2<-rbind(rowSums(tmp_a2*abMAT),rowSums(tmp_a2))
				tmp_b2<-tmp*(P2-c2)*D*(-a2)
				partialSIR_b2<-rbind(rowSums(tmp_b2*abMAT),rowSums(tmp_b2))
				partialSIR_c2<-rbind(rowSums(tmp*abMAT),rowSums(tmp))
				tmp1<--(P1-c1)/(1-c1)*(1-P1)*a1
				tmp2<-P2-P1
				tmp3<-(1+c1-2*P1)/(1-c1)*a1
				tmp4<-D*(abMAT-A*b1-B)/A
				tmp6<-(P1-c1)/(1-c1)*(1-P1)
				tmp_a1<-matx(colSums(tmp1),ni)*(tmp4*tmp6)+matx(colSums(tmp2),ni)*(tmp3*tmp4*tmp6)+matx(colSums(tmp2),ni)*(tmp6)
				tmp_a1<-tmp_a1*wtMAT
				partialSIR_a1<-rbind(rowSums(tmp_a1*abMAT),rowSums(tmp_a1))
				tmp61<- -tmp6*D*a1
				tmp_b1<-matx(colSums(tmp1),ni)*tmp61+matx(colSums(tmp2),ni)*(tmp3*tmp61)
				tmp_b1<-tmp_b1*wtMAT
				partialSIR_b1<-rbind(rowSums(tmp_b1*abMAT),rowSums(tmp_b1))
				tmp7<-(P1-c1)/(1-c1)*a1
				tmp8<-(1-P1)/(1-c1)
				tmp_c1<-matx(colSums(tmp1),ni)*tmp8-matx(colSums(tmp2),ni)*(tmp7*tmp8)
				tmp_c1<-tmp_c1*wtMAT
				partialSIR_c1<-rbind(rowSums(tmp_c1*abMAT),rowSums(tmp_c1))
			}

			if (method=="mean-mean" | method=="mean-sigma" | method=="mean-gmean") { 
				if (itmp==1) mat<-cbind(
					c(partialA_b1,partialA_b2),
					c(partialB_b1,partialB_b2))
				if (itmp==2) mat<-cbind(
					c(partialA_b1,partialA_a1,partialA_b2,partialA_a2),
					c(partialB_b1,partialB_a1,partialB_b2,partialB_a2))
				if (itmp==3) mat<-cbind(
					c(partialA_b1,partialA_a1,partialA_c1,partialA_b2,partialA_a2,partialA_c2),
					c(partialB_b1,partialB_a1,partialB_c1,partialB_b2,partialB_a2,partialB_c2))
			}
			if (method=="Haebara" | method=="Stocking-Lord") {
				if (itmp==1) partialSIR_gamma<-cbind(partialSIR_b1,partialSIR_b2)
				if (itmp==1) partialSIR_gamma[1,]<-0
				if (itmp==2) partialSIR_gamma<-cbind(partialSIR_b1,partialSIR_a1,partialSIR_b2,partialSIR_a2)
				if (itmp==3) partialSIR_gamma<-cbind(partialSIR_b1,partialSIR_a1,partialSIR_c1,partialSIR_b2,partialSIR_a2,partialSIR_c2)
				if (itmp==1) {invpartialSIR_AB<-partialSIR_AB
					invpartialSIR_AB[2,2]<-1/invpartialSIR_AB[2,2]}
				if (itmp==1) mat1<--invpartialSIR_AB%*%partialSIR_gamma
				if (itmp>1) mat1<--solve(partialSIR_AB)%*%partialSIR_gamma
				mat<-t(mat1)
			}
			if (ni==0) mat<-matrix(NA,2,2)
			colnames(mat)<-c("A","B")
			rownames(mat)<-c(paste(comuni,suff1,sep=""),paste(comuni,suff2,sep=""))
			if (ni>0) {
				var1<-var1[tab$Row.names,tab$Row.names]
				var2<-var2[tab$Row.names,tab$Row.names]
				rownames(var1)<-colnames(var1)<-paste(rownames(var1),suff1,sep="")
				rownames(var2)<-colnames(var2)<-paste(rownames(var2),suff2,sep="")
				var12<-blockdiag(var1,var2)
				varAB<-t(mat)%*%var12%*%mat
			}
			else { 
				varAB<-matrix(NA,2,2)
				var12<-NULL 
			}
		}
		if (is.null(var1) | is.null(var2)) {
			var12<-NULL
			mat<-NULL
			varAB<-matrix(NA,2,2)
		}
		if (!is.null(var12)) {
			var1<-mod1$var
			var2<-mod2$var
			rownames(var1)<-colnames(var1)<-paste(rownames(var1),suff1,sep="")
			rownames(var2)<-colnames(var2)<-paste(rownames(var2),suff2,sep="")
			varFull<-list(var1,var2)
			#names(varFull)<-c(name1,name2)
		}
		else varFull<-NULL
		taball$value12<-NA
		taball$value12[1:niall]<-A*taball$value1[1:niall]+B
		if (itmp>1) taball$value12[(niall+1):(2*niall)]<-taball$value1[(niall+1):(2*niall)]/A
		if (itmp==3) taball$value12[(2*niall+1):(3*niall)]<-taball$value1[(2*niall+1):(3*niall)]
		colnames(taball)<-c("Item",name1,name2,paste(name1,name2,sep=".as."))
		out<-list(tab1=tab1,tab2=tab2,tab=taball,var12=var12,varFull=varFull,
		partial=mat,A=A,B=B,varAB=varAB,commonitem=list(comuni),suffixes=c(suff1,suff2),ni=ni)
	}
	out$forms<-forms
	out$method<-method
	out$itmp<-itmp
	class(out) <- "eqc"
	return(out)
}


matx<-function(vect,n) {
	rep(1,n)%x%t(vect)
}


print.eqc<-function(x, ...)
{
	cat("Direct equating coefficients \n")
	cat("Method: ")
	cat(x$method,"\n")
	cat("Link: ")
	cat(x$forms,"\n")
}


summary.eqc <- function(object, ...)
{
	if (sum(object$ni>0)==length(object$ni)) {
		ct<-cbind(Estimate=c(object$A,object$B),StdErr=c(sqrt(diag(object$varAB))))
		rownames(ct)<-c("A","B")
	}
	else ct<-NULL
	out<-list(forms=object$forms,method=object$method,coefficients=ct)
	class(out)<-"summary.eqc"
	return(out)
}

print.summary.eqc <- function(x, ...)
{
	cat("Link: ")
	cat(x$forms,"\n")
	cat("Method: ")
	cat(x$method,"\n")
	if (!is.null(x$coefficients)) { cat("Equating coefficients:\n")
		print(x$coefficients,digits=5)}
	else cat("no common items\n")
}


alldirec<-function(mods,method="mean-mean",all=FALSE,quadrature=TRUE,nq=30,direction="both")
{
	options(warn=-1)
	nt<-length(mods)
	direclist<-list()
	k<-1
	if (direction=="both") {
		for (i in 1:nt) {
			for (j in 1:nt) {
				if (i!=j) {
					tmp<-direc(mods[i],mods[j],suff1=paste(".",i,sep=""),suff2=paste(".",j,sep=""),method=method,quadrature=quadrature,nq=nq)
					if (tmp$ni>0 | all) {
						direclist[[k]]<-tmp
						names(direclist)[[k]]<-tmp$forms
						k<-k+1
					}
				}
			}
		}
	}
	if (direction=="back") {
		for (i in 2:nt) {
			for (j in 1:(i-1)) {
				tmp<-direc(mods[i],mods[j],suff1=paste(".",i,sep=""),suff2=paste(".",j,sep=""),method=method,quadrature=quadrature,nq=nq)
				if (tmp$ni>0 | all) {
					direclist[[k]]<-tmp
					names(direclist)[[k]]<-tmp$forms
					k<-k+1
				}
			}
		}
	}
	if (direction=="forward") {
		for (j in 2:nt) {
			for (i in 1:(j-1)) {
				tmp<-direc(mods[i],mods[j],suff1=paste(".",i,sep=""),suff2=paste(".",j,sep=""),method=method,quadrature=quadrature,nq=nq)
				if (tmp$ni>0 | all) {
					direclist[[k]]<-tmp
					names(direclist)[[k]]<-tmp$forms
					k<-k+1
				}
			}
		}
	}
	class(direclist)<-"eqclist"
	options(warn=0)
	return(direclist)
}


print.eqclist<-function(x, ...)
{
	cat("Direct equating coefficients \n")
	cat("Method: ")
	cat(x[[1]]$method,"\n")
	cat("Links: \n")
	for (i in 1:length(x)) cat(x[[i]]$forms,"\n")
}

summary.eqclist <- function(object, link=NULL, ...)
{
	if(is.null(link)) link<-sapply(object,FUN=function(x) x$forms)
	out<-list()
	j<-1
	for (i in 1:length(object))
		if (object[[i]]$forms%in%link) {
			out[[j]]<-summary(object[[i]])
			j<-j+1
		}
	class(out)<-"summary.eqclist"
	return(out)
}

print.summary.eqclist<-function(x, ...)
{
	for (i in 1:length(x)) {
		print(x[[i]])
		cat("\n\n")
	}
}




chainec<-function(r=NULL,direclist,f1=NULL,f2=NULL,pths=NULL)
{
	if (is.null(r) & is.null(pths)) stop("argument \"r\" needs to be specified if argument \"pths\" is NULL.")
	if (is.null(r) & !is.null(pths)) r<-ncol(pths)
	if (r<3) stop("r should be at least 3.")
	if (is.null(pths)) {
		sel<-sapply(direclist,FUN= function(x)(x$ni!=0))
		nl<-names(direclist)[sel]
		nll<-strsplit(nl,split=".",fixed=TRUE)
		l<-data.frame(f1=sapply(nll,FUN=function(x) x[1]),f2=sapply(nll,FUN=function(x) x[2]),stringsAsFactors=FALSE)
		if (is.null(f1)) pths<-l
		if (!is.null(f1)) pths<-l[l$f1==f1,]
		colnames(pths)<-paste(colnames(pths),1,sep=".")
		if (r>3) {
			for (k in 1:(r-3)) {
				pths<-merge(pths,l,by.x=k+1,by.y=1)
				colnames(pths)<-paste(colnames(pths),k+1,sep=".")
				pths<-pths[,c(2:(k+1),1,k+2)]
				pths<-pths[pths[,k]!=pths[,k+2],]
			}
		}
		if (is.null(f2)) pths<-merge(pths,l,by.x=r-1,by.y=1)
		if (!is.null(f2)) pths<-merge(pths,l[l$f2==f2,],by.x=r-1,by.y=1)
		pths<-pths[,c(2:(r-1),1,r)]
		pths<-pths[pths[,r-2]!=pths[,r],] 
	}
   
	pths<-pths[pths[,1]!=pths[,r],]
	if (nrow(pths)==0) stop("There are not paths of length ", r, ".")
	
	nomi<-pths[,1]
	for (k in 2:r) nomi<-paste(nomi,pths[,k],sep=".")
	out<- vector("list", nrow(pths))
	for (j in 1:nrow(pths)) {
		name1<-as.character(pths[j,1])
		name2<-as.character(pths[j,r])
		ni<-c()
		A<-1
		B<-0
		partialA<-c()
		partialB<-c()
		varAll<-matrix(0,0,0)
		varFull<-list()
		suffixes<-c()
		comuni<-list()
		missing<-FALSE
		varNULL<-FALSE
		for (k in 1:(r-1)) {
			nome<-paste(pths[j,k],pths[j,k+1],sep=".")
			link<-direclist[[nome]]
			if (k==1)   tab1<-link$tab1
			if (k==r-1) tab2<-link$tab2
			if (!is.null(link)) {
				ni<-c(ni,link$ni)
				if (link$ni!=0) {
					partialAk<-A*link$partial[,1] 
					partialA<-partialA*link$A
					partialA<-c(partialA,partialAk)

					partialBk<-B*link$partial[,1]+link$partial[,2] 
					partialB<-partialB*link$A
					partialB<-c(partialB,partialBk)

					A<-link$A*A  #A1...k=Ak-1k*A1...k-1
					B<-link$B+link$A*B #B1...k=Bkk-1+Akk-1*B1...k-1
					
					if (!is.null(link$var12)) varAll<-blockdiag(varAll,link$var12)
					if (is.null(link$var12)) varNULL<-TRUE
					if (!is.null(link$varFull) & k==1) varFull[[k]]<-(link$varFull[[1]])
					if (!is.null(link$varFull)) varFull[[k+1]]<-(link$varFull[[2]])
					if (k==1) suffixes<-c(suffixes,link$suffixes[1])
					suffixes<-c(suffixes,link$suffixes[2])
					comuni[[k]]<-link$commonitem[[1]]
				}
				else {
					warning("forms ",nome," have no common items\n")
					missing<-TRUE
				}
			}
			if (is.null(link)) {
				warning("link of forms ",nome," is missing\n")
				missing<-TRUE
				ni<-c(ni,0)
			}
		}

		if (!missing) {
			partialA<-tapply(partialA,names(partialA),sum)
			partialB<-tapply(partialB,names(partialB),sum)

			mat<-merge(partialA,partialB,by=0)
			nom<-mat$Row.names
			mat<-mat[,-1]
			colnames(mat)<-c("A","B")
			rownames(mat)<-nom
			mat<-as.matrix(mat)
			sel<-unique(rownames(varAll))
			varAB<-t(mat[sel,])%*%varAll[sel,sel]%*%mat[sel,]
			if(varNULL) varAB<-matrix(NA,2,2)
			taball<-merge(tab1,tab2,by=0,suffixes=c(.1,.2),all=T)
			niall<-nrow(taball)/link$itmp
			taball$value12<-NA
			taball$value12[1:niall]<-A*taball$value1[1:niall]+B
			if (link$itmp>1) taball$value12[(niall+1):(2*niall)]<-taball$value1[(niall+1):(2*niall)]/A
			if (link$itmp==3) taball$value12[(2*niall+1):(3*niall)]<-taball$value1[(2*niall+1):(3*niall)]
			colnames(taball)<-c("Item",name1,name2,paste(name1,name2,sep=".as."))
			out[[j]]$tab1<-tab1
			out[[j]]$tab2<-tab2
			out[[j]]$tab<-taball
			out[[j]]$varAll<-varAll
			out[[j]]$varFull<-varFull
			out[[j]]$partial<-mat[sel,]
			out[[j]]$A<-A
			out[[j]]$B<-B
			out[[j]]$varAB<-varAB
			out[[j]]$commonitem<-comuni
			out[[j]]$suffixes<-suffixes
		}
		out[[j]]$ni<-ni
		out[[j]]$forms<-nomi[j]
		if (!is.null(link)) out[[j]]$method<-link$method
		if (!is.null(link)) out[[j]]$itmp<-link$itmp
		else out[[j]]$method<-""
		class(out[[j]])<-"ceqc"
	} 
	names(out)<-nomi
	class(out) <- "ceqclist"
	return(out)
}


print.ceqc<-function(x, ...)
{
	cat("Chain equating coefficients \n")
	cat("Method: ")
	cat(x$method,"\n")
	cat("Path: ")
	cat(x$forms,"\n")
}


summary.ceqc <- function(object, ...)
{
	if (sum(object$ni>0)==length(object$ni)) {
		ct<-cbind(Estimate=c(object$A,object$B),StdErr=c(sqrt(diag(object$varAB))))
		rownames(ct)<-c("A","B")
	}
	else ct<-NULL
	out<-list(forms=object$forms,method=object$method,coefficients=ct)
	class(out)<-"summary.ceqc"
	return(out)
}

print.summary.ceqc <- function(x, ...)
{
	cat("Path: ")
	cat(x$forms,"\n")
	cat("Method: ")
	cat(x$method,"\n")
	if (!is.null(x$coefficients)) { cat("Equating coefficients:\n")
		print(x$coefficients,digits=5)}
	else cat("no common items\n")
}

print.ceqclist<-function(x, ...)
{
	cat("Chain equating coefficients \n")
	cat("Method: ")
	cat(x[[1]]$method,"\n")
	cat("Paths: \n")
	for (i in 1:length(x)) cat(x[[i]]$forms,"\n")
}

summary.ceqclist <- function(object, path=NULL, ...)
{
	if(is.null(path)) path<-sapply(object,FUN=function(x) x$forms)
	out<-list()
	j<-1
	for (i in 1:length(object))
		if (object[[i]]$forms%in%path) {
			out[[j]]<-summary(object[[i]])
			j<-j+1
		}
	class(out)<-"summary.ceqclist"
	return(out)
}


print.summary.ceqclist<-function(x, ...)
{
	for (i in 1:length(x)) {
		print(x[[i]])
		cat("\n\n")
	}
}


bisectorec<-function(ecall,mods=NULL,weighted=TRUE,unweighted=TRUE)
{
	if (!is.null(mods)) cat("Note: from version 2.0 argument mods can be left unspecified.")
	if (length(table(sapply(ecall,FUN=function(x) x$method)))!=1) stop("ecall contains different methods.")
	itmp<-sapply(ecall,FUN=function(x) x$itmp)
	if (length(table(itmp))>1) stop("Mixed models not allowed. Number of item parameters differs.")
	varNULL<-FALSE
	if (any(sapply(ecall,FUN=function(x) is.na(x$varAB)))) varNULL<-TRUE
	if (varNULL & weighted) {
		stop("Weighted bisector is unfeasible with NULL covariance matrix")
		weighted<-FALSE
	}
	if (!varNULL) {
		part<-lapply(ecall,FUN=function(x) data.frame(A=x$partial[,1],
					B=x$partial[,2],stringsAsFactors = FALSE))
		for (i in 1:length(part)) {
			part[[i]]$path<-names(part)[i]
			part[[i]]$par<-rownames(part[[i]])
		}
		partall<-part[[1]]
		for (i in 2:length(part))  partall<-rbind(partall,part[[i]])
		partall$link<-path2link(partall$path)
	}
	else partall<-NULL
	coall<-data.frame(t(sapply(ecall,FUN=function(x) x[c("A","B")])))
	coall$seA<-sapply(ecall,FUN=function(x) x$varAB[1,1]^0.5)
	coall$seB<-sapply(ecall,FUN=function(x) x$varAB[2,2]^0.5)
	for (i in 1:4) coall[,i]<-unlist(coall[,i])
	coall$path<-rownames(coall)
	coall$link<-path2link(coall$path)

	coall$weights<-NA
	links<-sort(unique(coall$link))
	#if (!varNULL) varFull<-VarExt(mods)
	if (!varNULL) varFull<-lapply(ecall,FUN=function(x) x$varFull)
	else varFull<-NULL
	suffixes<-lapply(ecall,FUN=function(x) x$suffixes)
	if (unweighted) {
		coall$weights<-1
		bisl<-bisco(coall,varFull,partall)
		bis<-data.frame(link=sapply(bisl,FUN=function(x) x$link),path="bisector",
			A=sapply(bisl,FUN=function(x) x$A),B=sapply(bisl,FUN=function(x) x$B),
			seA=sapply(bisl,FUN=function(x) x$varAB[1,1]^0.5),
			seB=sapply(bisl,FUN=function(x) x$varAB[2,2]^0.5),weights=NA)
		for (i in links) {
			tab1<-ecall[coall$link==i][[1]]$tab1
			tab2<-ecall[coall$link==i][[1]]$tab2
			taball<-merge(tab1,tab2,by=0,all=T)
			taball$value12<-NA
			colnames(taball)<-colnames(ecall[coall$link==i][[1]]$tab)
			Item<-taball$Item
			A<-bisl[[i]]$A
			B<-bisl[[i]]$B
			taball[,4][substr(Item,1,6)=="Dffclt"]<-A*taball[,2][substr(Item,1,6)=="Dffclt"]+B
			taball[,4][substr(Item,1,6)=="Dscrmn"]<-taball[,2][substr(Item,1,6)=="Dscrmn"]/A
			taball[,4][substr(Item,1,6)=="Gussng"]<-taball[,2][substr(Item,1,6)=="Gussng"]
			bisl[[i]]$tab<-taball
			bisl[[i]]$itmp<-itmp[1]
			suff<-suffixes[i==coall$link][[1]]
			bisl[[i]]$suffixes<-suff[c(1,length(suff))]
		}
	}
	else bisl<-NULL
	if (weighted) {
		for (i in links) {
			colink<-coall[coall$link==i,]
			nl<-nrow(colink)
			weights<-rep(1,nl)
			partlink<-partall[partall$link==i,]
			linkvf<-path2link(names(varFull))
			varFull1<-varFull[linkvf==i]
			o<-optim(par=weights,fn=VarTrasf,colink=colink,varFull=varFull1,partlink=partlink,control=list(maxit=10000,reltol=1e-5))
			coall[coall$link==i,]$weights<-abs(o$par)
		}
		wbisl<-bisco(coall,varFull,partall)
		#wbis$path<-"weighted bisector"
		wbis<-data.frame(link=sapply(wbisl,FUN=function(x) x$link),path="weighted bisector",
			A=sapply(wbisl,FUN=function(x) x$A),B=sapply(wbisl,FUN=function(x) x$B),
			seA=sapply(wbisl,FUN=function(x) x$varAB[1,1]^0.5),
			seB=sapply(wbisl,FUN=function(x) x$varAB[2,2]^0.5),weights=NA)
		for (i in links) {
			tab1<-ecall[coall$link==i][[1]]$tab1
			tab2<-ecall[coall$link==i][[1]]$tab2
			taball<-merge(tab1,tab2,by=0,all=T)
			taball$value12<-NA
			colnames(taball)<-colnames(ecall[coall$link==i][[1]]$tab)
			Item<-taball$Item
			A<-wbisl[[i]]$A
			B<-wbisl[[i]]$B
			taball[,4][substr(Item,1,6)=="Dffclt"]<-A*taball[,2][substr(Item,1,6)=="Dffclt"]+B
			taball[,4][substr(Item,1,6)=="Dscrmn"]<-taball[,2][substr(Item,1,6)=="Dscrmn"]/A
			taball[,4][substr(Item,1,6)=="Gussng"]<-taball[,2][substr(Item,1,6)=="Gussng"]
			wbisl[[i]]$tab<-taball
			wbisl[[i]]$itmp<-itmp[1]
			suff<-suffixes[i==coall$link][[1]]
			wbisl[[i]]$suffixes<-suff[c(1,length(suff))]
		}
	}
	else wbisl<-NULL
	sel<-c("link","path","A","B","seA","seB","weights")
	if (unweighted) coall<-rbind(coall[,sel],bis[,sel])
	if (weighted) coall<-rbind(coall[,sel],wbis[,sel])
	coall<-coall[order(coall[,1]),]
	rownames(coall)<-NULL
	meq<-list(coef=coall,method=ecall[[1]]$method,bis=bisl,wbis=wbisl)
	class(meq)<-"meqc"
	return(meq)
}


bisco<-function(coall,varFull,partall)
{
	coall$w<-BisW(coall$A)*coall$weights
	W<-tapply(coall$w,coall$link,FUN=sum)
	mA<-tapply(coall$A*coall$w,coall$link,sum)/W
	mB<-tapply(coall$B*coall$w,coall$link,sum)/W
	coall$W<-W[coall$link]
	coall$mA<-mA[coall$link]
	coall$mB<-mB[coall$link]
	#out<-data.frame(link=names(mA),A=mA,B=mB,seA=NA,seB=NA,corAB=NA,path="bisector",weights=NA)
	out1<-list()
	if (!is.null(varFull)) {
		coall$partAA<-((1-coall$A^2/(1+coall$A^2))*coall$w*coall$W+
						coall$mA*coall$W*coall$weights*coall$A*(1+coall$A^2)^(-1.5))/
						coall$W^2
		coall$partBA<- -(coall$A*coall$B/(1+coall$A^2)*coall$w*coall$W+
						coall$mB*coall$W*coall$weights*coall$A*(1+coall$A^2)^(-1.5))/
						coall$W^2
		coall$partBB<- coall$w/coall$W

		links<-unique(coall$link)
		for (i in 1:length(links)){
			coi<-coall[coall$link==links[i],]
			parti<-partall[partall$link==links[i],]
			parti$partAA<-coi[parti$path,]$partAA
			parti$partBA<-coi[parti$path,]$partBA
			parti$partBB<-coi[parti$path,]$partBB

			partialAmean<-tapply(parti$partAA*parti$A,parti$par,FUN=sum)
			partialBmean<-tapply(parti$partBA*parti$A+parti$partBB*parti$B,parti$par,FUN=sum)
			partialmean<-cbind(partialAmean,partialBmean)
			colnames(partialmean)<-c("A","B")
			sel<-rownames(partialmean)
			#strsplit(sel,".",fixed=TRUE)
			linkvf<-path2link(names(varFull))
			varFull1<-varFull[linkvf==links[i]]
			varFull2<-list()
			varAll<-matrix(0,0,0)
			ptt<-c()
			for (k in 1:length(varFull1)) {
				for (j in 1:length(varFull1[[k]])) {
					v1<-varFull1[[k]][[j]]
					varAll<-blockdiag(varAll,v1)
					if (any(!rownames(v1)%in%ptt)) varFull2<-c(varFull2,varFull1[[k]][j])
					ptt<-c(ptt,rownames(v1))
			}}
			varAB<-t(partialmean)%*%varAll[sel,sel]%*%partialmean
			#seA<-varAB[1,1]^0.5
			#seB<-varAB[2,2]^0.5
			#out[i,]$seA<-seA
			#out[i,]$seB<-seB
			#out[i,]$corAB<-varAB[1,2]/(seA*seB)
			out1[[i]]<-list(link=links[i],partial=partialmean,A=mA[i],B=mB[i],varAB=varAB,varFull=varFull2)
			names(out1)[i]<-links[i]
		}
	}
	return(out1)
}


VarExt<-function(mods)
{
	varFull<-list()
	for (i in 1:length(mods)) {
		var1<-mods[[i]]$var
		rownames(var1)<-colnames(var1)<-paste(rownames(var1),i,sep=".")
		varFull[[i]]<-var1
	}
	return(varFull)
}


VarTrasf<-function(weights,colink,varFull,partlink)
{
	colink$weights<-abs(weights)
	colink$w<-BisW(colink$A)*colink$weights
	W<-sum(colink$w)
	mA<-sum(colink$A*colink$w)/W
	mB<-sum(colink$B*colink$w)/W
	colink$W<-W
	colink$mA<-mA
	colink$mB<-mB

	colink$partAA<-((1-colink$A^2/(1+colink$A^2))*colink$w*colink$W+
					colink$mA*colink$W*colink$weights*colink$A*(1+colink$A^2)^(-1.5))/
					colink$W^2
	colink$partBA<- -(colink$A*colink$B/(1+colink$A^2)*colink$w*colink$W+
					colink$mB*colink$W*colink$weights*colink$A*(1+colink$A^2)^(-1.5))/
					colink$W^2
	colink$partBB<- colink$w/colink$W

	partlink$partAA<-colink[partlink$path,]$partAA
	partlink$partBA<-colink[partlink$path,]$partBA
	partlink$partBB<-colink[partlink$path,]$partBB

	partialAmean<-tapply(partlink$partAA*partlink$A,partlink$par,FUN=sum)
	partialBmean<-tapply(partlink$partBA*partlink$A+partlink$partBB*partlink$B,partlink$par,FUN=sum)
	partialmean<-cbind(partialAmean,partialBmean)
	sel<-rownames(partialmean)

	varAll<-matrix(0,0,0)
	for (i in 1:length(varFull)) {
		for (j in 1:length(varFull[[i]])) {
			varAll<-blockdiag(varAll,varFull[[i]][[j]])
	}}

	varAB<-t(partialmean)%*%varAll[sel,sel]%*%partialmean
	out<-sum(diag(varAB))
	return(out)
}



BisW<-function(x) (1+x^2)^(-0.5)


print.meqc<-function(x, ...)
{
	bis<-any(x$coef$path=="bisector")
	wbis<-any(x$coef$path=="weighted bisector")
	if (bis & !wbis) cat("Bisector equating coefficients \n")
	if (!bis & wbis) cat("Weighted bisector equating coefficients \n")
	if (bis & wbis) cat("Bisector and weighted bisector equating coefficients \n")
	cat("Method: ")
	cat(x$method,"\n")
	links<-unique(x$coef$link)
	for (i in links) {
		cat("\n")
		cat("Link: ")
		cat(i, "\n")
		cat("  Paths: \n")
		paths<-x$coef$path[x$coef$link==i]
		paths<-paths[paths!="bisector" & paths!="weighted bisector"]
		for (j in paths) cat(" ",j, "\n")
	}
}

summary.meqc <- function(object, ...)
{
	method<-object$method
	object<-object$coef
	link<-sort(unique(object$link))
	tab<-list()
	for (i in 1:length(link)) {
		objecti<-object[object$link==link[i],]
		objecti$coefA<-"A"
		objecti$coefB<-"B"
		objecti1<-objecti[,c("coefA","path","A","seA")]
		objecti2<-objecti[,c("coefB","path","B","seB")]
		colnames(objecti1)<-colnames(objecti2)<-c("","Path","Estimate","StdErr")
		tab[[i]]<-rbind(objecti1,objecti2)
	}
	out<-list(link=link,method=method,coefficients=tab)
	class(out)<-"summary.meqc"
	return(out)
}

print.summary.meqc <- function(x, ...)
{
	for (i in 1:length(x$link)) {
		cat("Link: ")
		cat(x$link[i],"\n")
		cat("Method: ")
		cat(x$method,"\n")
		cat("Equating coefficients:\n")
		print(x$coefficients[[i]],row.names=F,digits=5)
		cat("\n")
	}
}


path2link<-function(x)
{
	tt<-strsplit(x,".",fixed=TRUE)
	tmp1<-sapply(tt,FUN=function(x) x[1])
	tmp2<-sapply(tt,FUN=function(x) x[length(x)])
	return(paste(tmp1,tmp2,sep="."))
}



blockdiag<-function(m1,m2)
{
	r1<-nrow(m1)
	r2<-nrow(m2)
	c1<-ncol(m1)
	c2<-ncol(m2)
	out<-matrix(0,r1+r2,c1+c2)
	if (r1>0) out[1:r1,1:c1]<-m1
	out[(r1+1):(r1+r2),(c1+1):(c1+c2)]<-m2
	rownames(out)<-c(rownames(m1),rownames(m2))
	colnames(out)<-c(colnames(m1),colnames(m2))
	return(out)
} 


convert<-function(A,B,coef=NULL,person.par=NULL)
{
	if (!is.null(coef)) {
		itms<-names(coef)
		Dffclt<-coef[substr(itms,1,6)=="Dffclt"]
		Dscrmn<-coef[substr(itms,1,6)=="Dscrmn"]
		Gussng<-coef[substr(itms,1,6)=="Gussng"]
		Dffclt<-Dffclt*A+B
		if (length(Dscrmn)>0) Dscrmn<-Dscrmn/A
		coef1<-c(Dffclt,Dscrmn,Gussng)
	}
	else coef1<-NULL
	if (!is.null(person.par)) person.par1<-person.par*A+B
	else person.par1<-NULL
	return(list(coef=coef1,person.par=person.par1))
}


import.flexmirt<-function(fnamep,fnamev=NULL,fnameirt=NULL,display=TRUE,digits=2) {
	par<-read.table(fnamep,fill=TRUE)
	ngr<-sum(par$V1==0,na.rm=TRUE)
	if (ngr>1) stop("Cannot handle multiple groups.")
	if(max(par$V4,na.rm=TRUE)>1) stop("Cannot handle multiple factors.")
	if(max(par$V6,na.rm=TRUE)>2) stop("Cannot handle multiple response models.")
	par<-par[par$V1!=0 & !is.na(par$V1),]
	if(length(unique(par$V5))>1) stop("Cannot handle mixed item types.")
	if(unique(par$V5)==3) stop("Cannot handle nominal categories models.")
	if(unique(par$V5)==2) itmp=2
	if(unique(par$V5)==1) itmp=3
	if(length(unique(par$V8))==1) itmp=1
	ifelse(all(par$V8==1),Rasch<-TRUE,Rasch<-FALSE)
	if (itmp==1 | itmp==2) {p<-par[,c(7,8)]
		colnames(p)<-c("c","a")}
	if (itmp==3) {p<-par[,7:9]
		colnames(p)<-c("logit-g","c","a")}
	rownames(p)<-par$V2
	p<-as.matrix(p)

	if (!is.null(fnamev)) {
		if (is.null(fnameirt)) stop("fnameirt required if fnamev is not NULL.")
		vm<-read.table(fnamev,sep=",")
		vm<-as.matrix(vm[,colSums(is.na(vm))!=nrow(vm)])
		colnames(vm)<-NULL
		n<-nrow(p)
		if (itmp==1 & !Rasch) reord<-c((n+1):(2*n),1:n)
		if (itmp==1 & Rasch) reord<-1:n
		if (itmp==2) reord<-c((n+1):(2*n),1:n)
		if (itmp==3) reord<-c((2*n+1):(3*n),(n+1):(2*n),1:n)
		
		irt <- readLines(fnameirt)
		startRow <- grep(" *Item *Label *P# ", irt)
		nm <- unlist(strsplit(irt[startRow], " +"))[-1]
		widths <- strsplit(irt[startRow], "(?<=[A-Za-z#.])(?=\\s)", perl = TRUE)[[1]]
		widths <- sapply(widths, nchar)
		tab <- read.fwf(fnameirt, widths, skip = startRow, na.strings = "----",
			col.names = nm, comment.char = "", sep = ",", 
			nrows = as.numeric(n), strip.white = TRUE, check.names = FALSE)
		param.number <- unlist(tab[, names(tab) == "P#"])
		pn.max<-max(param.number,na.rm=TRUE)
		num.na<-sum(is.na(param.number))
		if (num.na>0 & !Rasch) {
			fixed<-(pn.max+1):(pn.max+num.na)
			param.number[is.na(param.number)]<-fixed
			nf<-length(fixed)
			nvm1<-nrow(vm)+nf
			vm1<-matrix(0,nvm1,nvm1)
			vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
			vm<-vm1
		}
		if (num.na>0 & Rasch) param.number<-param.number[!is.na(param.number)]
		vm<-vm[param.number,param.number]
		if (!Rasch & n*ncol(p)!=nrow(vm)) stop("Number of parameters and dimension of the covariance matrix do not match.")
		if (Rasch & n!=nrow(vm)) stop("Number of parameters and dimension of the covariance matrix do not match.")
		vm<-vm[reord,reord]
	}
	else vm<-NULL
	if (display) {
		out<-matrix(NA,nrow(p),ncol(p)*3)
		out[,seq(2,ncol(out),by=3)]<-p
		se<-diag(vm)^0.5
		if (!is.null(vm)) out[,seq(3,ncol(out),by=3)]<-se
		if (!is.null(vm)) out[,seq(1,ncol(out),by=3)]<-param.number[reord]
		if (itmp==1) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==2) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==3) colnames(out)<-c("par.num.","logit-g","s.e.","par.num.","c","s.e.","par.num.","a","s.e.")
		if (Rasch) out[,c(4,6)]<-NA
		rownames(out)<-par$V2
		out<-round(out,digits)
		print(out)
	}
	return(list(coef=p,var=vm))
}


import.irtpro<-function(fnamep,fnamev=NULL,fnameirt=NULL,display=TRUE,digits=2) {
	par<-read.table(fnamep,fill=TRUE)
	ngr<-sum(par$V3==0,na.rm=TRUE)
	if (ngr>1) stop("Cannot handle multiple groups.")
	if(max(par$V2,na.rm=TRUE)>1) stop("Cannot handle multiple factors.")
	if(max(par$V4,na.rm=TRUE)>2) stop("Cannot handle multiple response models.")
	par<-par[par$V3!=0 & !is.na(par$V3),]
	if(length(unique(par$V3))>1) stop("Cannot handle mixed item types.")
	if(unique(par$V3)==2) itmp=2
	if(unique(par$V3)==1) itmp=3
	if(length(unique(par$V5))==1) itmp=1
	ifelse(all(par$V5==1),Rasch<-TRUE,Rasch<-FALSE)
	if (itmp==1 | itmp==2) {p<-par[,c(6,5)]
		colnames(p)<-c("c","a")}
	if (itmp==3) {p<-par[,c(7,6,5)]
		colnames(p)<-c("g","c","a")}
	rownames(p)<-par$V1
	p<-as.matrix(p)

	if (!is.null(fnamev)) {
		if (is.null(fnameirt)) stop("fnameirt required if fnamev is not NULL.")
		vm<-read.table(fnamev,sep=",")
		vm<-as.matrix(vm[,colSums(is.na(vm))!=nrow(vm)])
		colnames(vm)<-NULL
		n<-nrow(p)
		if (itmp==1 & !Rasch) reord<-c((n+1):(2*n),1:n)
		if (itmp==1 & Rasch) reord<-1:n
		if (itmp==2) reord<-c((n+1):(2*n),1:n)
		if (itmp==3) reord<-c((2*n+1):(3*n),(n+1):(2*n),1:n)
		irt <- readLines(fnameirt)
		startRow <- grep(" *Item *Label *P# ", irt)
		nm <- unlist(strsplit(irt[startRow], " +"))[-1]
		widths <- strsplit(irt[startRow], "(?<=[A-Za-z#.])(?=\\s)", perl = TRUE)[[1]]
		widths <- sapply(widths, nchar)
		tab <- read.fwf(fnameirt, widths, skip = startRow, na.strings = "----",
			col.names = nm, comment.char = "", sep = ",", 
			nrows = as.numeric(n), strip.white = TRUE, check.names = FALSE)
		param.number <- unlist(tab[, names(tab) == "P#"])
		pn.max<-max(param.number,na.rm=TRUE)
		num.na<-sum(is.na(param.number))
		if (num.na>0 & !Rasch) {
			fixed<-(pn.max+1):(pn.max+num.na)
			param.number[is.na(param.number)]<-fixed
			nf<-length(fixed)
			nvm1<-nrow(vm)+nf
			vm1<-matrix(0,nvm1,nvm1)
			vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
			vm<-vm1
		}
		if (num.na>0 & Rasch) param.number<-param.number[!is.na(param.number)]
		vm<-vm[param.number,param.number]
		if (!Rasch & n*ncol(p)!=nrow(vm)) stop("Number of parameters and dimension of the covariance matrix do not match.")
		if (Rasch & n!=nrow(vm)) stop("Number of parameters and dimension of the covariance matrix do not match.")
		vm<-vm[reord,reord]
	}
	else vm<-NULL
	if (display) {
		out<-matrix(NA,nrow(p),ncol(p)*3)
		out[,seq(2,ncol(out),by=3)]<-p
		se<-diag(vm)^0.5
		if (!is.null(vm)) out[,seq(3,ncol(out),by=3)]<-se
		if (!is.null(vm)) out[,seq(1,ncol(out),by=3)]<-param.number[reord]
		if (itmp==1) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==2) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==3) colnames(out)<-c("par.num.","g","s.e.","par.num.","c","s.e.","par.num.","a","s.e.")
		if (Rasch) out[,c(4,6)]<-NA
		rownames(out)<-par$V1
		out<-round(out,digits)
		print(out)
	}
	return(list(coef=p,var=vm))
}




import.ltm<-function(mod,display=TRUE,digits=4) {
	if (class(mod)=="grm") stop("Cannot handle multiple response models.")
	if (class(mod)=="gpcm") stop("Cannot handle multiple response models.")
	if (class(mod)=="ltm") if (mod$ltst$factors>1) stop("Cannot handle multiple factors.")
	if (class(mod)=="ltm") if (ncol(mod$coef)>2) stop("Cannot handle not IRT models.")
	p<-mod$coef
	vm<-solve(mod$hessian)
	if (!is.null(mod$constraint) & class(mod)=="rasch") {
		if (!all(p[,ncol(p)]==1)) {
			cstr<-mod$constraint
			nf<-nrow(cstr)
			fixed<-cstr[,1]
			nvm1<-nrow(vm)+nf
			vm1<-matrix(0,nvm1,nvm1)
			vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
			vm<-vm1
		}
	}
	if (!is.null(mod$constraint) & (class(mod)=="ltm" | class(mod)=="tpm")) {
		cstr<-mod$constraint
		nf<-nrow(cstr)
		fixed<-cstr[,1]+(cstr[,2]-1)*nrow(p)
		nvm1<-nrow(vm)+nf
		vm1<-matrix(0,nvm1,nvm1)
		vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
		vm<-vm1
	}
	if (length(unique(p[,ncol(p)]))==1 & !all(p[,ncol(p)]==1)) {
		nv<-nrow(vm)-1
		vm<-vm[c(1:nv,rep(nv+1,nrow(p))),c(1:nv,rep(nv+1,nrow(p)))]
	}
	if (display) {
		out<-matrix(NA,nrow(p),ncol(p)*2)
		out[,seq(1,ncol(out),by=2)]<-p
		se<-diag(vm)^0.5
		if (!is.null(vm)) out[,seq(2,ncol(out),by=2)]<-se
		if (all(p[,ncol(p)]==1)) out[,ncol(out)]<-0
		rownames(out)<-rownames(p)
		colnames(out)<-1:ncol(out)
		colnames(out)[seq(1,ncol(out),by=2)]<-colnames(p)
		colnames(out)[seq(2,ncol(out),by=2)]<-rep("s.e.",ncol(p))
		out<-round(out,digits)
		print(out)
		cat("\nNOTE: Use the modIRT function to transform parameters in the usual IRT parameterization")
		cat("\n")
	}
	return(list(coef=p,var=vm))
}



import.mirt<-function(mod,display=TRUE,digits=3) {
	ngroups<-extract.mirt(mod,what="ngroups")
	if (ngroups>1) stop("Cannot handle multiple groups.")
	nfact<-extract.mirt(mod,what="nfact")
	if (nfact>1) stop("Cannot handle multiple factors.")
	itemtype<-extract.mirt(mod,what="itemtype")	
	if (length(unique(itemtype))>1) stop("Cannot handle mixed item types.")
	itemtype<-itemtype[1]
	if (itemtype!="Rasch" & itemtype!="2PL" & itemtype!="3PL") stop(paste("Cannot handle the",itemtype,"model."))
	constrain<-extract.mirt(mod,what="constrain")
	ni<-extract.mirt(mod,what="nitems")
	itemnames<-extract.mirt(mod,"itemnames")
	par<-mod2values(mod) 
	par<-par[par$item%in%itemnames,]
	if (any(par$class!="dich")) stop("Cannot handle multiple response models.")
	par<-par[,c("item","name","value","est","parnum")]

	par$fixed<-FALSE
	if (any(par$name=="g" & par$value> 0)) par[par$name=="g",]$fixed<-TRUE
	if (any(par$fixed)) itemtype<-"3PL"
	par<-par[par$est==TRUE | par$fixed==TRUE,]

	par$value[par$name=="g"]<-qlogis(par$value[par$name=="g"])
	
	par$parnum<-as.character(par$parnum)

	vcov<-extract.mirt(mod,what="vcov")
	if (all(dim(vcov)>c(1,1))) {
		vcov<-vcov[!apply(is.na(vcov),1,all),]
		vcov<-vcov[,!apply(is.na(vcov),2,all)]
		rnvm<-rownames(vcov)
		rnvm<-sub(pattern="a1.",x=rnvm,replacement="")
		rnvm<-sub(pattern="d.",x=rnvm,replacement="")
		rnvm<-sub(pattern="g.",x=rnvm,replacement="")
		rnvm<-sub(pattern="u.",x=rnvm,replacement="")
		vm<-vcov
		rownames(vm)<-colnames(vm)<-rnvm
		if (any(par$fixed)) {
			rnvm<-rownames(vm)
			nf<-sum(par$fixed)
			rnvm1<-c(rnvm,par$parnum[par$fixed])
			nrvm<-nrow(vm)
			vm1<-matrix(0,nrvm+nf,nrvm+nf)
			vm1[1:nrvm,1:nrvm]<-vm
			colnames(vm1)<-rownames(vm1)<-rnvm1
			vm<-vm1
		}
		if (length(constrain)>0) {
			for (i in 1:length(constrain)) {
				rnvm<-rownames(vm)
				cstr<-constrain[[i]]
				nf<-length(cstr)
				nc<-cstr[1]
				for (j in 2:nf) nc<-paste(nc,cstr[j],sep=".")
				rnvm1<-c(rnvm,rep(nc,nf-1))
				vm<-vm[rnvm1,rnvm1]
				par[par$parnum%in%cstr,]$parnum<-nc
			}
		}
		vm<-vm[par$parnum,par$parnum]
	}
	else vm<-NULL
	if (!is.null(vm)) par$SE<-diag(vm)^0.5
	else par$SE<-NA
	out<-reshape(par,direction="wide",idvar="item",timevar="name",drop=c("est"))
	rownames(out)<-out$item
	out$item<-c()
	pn<-c(out$parnum.g,out$parnum.d,out$parnum.a1)
	vm<-vm[pn,pn]
	rownames(vm)<-colnames(vm)<-NULL
	if (itemtype=="Rasch") out<-out[,c("value.d","SE.d")]
	if (itemtype=="2PL") out<-out[,c("value.d","SE.d","value.a1","SE.a1")]
	if (itemtype=="3PL") out<-out[,c("value.g","SE.g","value.d","SE.d","value.a1","SE.a1")]
	p<-out
	p$SE.g<-c()
	p$SE.d<-c()
	p$SE.a1<-c()
	p<-as.matrix(p)
	if (display) {
		if (!is.null(vm)) out<-round(out,digits)
		if (is.null(vm)) out[,seq(1,ncol(out),by=2)]<-round(out[,seq(1,ncol(out),by=2)],digits)
		print(out)
		cat("\nNOTE: Use the modIRT function to transform parameters in the usual IRT parameterization")
		cat("\n")
	}
	return(list(coef=p,var=vm))
}




itm <- function(x, ...) UseMethod("itm")

itm.eqc<-function(x, ...) x$tab

itm.eqclist<-function(x,link=NULL, ...)
{
	if (is.null(link)) stop("Argument link can not be NULL.")
	return(x[[link]]$tab)
}

itm.ceqc<-function(x, ...) x$tab

itm.ceqclist<-function(x,path=NULL, ...)
{
	if (is.null(path)) stop("Argument path can not be NULL.")
	return(x[[path]]$tab)
}

itm.meqc<-function(x, link=NULL, bistype=NULL, ...)
{
	tmp<-NULL
	one<-TRUE
	bis<-wbis<-FALSE
	if (!is.null(x$bis)) bis<-TRUE
	if (!is.null(x$wbis)) wbis<-TRUE
	if (bis & wbis) one<-FALSE
	if (one) {
		if (is.null(bistype)) {if (bis) tmp<-x$bis else tmp<-x$wbis}
		if (!is.null(bistype)) {if (bistype=="unweighted") tmp<-x$bis; if (bistype=="weighted") tmp<-x$wbis}
		if (is.null(tmp)) stop("bistype does not match.")
	}
	if (!one){
		if (is.null(bistype)) stop("Speficy argument bistype.")
		else {if (bistype=="unweighted") tmp<-x$bis; if (bistype=="weighted") tmp<-x$wbis}
	}
	if (length(tmp)>1 & is.null(link)) stop("Argument link can not be NULL.")
	if (is.null(link)) return(tmp[[1]]$tab)
	else return(tmp[[link]]$tab)
}


eqc <- function(x, ...) UseMethod("eqc")

eqc.eqc<-function(x, ...)
{
	link<-x$forms
	A1<-x$A
	B1<-x$B
	out<-data.frame(link=link,A=A1,B=B1)
	rownames(out)<-NULL
	return(out)
}

eqc.eqclist<-function(x,link=NULL, ...)
{
	link1<-sapply(x,FUN=function(x) x$forms)
	A1<-sapply(x,FUN=function(x) x$A)
	B1<-sapply(x,FUN=function(x) x$B)
	out<-data.frame(link=link1,A=A1,B=B1)
	if (is.null(link)) link<-out$link
	out<-out[out$link%in%link,]
	rownames(out)<-NULL
	return(out)
}

eqc.ceqc<-function(x, ...)
{
	path<-x$forms
	A1<-x$A
	B1<-x$B
	out<-data.frame(path=path,A=A1,B=B1)
	rownames(out)<-NULL
	return(out)
}

eqc.ceqclist<-function(x,link=NULL,path=NULL, ...)
{
	path1<-sapply(x,FUN=function(x) x$forms)
	link1<-path2link(path1)
	A1<-sapply(x,FUN=function(x) x$A)
	B1<-sapply(x,FUN=function(x) x$B)
	out<-data.frame(link=link1,path=path1,A=A1,B=B1)
	if (is.null(link)) link<-out$link
	if (is.null(path)) path<-out$path
	out<-out[out$link%in%link & out$path%in%path,]
	rownames(out)<-NULL
	return(out)
}

eqc.meqc<-function(x,link=NULL,path=NULL, ...)
{
	if (is.null(link)) link<-x$coef$link
	if (is.null(path)) path<-x$coef$path
	out<-x$coef[x$coef$link%in%link & x$coef$path%in%path,1:4]
	rownames(out)<-NULL
	return(out)
}


score<-function(obj, link=NULL, path=NULL, method="TSE", D=1, scores=NULL, se=TRUE, bistype=NULL, nq=30, w=0.5, theta=NULL, weights=NULL)
{
	if (class(obj)=="ceqclist" & is.null(path) & length(obj)==1) path<-obj[[1]]$forms
	if (class(obj)=="ceqclist" & is.null(path)) stop("Specify path.")
	if (class(obj)=="eqclist" & is.null(link)) stop("Specify link.")
	if (class(obj)=="meqc" & is.null(bistype) & is.null(obj$wbis) & !is.null(obj$bis)) bistype<-"unweighted"
	if (class(obj)=="meqc" & is.null(bistype) & !is.null(obj$wbis) & !is.null(obj$bis)) bistype<-"weighted"
	if (class(obj)=="meqc" & is.null(link) & length(unique(obj$coef$link))==1) link<-unique(obj$coef$link)
	if (class(obj)=="meqc" & is.null(link)) stop("Specify link.")
	if (!is.null(scores)) if (any(round(scores)!=scores)) stop("Scores should be integer values.")
	itmpar<-itm(obj,link=link,path=path,bistype=bistype)
	Item<-itmpar$Item
	# conversion from Y to X
	diffY<-itmpar[,2][substr(Item,1,6)=="Dffclt"]
	discrY<-itmpar[,2][substr(Item,1,6)=="Dscrmn"]
	guessY<-itmpar[,2][substr(Item,1,6)=="Gussng"]
	
	diffX<-itmpar[,3][substr(Item,1,6)=="Dffclt"]
	discrX<-itmpar[,3][substr(Item,1,6)=="Dscrmn"]
	guessX<-itmpar[,3][substr(Item,1,6)=="Gussng"]

	diffY2X<-itmpar[,4][substr(Item,1,6)=="Dffclt"]
	discrY2X<-itmpar[,4][substr(Item,1,6)=="Dscrmn"]
	guessY2X<-itmpar[,4][substr(Item,1,6)=="Gussng"]

	if (class(obj)=="eqc") obj1<-obj
	if (class(obj)=="ceqc") obj1<-obj
	if (class(obj)=="ceqclist") obj1<-obj[[path]]
	if (class(obj)=="eqclist") obj1<-obj[[link]]
	if (class(obj)=="meqc") {
		if (is.null(bistype)) {
			if (!is.null(obj$wbis[[link]])) obj1<-obj$wbis[[link]] else obj1<-obj$bis[[link]]
		}
		if (!is.null(bistype)) {
			if (bistype=="weighted") obj1<-obj$wbis[[link]]
			if (bistype=="unweighted") obj1<-obj$bis[[link]]
		}
	}
	names(diffY)<-paste(Item[substr(Item,1,6)=="Dffclt"],obj1$suff[1],sep="")
	if (length(discrY)>0) names(discrY)<-paste(Item[substr(Item,1,6)=="Dscrmn"],obj1$suff[1],sep="")
	if (length(guessY)>0) names(guessY)<-paste(Item[substr(Item,1,6)=="Gussng"],obj1$suff[1],sep="")
	names(diffX)<-paste(Item[substr(Item,1,6)=="Dffclt"],obj1$suff[length(obj1$suff)],sep="")
	if (length(discrX)>0) names(discrX)<-paste(Item[substr(Item,1,6)=="Dscrmn"],obj1$suff[length(obj1$suff)],sep="")
	if (length(guessX)>0) names(guessX)<-paste(Item[substr(Item,1,6)=="Gussng"],obj1$suff[length(obj1$suff)],sep="")
	
	diffY<-diffY[!is.na(diffY)]
	discrY<-discrY[!is.na(discrY)]
	guessY<-guessY[!is.na(guessY)]
	diffX<-diffX[!is.na(diffX)]
	discrX<-discrX[!is.na(discrX)]
	guessX<-guessX[!is.na(guessX)]
	diffY2X<-diffY2X[!is.na(diffY2X)]
	discrY2X<-discrY2X[!is.na(discrY2X)]
	guessY2X<-guessY2X[!is.na(guessY2X)]

	rX<-length(diffX)
	rY<-length(diffY)

	if (length(discrY)==0) discrY<-rep(1,rY)
	if (length(guessY)==0) guessY<-rep(0,rY)
	if (length(discrX)==0) discrX<-rep(1,rX)
	if (length(guessX)==0) guessX<-rep(0,rX)
	if (length(discrY2X)==0) discrY2X<-rep(1,rY)
	if (length(guessY2X)==0) guessY2X<-rep(0,rY)
	
	if (method=="TSE") {
		scoresall<-0:(rX)
		if (!is.null(scores)) scores<-sort(scores)
		if (is.null(scores)) scores<-scoresall
		outbound<-scores[!scores%in%scoresall]
		if (any(!scores%in%scoresall)) cat("The following scores are out of range:", outbound, "\n")
		scores<-scores[scores%in%scoresall]
		smin<-truescore(person.par=-100,diff=diffX,discr=discrX,guess=guessX,Item=Item,D=D)
		undermin<-scores[scores<smin]
		if (any(scores<smin)) cat("The following scores are not attainable:", undermin, "\n")
		scores<-scores[scores>=smin]
		st<-seq(-3,3,l=(rX+1)) # starting values
		st<-st[scoresall%in%scores]
		theta1<-st
		
		cond<-TRUE
		while (cond) {
			theta2<-theta1-(scores-truescore(person.par=theta1,diff=diffX,discr=discrX,guess=guessX,Item=Item,D=D))/
				partialtruescore(person.par=theta1,diff=diffX,discr=discrX,guess=guessX,Item=Item,D=D)
			theta2[theta2==Inf]<-2
			theta2[theta2==-Inf]<- -2
			cond<-max(abs(theta2-theta1))>0.000001
			theta1<-theta2
		}
		theta<-theta2
		ts<-truescore(person.par=theta,diff=diffY2X,discr=discrY2X,guess=guessY2X,Item=Item,D=D)
		out<-data.frame(theta,scores,ts)
		colnames(out)<-c("theta",colnames(itmpar)[3:4])
		rownames(out)<-NULL
		
		if (se) {
			# computation of standar errors (Ogasawara, JEBS, 2001)
			varFull<-obj1$varFull
			if (is.null(varFull)) stop("The asymptotic covariance matrix of item parameters is necessary to compute standard errors. Set se=FALSE.")
			acovalpha<-varFull[[1]]
			for (j in 2:length(varFull)) {
				acovalpha<-blockdiag(acovalpha,varFull[[j]])
			}
			dAB_dalpha<-matrix(0,nrow(acovalpha),2)
			rownames(dAB_dalpha)<-rownames(acovalpha)
			sel<-rownames(dAB_dalpha)%in%rownames(obj1$partial)
			dAB_dalpha[sel,]<-obj1$partial[rownames(dAB_dalpha)[sel],]
			acovAB<-obj1$varAB
			acovABalpha<-t(dAB_dalpha)%*%acovalpha
			acovbeta<-rbind(cbind(acovalpha,t(acovABalpha)),cbind(acovABalpha,acovAB))
			
			K<-length(theta)
			se<-rep(NA,length(theta))
			for (k in 1:K) {
				deta_dtheta<-0
				deta_dalphaY_a<-matrix(NA,rY,1)
				deta_dalphaY_b<-matrix(NA,rY,1)
				deta_dalphaY_c<-matrix(NA,rY,1)
				deta_dAB<-matrix(0,2,1)
				for (i in 1:rY) {
					PYg<-irtp1(ab=theta[k],diff=diffY2X[i],discr=discrY2X[i],guess=guessY2X[i],D=D)
					deta_dtheta<-deta_dtheta+(1-PYg)*(PYg-guessY[i])*D*discrY[i]/(obj1$A*(1-guessY[i]))
					deta_dalphaY_a[i,1]<- (1-PYg)/(1-guessY[i])*(PYg-guessY[i])*D*(theta[k]-obj1$A*diffY[i]-obj1$B)/obj1$A
					deta_dalphaY_b[i,1]<- (1-PYg)/(1-guessY[i])*(PYg-guessY[i])*D*(-discrY[i])
					deta_dalphaY_c[i,1]<- (1-PYg)/(1-guessY[i])
					deta_dAB[1,1]<- deta_dAB[1,1]-(1-PYg)*(PYg-guessY[i])*D*discrY[i]/(1-guessY[i])*(theta[k]-obj1$B)/obj1$A^2
					deta_dAB[2,1]<- deta_dAB[2,1]-(1-PYg)*(PYg-guessY[i])*D*discrY[i]/(1-guessY[i])/obj1$A
				}
				if (obj1$itmp==3) deta_dalphaY<-rbind(deta_dalphaY_b,deta_dalphaY_a,deta_dalphaY_c)
				if (obj1$itmp==2) deta_dalphaY<-rbind(deta_dalphaY_b,deta_dalphaY_a)
				if (obj1$itmp==1) deta_dalphaY<-deta_dalphaY_b
				rownames(deta_dalphaY)<-c(names(diffY),names(discrY),names(guessY))
				rownames(deta_dAB)<-c("A","B")

				denom<-0 # denom of dtheta_dalphaX
				dtheta_dalphaX_a<-matrix(NA,rX,1)
				dtheta_dalphaX_b<-matrix(NA,rX,1)
				dtheta_dalphaX_c<-matrix(NA,rX,1)
				for (i in 1:rX) {
					PXg<-irtp1(ab=theta[k],diff=diffX[i],discr=discrX[i],guess=guessX[i],D=D)
					denom<-denom+(1-PXg)*(PXg-guessX[i])*D*discrX[i]/((1-guessX[i]))
					dtheta_dalphaX_a[i,1]<- -(1-PXg)/(1-guessX[i])*(PXg-guessX[i])*D*(theta[k]-diffX[i])
					dtheta_dalphaX_b[i,1]<- -(1-PXg)/(1-guessX[i])*(PXg-guessX[i])*D*(-discrX[i])
					dtheta_dalphaX_c[i,1]<- -(1-PXg)/(1-guessX[i])
				}
				dtheta_dalphaX_a<-dtheta_dalphaX_a/denom
				dtheta_dalphaX_b<-dtheta_dalphaX_b/denom
				dtheta_dalphaX_c<-dtheta_dalphaX_c/denom
				if (obj1$itmp==3) dtheta_dalphaX<-rbind(dtheta_dalphaX_b,dtheta_dalphaX_a,dtheta_dalphaX_c)
				if (obj1$itmp==2) dtheta_dalphaX<-rbind(dtheta_dalphaX_b,dtheta_dalphaX_a)
				if (obj1$itmp==1) dtheta_dalphaX<-dtheta_dalphaX_b
				deta_dalphaX<-deta_dtheta*dtheta_dalphaX
				if (obj1$itmp==3) rownames(deta_dalphaX)<-c(names(diffX),names(discrX),names(guessX))
				if (obj1$itmp==2) rownames(deta_dalphaX)<-c(names(diffX),names(discrX))
				if (obj1$itmp==1) rownames(deta_dalphaX)<-names(diffX)
				deta_dbeta<-rbind(deta_dalphaY,deta_dalphaX,deta_dAB)
				if(any(!rownames(deta_dbeta)[1:(nrow(deta_dbeta)-2)]%in%rownames(acovbeta))) stop("Names mismatch.")
				deta_dbeta1<-matrix(0,nrow(acovbeta),1)
				rownames(deta_dbeta1)<-rownames(acovbeta)
				deta_dbeta1[rownames(deta_dbeta),]<-deta_dbeta
				se[k]<-sqrt(t(deta_dbeta1)%*%acovbeta%*%deta_dbeta1)
			}
			out$StdErr<-se
		}
	}
	
	if (method=="OSE") {
		if (length(theta)!=length(weights)) stop("Theta and weights should have the same length.")
		if (!is.null(theta) & is.null(weights)) stop("Specify weights or delete theta.")
		if (is.null(theta) & !is.null(weights)) stop("Specify theta or delete weights.")
		if (is.null(theta) & is.null(weights)) {
			gq<-gauss.quad.prob(n=nq,dist="normal")
			theta<-gq$nodes
			weights<-gq$weights
		}
		lt<-length(theta)
		
		# pop 1
		prX1<-matrix(NA,lt,rX)
		prY1<-matrix(NA,lt,rY)
		for (i in 1:rX) prX1[,i]<-irtp1(ab=theta,diff=diffX[i],discr=discrX[i],guess=guessX[i],D=D)
		for (i in 1:rY) prY1[,i]<-irtp1(ab=theta,diff=diffY2X[i],discr=discrY2X[i],guess=guessY2X[i],D=D)
		f1X_theta<-matrix(NA,lt,rX+1)
		f1Y_theta<-matrix(NA,lt,rY+1)
		for (i in 1:length(theta)) f1X_theta[i,]<-fr(r=rX,pr=prX1[i,])
		for (i in 1:length(theta)) f1Y_theta[i,]<-fr(r=rY,pr=prY1[i,])
		mweightX<-matrix(rep(weights,rX+1),ncol=rX+1)
		mweightY<-matrix(rep(weights,rY+1),ncol=rY+1)
		f1X<-colSums(f1X_theta*mweightX)
		f1Y<-colSums(f1Y_theta*mweightY)

		# pop 2
		prX2<-matrix(NA,lt,rX)
		prY2<-matrix(NA,lt,rY)
		for (i in 1:rX) prX2[,i]<-irtp1(ab=theta*obj1$A+obj1$B,diff=diffX[i],discr=discrX[i],guess=guessX[i],D=D)
		for (i in 1:rY) prY2[,i]<-irtp1(ab=theta*obj1$A+obj1$B,diff=diffY2X[i],discr=discrY2X[i],guess=guessY2X[i],D=D) #sarebbe lo stesso usare theta non trasformato e item parameter non trasformati
		rX<-ncol(prX2)
		rY<-ncol(prY2)
		f2X_theta<-matrix(NA,lt,rX+1)
		f2Y_theta<-matrix(NA,lt,rY+1)
		for (i in 1:length(theta)) f2X_theta[i,]<-fr(r=rX,pr=prX2[i,])
		for (i in 1:length(theta)) f2Y_theta[i,]<-fr(r=rY,pr=prY2[i,])
		mweightX<-matrix(rep(weights,rX+1),ncol=rX+1)
		mweightY<-matrix(rep(weights,rY+1),ncol=rY+1)
		f2X<-colSums(f2X_theta*mweightX)
		f2Y<-colSums(f2Y_theta*mweightY)

		# synthetic population
		fX<-f1X*w+f2X*(1-w)
		fY<-f1Y*w+f2Y*(1-w)
		
		FX <- cumsum(fX)
		FY <- cumsum(fY)
		FX<-c(0,FX)
		FY<-c(0,FY)
		
		Px<-c()
		for (j in 2:(rX+2)) Px <- c(Px, FX[j-1]+(j-(j-.5))*(FX[j]-FX[j-1]))
		
		FY<-FY[-1]

		eY<-c()
		eY1<-c()
		for (xx in 1:(rX+1)) {
			Ps<-Px[xx]
			yU<-min((0:rY)[FY>Ps])
			GyU<-FY[yU+1]
			if (yU==0) GyUm1<-0 else GyUm1<-FY[yU]
			eY<-c(eY,(Ps-GyUm1)/(GyU-GyUm1)+(yU-0.5))
			if (yU==0) yL<- -1 else yL<-max((0:rY)[FY<Ps])
			GyLp1<-FY[yL+2]
			if (yU==0) GyL<-0 else GyL<-FY[yL+1]
			eY1<-c(eY1,(Ps-GyL)/(GyLp1-GyL)+(yL+0.5))
		}
		eY<-apply(cbind(eY,eY1),1,mean)
		xx<-0:rX
		
		out<-data.frame(xx,eY)
		colnames(out)<-c(colnames(itmpar)[3:4])
		if (is.null(scores)) scores<-xx

		# computation of standar errors (Ogasawara, Psychometrika, 2003)
		if (se) {
			varFull<-obj1$varFull
			if (is.null(varFull)) stop("The asymptotic covariance matrix of item parameters is necessary to compute standard errors. Set se=FALSE.")
			acovalpha<-varFull[[1]]
			for (j in 2:length(varFull)) {
				acovalpha<-blockdiag(acovalpha,varFull[[j]])
			}
			dAB_dalpha<-matrix(0,nrow(acovalpha),2)
			rownames(dAB_dalpha)<-rownames(acovalpha)
			sel<-rownames(dAB_dalpha)%in%rownames(obj1$partial)
			dAB_dalpha[sel,]<-obj1$partial[rownames(dAB_dalpha)[sel],]
			acovAB<-obj1$varAB
			acovABalpha<-t(dAB_dalpha)%*%acovalpha
			acovbeta<-rbind(cbind(acovalpha,t(acovABalpha)),cbind(acovABalpha,acovAB))
		
			Da1X<-Db1X<-Dc1X<-list()
			Da1Y<-Db1Y<-Dc1Y<-list()
			Da2X<-Db2X<-Dc2X<-list()
			Da2Y<-Db2Y<-Dc2Y<-list()
			DA1Y<-DB1Y<-DA2X<-DB2X<-0
			for (j in 1:rX) {
				f1X_theta_uno<-matrix(NA,lt,rX+1)
				for (i in 1:length(theta)) f1X_theta_uno[i,]<-fr_cond(r=rX,pr=prX1[i,],which=j,value=1)
				f1X_theta_zero<-matrix(NA,lt,rX+1)
				for (i in 1:length(theta)) f1X_theta_zero[i,]<-fr_cond(r=rX,pr=prX1[i,],which=j,value=0)
				f1X_theta_diff<-f1X_theta_uno-f1X_theta_zero
				
				f2X_theta_uno<-matrix(NA,lt,rX+1)
				for (i in 1:length(theta)) f2X_theta_uno[i,]<-fr_cond(r=rX,pr=prX2[i,],which=j,value=1)
				f2X_theta_zero<-matrix(NA,lt,rX+1)
				for (i in 1:length(theta)) f2X_theta_zero[i,]<-fr_cond(r=rX,pr=prX2[i,],which=j,value=0)
				f2X_theta_diff<-f2X_theta_uno-f2X_theta_zero
				
				da1X<-(prX1[,j]-guessX[j])*(1-prX1[,j])*D*(theta-diffX[j])/(1-guessX[j])
				db1X<-(prX1[,j]-guessX[j])*(1-prX1[,j])*D*(-discrX[j])/(1-guessX[j])
				dc1X<-(1-prX1[,j])/(1-guessX[j])
				da1X<-matrix(rep(da1X,rX+1),ncol=rX+1)
				db1X<-matrix(rep(db1X,rX+1),ncol=rX+1)
				dc1X<-matrix(rep(dc1X,rX+1),ncol=rX+1)
				da1X<-da1X*f1X_theta_diff
				db1X<-db1X*f1X_theta_diff
				dc1X<-dc1X*f1X_theta_diff
				Da1X[[j]]<-da1X
				Db1X[[j]]<-db1X
				Dc1X[[j]]<-dc1X
				
				da2X<-(prX2[,j]-guessX[j])*(1-prX2[,j])*D*(obj1$A*theta+obj1$B-diffX[j])/(1-guessX[j])
				db2X<-(prX2[,j]-guessX[j])*(1-prX2[,j])*D*(-discrX[j])/(1-guessX[j])
				dc2X<-(1-prX2[,j])/(1-guessX[j])
				da2X<-matrix(rep(da2X,rX+1),ncol=rX+1)
				db2X<-matrix(rep(db2X,rX+1),ncol=rX+1)
				dc2X<-matrix(rep(dc2X,rX+1),ncol=rX+1)
				da2X<-da2X*f2X_theta_diff
				db2X<-db2X*f2X_theta_diff
				dc2X<-dc2X*f2X_theta_diff
				dA2X<-(prX2[,j]-guessX[j])*(1-prX2[,j])*D*(discrX[j]*theta)/(1-guessX[j])
				dB2X<-(prX2[,j]-guessX[j])*(1-prX2[,j])*D*(discrX[j])/(1-guessX[j])
				dA2X<-matrix(rep(dA2X,rX+1),ncol=rX+1)
				dB2X<-matrix(rep(dB2X,rX+1),ncol=rX+1)
				dA2X<-dA2X*f2X_theta_diff
				dB2X<-dB2X*f2X_theta_diff
				DA2X<-DA2X+dA2X
				DB2X<-DB2X+dB2X
				Da2X[[j]]<-da2X
				Db2X[[j]]<-db2X
				Dc2X[[j]]<-dc2X
			}
			for (j in 1:rY) {
				f1Y_theta_uno<-matrix(NA,lt,rY+1)
				for (i in 1:length(theta)) f1Y_theta_uno[i,]<-fr_cond(r=rY,pr=prY1[i,],which=j,value=1)
				f1Y_theta_zero<-matrix(NA,lt,rY+1)
				for (i in 1:length(theta)) f1Y_theta_zero[i,]<-fr_cond(r=rY,pr=prY1[i,],which=j,value=0)
				f1Y_theta_diff<-f1Y_theta_uno-f1Y_theta_zero
				
				f2Y_theta_uno<-matrix(NA,lt,rY+1)
				for (i in 1:length(theta)) f2Y_theta_uno[i,]<-fr_cond(r=rY,pr=prY2[i,],which=j,value=1)
				f2Y_theta_zero<-matrix(NA,lt,rY+1)
				for (i in 1:length(theta)) f2Y_theta_zero[i,]<-fr_cond(r=rY,pr=prY2[i,],which=j,value=0)
				f2Y_theta_diff<-f2Y_theta_uno-f2Y_theta_zero
				
				da1Y<-(prY1[,j]-guessY[j])*(1-prY1[,j])*D*((theta-obj1$B)/obj1$A-diffY[j])/(1-guessY[j])
				db1Y<-(prY1[,j]-guessY[j])*(1-prY1[,j])*D*(-discrY[j])/(1-guessY[j])
				dc1Y<-(1-prY1[,j])/(1-guessY[j])
				da1Y<-matrix(rep(da1Y,rY+1),ncol=rY+1)
				db1Y<-matrix(rep(db1Y,rY+1),ncol=rY+1)
				dc1Y<-matrix(rep(dc1Y,rY+1),ncol=rY+1)
				da1Y<-da1Y*f1Y_theta_diff
				db1Y<-db1Y*f1Y_theta_diff
				dc1Y<-dc1Y*f1Y_theta_diff
				dA1Y<-(prY1[,j]-guessY[j])*(1-prY1[,j])*D*(discrY[j]/obj1$A^2*(obj1$B-theta))/(1-guessY[j])
				dB1Y<-(prY1[,j]-guessY[j])*(1-prY1[,j])*D*(-discrY[j]/obj1$A)/(1-guessY[j])
				dA1Y<-matrix(rep(dA1Y,rY+1),ncol=rY+1)
				dB1Y<-matrix(rep(dB1Y,rY+1),ncol=rY+1)
				dA1Y<-dA1Y*f1Y_theta_diff
				dB1Y<-dB1Y*f1Y_theta_diff
				DA1Y<-DA1Y+dA1Y
				DB1Y<-DB1Y+dB1Y
				Da1Y[[j]]<-da1Y
				Db1Y[[j]]<-db1Y
				Dc1Y[[j]]<-dc1Y
				
				da2Y<-(prY2[,j]-guessY[j])*(1-prY2[,j])*D*(theta-diffY[j])/(1-guessY[j])
				db2Y<-(prY2[,j]-guessY[j])*(1-prY2[,j])*D*(-discrY[j])/(1-guessY[j])
				dc2Y<-(1-prY2[,j])/(1-guessY[j])
				da2Y<-matrix(rep(da2Y,rY+1),ncol=rY+1)
				db2Y<-matrix(rep(db2Y,rY+1),ncol=rY+1)
				dc2Y<-matrix(rep(dc2Y,rY+1),ncol=rY+1)
				da2Y<-da2Y*f2Y_theta_diff
				db2Y<-db2Y*f2Y_theta_diff
				dc2Y<-dc2Y*f2Y_theta_diff
				Da2Y[[j]]<-da2Y
				Db2Y[[j]]<-db2Y
				Dc2Y[[j]]<-dc2Y
			}
			evalX<-out[,1]
			evalY<-out[,2]
			DPxa<-DPx(D1=Da1X,D2=Da2X,mweight=mweightX,w=w,xp=evalX)
			DPxb<-DPx(D1=Db1X,D2=Db2X,mweight=mweightX,w=w,xp=evalX)
			DPxc<-DPx(D1=Dc1X,D2=Dc2X,mweight=mweightX,w=w,xp=evalX)
			DPya<-DPx(D1=Da1Y,D2=Da2Y,mweight=mweightY,w=w,xp=evalY)
			DPyb<-DPx(D1=Db1Y,D2=Db2Y,mweight=mweightY,w=w,xp=evalY)
			DPyc<-DPx(D1=Dc1Y,D2=Dc2Y,mweight=mweightY,w=w,xp=evalY)
			DPyA<-DPx(D1=list(DA1Y),D2=list(matrix(0,nrow(DA1Y),ncol(DA1Y))),mweight=mweightY,w=w,xp=evalY)
			DPyB<-DPx(D1=list(DB1Y),D2=list(matrix(0,nrow(DB1Y),ncol(DB1Y))),mweight=mweightY,w=w,xp=evalY)
			DPxA<-DPx(D1=list(matrix(0,nrow(DA2X),ncol(DA2X))),D2=list(DA2X),mweight=mweightX,w=w,xp=evalX)
			DPxB<-DPx(D1=list(matrix(0,nrow(DB2X),ncol(DB2X))),D2=list(DB2X),mweight=mweightX,w=w,xp=evalX)
			colnames(DPxa)<-names(discrX)
			colnames(DPxb)<-names(diffX)
			colnames(DPxc)<-names(guessX)
			colnames(DPya)<-names(discrY)
			colnames(DPyb)<-names(diffY)
			colnames(DPyc)<-names(guessY)
			colnames(DPyA)<-"A"
			colnames(DPyB)<-"B"
			colnames(DPxA)<-"A"
			colnames(DPxB)<-"B"
			
			se<-rep(NA,rX+1)
			for (k in 1:(rX+1)) {
				yp<-out[k,2]
				if (obj1$itmp==3) Der<-c(-DPyc[k,],-DPyb[k,],-DPya[k,],DPxc[k,],DPxb[k,],DPxa[k,],DPxA[k]-DPyA[k],DPxB[k]-DPyB[k])
				if (obj1$itmp==2) Der<-c(-DPyb[k,],-DPya[k,],DPxb[k,],DPxa[k,],DPxA[k]-DPyA[k],DPxB[k]-DPyB[k])
				if (obj1$itmp==1) Der<-c(-DPyb[k,],DPxb[k,],DPxB[k]-DPyB[k])
				Der<-Der/fY[round(yp)+1] # +1 because it starts from zero
				if (obj1$itmp==1) names(Der)[length(Der)]<-"B"
				else names(Der)[(length(Der)-1):length(Der)]<-c("A","B")
				if(any(!names(Der)%in%rownames(acovbeta))) stop("Names mismatch.")
				Der1<-matrix(0,nrow(acovbeta),1)
				rownames(Der1)<-rownames(acovbeta)
				Der1[names(Der),]<-Der
				se[k]<-(t(Der1)%*%acovbeta%*%(Der1))^0.5
			}
			out$StdErr<-se
		}
		outbound<-scores[!scores%in%xx]
		if (any(!scores%in%xx)) cat("The following scores are out of range:", outbound, "\n")
		out<-out[xx%in%scores,]
	}
	return(out)
}

# recursion formula 
fr<-function(r,pr)
{
	if (r==0) stop("r should be at least 1.")
	if (r==1) out<-c(1-pr[r],pr[r])
	else {
		frm1<-fr(r-1,pr=pr) # vector with dimension r
		p1<-frm1[1]*(1-pr[r])
		if (r==2) p2<-frm1[2]*(1-pr[r])+frm1[1]*pr[r]
		else p2<-frm1[2:r]*(1-pr[r])+frm1[1:(r-1)]*pr[r]
		p3<-frm1[r]*pr[r]
		out<-c(p1,p2,p3)
	}
	return(out)
}

# compute true score
truescore<-function(person.par,diff,discr,guess,Item,D)
{
	out<-rep(0,length(person.par))
	ni<-length(diff)
	for (i in 1:ni)
		if (!is.na(diff[i]))
			out<-out+irtp1(ab=person.par,diff=diff[i],discr=discr[i],guess=guess[i],D=D)
	return(out)
}

partialtruescore<-function(person.par,diff,discr,guess,Item,D)
{
	out<-rep(0,length(person.par))
	ni<-length(diff)
	for (i in 1:ni)
		if (!is.na(diff[i])) {
			p<-irtp1(ab=person.par,diff=diff[i],discr=discr[i],guess=guess[i],D=D)
			out<-out+D*discr[i]*(1-p)*(p-guess[i])/(1-guess[i])
		}
	return(-out)
}

fr_cond<-function(r,pr,which,value)
{
	pr[which]<-value
	fr(r,pr)
}

DPx<-function(D1,D2,mweight,w,xp)
{
	out<-c()
	for (k in 1:length(D1)) {
		DfX<-colSums(D1[[k]]*mweight)*w+colSums(D2[[k]]*mweight)*(1-w)
		DPxk<-c()
		for (j in 1:(length(xp))) DPxk<-c(DPxk,Fp(xp=xp[j],fX=DfX))
		out<-cbind(out,DPxk)
	}
	return(out)
}

Fp<-function(xp,fX)
{
	xx<-round(xp)
	wx<-xp-xx+0.5
	if (xx==0) res<-wx*fX[1]
	if (xx>0) res<-sum(fX[1:xx])+wx*fX[xx+1]
	return(res)
}


