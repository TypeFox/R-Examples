LBRecap=function(data,
 last.column.count=FALSE,
 neval=1000,
 startadd=0,
 by.incr=1,
 mbc.function=c("standard","markov","counts","integer","counts.integer"),
 mod=c("linear.logistic","M0","Mb","Mc","Mcb","Mt","Msubjective.cut","Msubjective"),
 heterogeneity = FALSE,
 markov.ord=NULL,
 z.cut=c(),
 meaningful.mat.subjective=NULL,
 meaningful.mat.new.value.subjective=NULL,
 td.cov=NULL,
 td.cov.formula="",
 verbose=FALSE,
 graph=FALSE,
 output=c("base","complete")){

 mod.est=NULL
 pH.hat=NULL
 logistic=0
 failure.conditions=T
 mod=match.arg(mod)
 output=match.arg(output)
 mbc.function=match.arg(mbc.function)

 model=switch(mod,
 linear.logistic="mod.linear.logistic",
 M0="mod.M0",
 Mb="mod.Mb",
 Mc="mod.Mc",
 Mcb="mod.Mcb",
 Mt="mod.Mt",
 Msubjective.cut="mod.Msubjective.cut",
 Msubjective="mod.Msubjective")	
 qb.function=switch(mbc.function,
  					standard="qb.standard",
  					markov="qb.markov",
  					integer="qb.integer",
  					counts="qb.count",
  					counts.integer="qb.count.integer")

if((mod!="Msubjective.cut") & length(z.cut)>0){
 	warning("z.cut without mod='Msubjective.cut' will be IGNORED!")
}
   if(!(class(data)=="data.frame" | class(data)=="matrix" | class(data)=="array" | class(data)=="table")){
       stop("input data must be a data.frame or a matrix object or an array")
   }
   if(class(data)=="table" | class(data)=="array"){
   	n.occ=length(dim(data))
mm=matrix(ncol=n.occ,nrow=2^n.occ)

for(i in 1:2^n.occ){
mm[i,]=as.numeric(intToBits(i-1)>0)[1:n.occ]
}

data.vec=c(data)
data[1]=0  # the array count is set to 0 by default since it is unobserved
dd=c()
for(i in 1:length(data)){
	dd=c(dd,rep(mm[i,],data.vec[i]))
		}

data.matrix=matrix(dd,ncol=length(dim(data)),byrow=T)
   }
   if(class(data.matrix)!="matrix"){
       data.matrix=as.matrix(data) 
       }
  if(last.column.count==TRUE){

if (any(data[,ncol(data)]<0)){
    stop("Last column must contain non negative frequencies/counts")
}
  	data=as.matrix(data)
  	data.matrix=matrix(ncol=(ncol(data)-1))

for(i in 1:nrow(data)){

d=rep(data[i,1:(ncol(data)-1)],(data[i,ncol(data)]))
dd=matrix(d,ncol=(ncol(data)-1),byrow=T)
data.matrix=rbind(data.matrix,dd)

}
data.matrix=data.matrix[-1,]
  }

  if(any(data.matrix!=0 & data.matrix!=1)) stop("data must be a binary matrix")
 if(sum(apply(data.matrix,1,sum)==0)){
 	warning("input data argument contains rows with all zeros: these rows will be removed and ignored")
        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }

 t=ncol(data.matrix)
 M=nrow(data.matrix)     

 ### add startadd rows corresponding to *startadd* unobserved units
 ### it is effective only when the input argument startadd is greater than 0
   data.matrix=rbind(data.matrix,matrix(0,nrow=startadd,ncol=t))
 if(model=="mod.Mc" | model=="mod.Mcb"){
 if (markov.ord != round(markov.ord) || markov.ord <= 0 || markov.ord > t) 
        stop("Error: The Markov order must be a non-negative integer smaller than t!")
 mbc.function="markov"
 z.cut=seq(0,1,by=1/(2^markov.ord))
 z.cut[1]=-0.1
      if(model=="mod.Mcb"){
      	z.cut=c(z.cut[1],0,z.cut[-1])
      }
        }
if(qb.function=="qb.standard") recodebinary=quant.binary
if(qb.function=="qb.integer") recodebinary=quant.binary.integer
if(qb.function=="qb.markov") recodebinary=quant.binary.markov
if(qb.function=="qb.count") recodebinary=quant.binary.counts
if(qb.function=="qb.count.integer") recodebinary=quant.binary.counts.integer
 v=c()
 meaningful.mat=matrix(NA,ncol=ncol(data.matrix),nrow=(M+startadd))
 for(i in 1:ncol(data.matrix)){
 	for(j in 1:(M+startadd)){
 	v=data.matrix[j,1:i]
 	v=paste(v,collapse="")
	if(qb.function=="qb.markov"){
      		meaningful.mat[j,i]=recodebinary(v,markov.ord)}
      else{meaningful.mat[j,i]=recodebinary(v)} 	}
 }

 meaningful.mat=cbind(rep(0,(M+startadd)),meaningful.mat[,1:t-1])
 meaningful.mat.new.value=rep(0,t*by.incr)
 if(model=="mod.Msubjective"){
   if(is.matrix(meaningful.mat.subjective)){
     if(ncol(meaningful.mat.subjective)!=ncol(data.matrix)||
          nrow(meaningful.mat.subjective)!=nrow(data.matrix))
       stop("meaningful.mat.subjective has different dimensions from data matrix")
     meaningful.mat=meaningful.mat.subjective  
   }
   else{ stop("meaningful.mat.subjective must be a numerical matrix") }
   if(is.numeric(meaningful.mat.new.value.subjective)){
     meaningful.mat.vec.new=meaningful.mat.new.value.subjective
   } 
 }
   if(model=="mod.M0"){
 	z.cut=c(-0.1,1)
 	}
 if(model=="mod.Mb"){
 	z.cut=c(-0.1,0,1)
 	}

 if(model=="mod.Mt"){
 	failure.conditions=F
 	z.cut=seq(0,t)
 	meaningful.mat=matrix(rep(1:t),nrow=(M+startadd),ncol=t,byrow=T)
	meaningful.mat.new.value=1:t
 	}
 if(model=="mod.linear.logistic" | model=="mod.Msubjective" | heterogeneity  ){logistic=1}
 meaningful.mat.vec=c(meaningful.mat)
 y=c(data.matrix)

if(logistic==0){
 p.1=matrix(ncol = (length(z.cut) - 1), nrow = neval)
 }
 nn.1=c()
 nn.0=c()
 n.11=c()
 n.10=c()
 n.01=c()
 n.00=c()
 log.lik=c()
 if(!is.null(td.cov)){
 	logistic=1
        model=paste(paste(model,"+",paste(names(td.cov),collapse=" +")),"[time-dependent covariates included]")
}

  if(heterogeneity){
        model=paste(model,"+ Individual heterogeneity effect")
}

  for(k in 1:neval){
  	dd=NULL
 meaningful.mat.new=rbind(meaningful.mat,matrix(rep(meaningful.mat.new.value,(k>1)),nrow=((k-1) *by.incr),ncol=t,byrow=T))
 dati.new=rbind(data.matrix,matrix(0,nrow=((k-1)*by.incr),ncol=t))
 meaningful.mat.vec.new=c(meaningful.mat.new)
 y.new=c(dati.new) 

if(length(z.cut)>0){	
	if(sum(is.na(cut(meaningful.mat.vec.new,z.cut)))>0){
		stop("z.cut is not correctly specified: NA's produced")	
	}
}

 if(logistic==0){

     for (j in 1:(length(z.cut) - 1)) {
       nn.0[j] = sum(meaningful.mat.vec.new > z.cut[j] & 
       meaningful.mat.vec.new <= z.cut[j + 1] & y.new == 0)
       nn.1[j] = sum(meaningful.mat.vec.new > z.cut[j] & 
       meaningful.mat.vec.new <= z.cut[j + 1] & y.new == 1)
            }

			p.1[k,]=nn.1/(nn.0+nn.1)

            n.00 = nn.0[nn.0 > 0]
            n.10 = nn.1[nn.0 > 0]
            n.01 = nn.0[nn.1 > 0]
            n.11 = nn.1[nn.1 > 0]
            log.lik[k] = lchoose((sum(nn.1) + sum(nn.0))/t, M) + 
                (sum(n.11 * log(n.11/(n.11 + n.01))) + (sum(n.00 * 
                  log(n.00/(n.10 + n.00)))))

}

 if(logistic==1){

	 if(!is.null(td.cov)){
	 	 	if(k>1){
	 	 		names(td.cov.mat) = names(td.cov)
 	         }
	   td.cov.mat=matrix(NA,ncol=ncol(td.cov),nrow=(nrow(td.cov)*(M+startadd+(k-1)*by.incr)))

 for(i in 1:ncol(td.cov.mat)){
	td.cov.mat[,i]=rep(td.cov[,i],each=(M+startadd+(k-1)*by.incr))
	}

	 failure.conditions=F

	 dd=data.frame(y.new,meaningful.mat.vec.new,td.cov.mat)
	 names(dd)=c("Y","Z",paste("X",rep(1:ncol(td.cov.mat)),sep=""))
         formula=paste("Y~Z+",td.cov.formula)

                    }else{
	 failure.conditions=F	
	 dd=data.frame(y.new,meaningful.mat.vec.new)
	 names(dd)=c("Y","Z")
	 formula="Y~Z"
     }
## use the mbc Z as a factor partitioning in terms of z.cut
## cut(dd$Z,breaks=z.cut)
         
         if(length(z.cut>0)){
                        dd$Z=cut(dd$Z,breaks=z.cut)
                        if(mod=="M0") {
                            formula=paste("Y~",td.cov.formula)
                        }
                    }

         if(heterogeneity){
             rand.vec=factor(rep(1:(nrow(dd)/t),t))
             dd=data.frame(dd,rand.vec)
             formula=paste(formula, "+ (1|rand.vec)")
             }

if(is.null(findbars(as.formula(formula)))){
 out.model=glm(formula,binomial(link = "logit"),data=dd)
}else{
     out.model=glmer(formula,binomial(link = "logit"),data=dd,nAGQ=10)
 }

 if( verbose & k%%ceiling(neval/10)==0){
     cat(paste("===>  ", round(k/neval*100),"% of likelihood evaluations\n",sep=""))
 }
 log.lik[k]=as.numeric(logLik(out.model))+lchoose(M+startadd+(k-1)*by.incr,M)
	 	}
}	 	
# closed k cycle

  k.hat=which(log.lik==max(log.lik))
  k=k.hat
  N.hat=M+startadd+(k-1)*by.incr

## redo model fit for k=the optimal k.hat so that
## the model output can be returned

 meaningful.mat.new=rbind(meaningful.mat,matrix(rep(meaningful.mat.new.value,(k>1)),nrow=((k-1) *by.incr),ncol=t,byrow=T))
 dati.new=rbind(data.matrix,matrix(0,nrow=((k-1)*by.incr),ncol=t))
 meaningful.mat.vec.new=c(meaningful.mat.new)
 y.new=c(dati.new) 

if(length(z.cut)>0){	
	if(sum(is.na(cut(meaningful.mat.vec.new,z.cut)))>0){
		stop("z.cut is not correctly specified: NA's produced")	
	}
}

 if(logistic==0){

     for (j in 1:(length(z.cut) - 1)) {
     nn.0[j] = sum(meaningful.mat.vec.new > z.cut[j] & 
     meaningful.mat.vec.new <= z.cut[j + 1] & y.new == 0)
     nn.1[j] = sum(meaningful.mat.vec.new > z.cut[j] & 
     meaningful.mat.vec.new <= z.cut[j + 1] & y.new == 1)
            }

            n.00 = nn.0[nn.0 > 0]
            n.10 = nn.1[nn.0 > 0]
            n.01 = nn.0[nn.1 > 0]
            n.11 = nn.1[nn.1 > 0]
            log.lik[k] = lchoose((sum(nn.1) + sum(nn.0))/t, M) + 
                (sum(n.11 * log(n.11/(n.11 + n.01))) + (sum(n.00 * 
                  log(n.00/(n.10 + n.00)))))

}

 if(logistic==1){

	 if(!is.null(td.cov)){
	 	 	if(k>1){
	 	 		names(td.cov.mat) = names(td.cov)
 	         }
	   td.cov.mat=matrix(NA,ncol=ncol(td.cov),nrow=(nrow(td.cov)*(M+startadd+(k-1)*by.incr)))

 for(i in 1:ncol(td.cov.mat)){
	td.cov.mat[,i]=rep(td.cov[,i],each=(M+startadd+(k-1)*by.incr))
	}

	 failure.conditions=F

	 dd=data.frame(y.new,meaningful.mat.vec.new,td.cov.mat)
	 names(dd)=c("Y","Z",paste("X",rep(1:ncol(td.cov.mat)),sep=""))
         formula=paste("Y~Z+",td.cov.formula)

 if(length(z.cut>0)){
                        dd$Z=cut(dd$Z,breaks=z.cut)
                        if(mod=="M0") {
                            formula=paste("Y~",td.cov.formula)
                        }
                    }

                    }else{
	 failure.conditions=F	
	 dd=data.frame(y.new,meaningful.mat.vec.new)
	 names(dd)=c("Y","Z")
	 formula="Y~Z"
     }

         if(heterogeneity){
             rand.vec=factor(rep(1:(nrow(dd)/t),t))
             dd=data.frame(dd,rand.vec)
             formula=paste(formula, "+ (1|rand.vec)")
             }
if(is.null(findbars(as.formula(formula)))){
 out.model=glm(formula,binomial(link = "logit"),data=dd)
}else{
     out.model=glmer(formula,binomial(link = "logit"),data=dd,nAGQ=10)
 }

 if( verbose & k%%ceiling(neval/10)==0){
     cat(paste("===>  Refit the model in N.hat"))
 }
 mod.est=out.model
 if(any(is.na(coef(mod.est)))){
     warning("There are some NA's in the model coefficients \nIt looks like there are some reduntant model components!")
 }
 log.lik[k]=as.numeric(logLik(out.model))+lchoose(M+(k-1)*by.incr,M)

            npar=length(na.omit(coef(mod.est)))
         if(class(mod.est)[1]=="glmerMod"){
            npar=length(na.omit(fixef(mod.est)))+ncol(ranef(mod.est))
            }
  AIC=-2*as.numeric(logLik(mod.est))+2*(npar+1)-2*lchoose(M+(k-1)*by.incr,M)
	 	}

 if(logistic==0){

 pH.hat=p.1[which(log.lik==max(log.lik)),]
 names(pH.hat)=paste("pH",1:length(pH.hat),sep="")

 AIC=-2*max(log.lik)+2*(length(z.cut))   #### +1 included

}else{

	max.log.lik=as.numeric(logLik(out.model))+lchoose(M+startadd+(k-1)*by.incr,M)
        AIC=-2*as.numeric(logLik(mod.est))+2*(length(na.omit(coef(mod.est)))+1)-2*lchoose(M+startadd+(k-1)*by.incr,M)
	}

  ## CONFIDENCE INTERVAL (Normal approximation)
 z2=(qnorm(0.975))^2
 Nplus=N.hat

 	if(N.hat<=(M+startadd+(neval-1)*(by.incr))){
 i=1

 while((2*(max(log.lik)-log.lik[which(log.lik==(max(log.lik)))+i])<z2 & (N.hat+i*(by.incr))< M+(neval-1)*(by.incr)))
 {
 	i=i+1 
 	}
 	 Nplus=N.hat+(i-1)*by.incr
 	}
 Nminus=N.hat
 	if(N.hat>M+startadd){
 i=1	
 while(2*(max(log.lik)-log.lik[(which(log.lik==(max(log.lik)))-i)])<z2 & (N.hat-i*(by.incr))>M+startadd){
 	i=i+1
 	}
 	 Nminus=N.hat-(i-1)*by.incr
	}

if(model=="mod.Mc"){
	model=paste("Mc",markov.ord,sep="")}

if(model=="mod.Mcb"){
	model=paste("Mc",markov.ord,"b",sep="")}

 if(by.incr>1){
     cat(paste("Likelihood evaluated in the following range: [N.min=",M+startadd,"N.max=",M+startadd,(neval-1)*by.incr,"] namely in 
             seq( from=",M+startadd," , to=",M+startadd+(neval-1)*by.incr," , by=",by.incr," )",sep=""))
     if(startadd>0) cat(paste("\n NOTICE that likelihood evaluations startadded from M+startadd =",M+startadd,"\nrather than from the number M =",M," of observed units",sep=""))
}

Nrange=seq(from=M+startadd,to=M+startadd+(neval-1)*by.incr,by=by.incr)
names(log.lik)=Nrange

 if(graph){
plot(Nrange,log.lik,xlab="N",ylab="logLikelihood(N)",main="logLikelihood evaluations")
points(N.hat,log.lik[as.character(N.hat)],col="red",pch=16)
axis(1,at=N.hat,padj=1)
abline(v=N.hat,lty=5)
 }

 print(paste("Logistic",logistic))
   if(logistic==1) {
 out=switch(output,
 base=list(Model=model,N.hat=N.hat,CI=c(Nminus,Nplus),pH.hat=pH.hat,AIC=AIC,L.Failure=FALSE),
 complete=list(Model=model,N.hat=N.hat,CI=c(Nminus,Nplus),AIC=AIC,L.Failure=FALSE,log.lik=log.lik,N.range=Nrange,meaningful.mat.matrix=meaningful.mat,vec.cut=z.cut,mod.est=mod.est))
 }else{
 n1=sum(y[meaningful.mat.vec<=z.cut[2]])
 n0=length(y[meaningful.mat.vec<=z.cut[2]])-n1
  if((M*(t-1)-n0<=((M-1)*(t-1))/2-1) & n1==M) { 
 out=switch(output,
 base=list(Model=model,N.hat=NA,CI=c(NA,NA),AIC=-Inf,L.Failure=TRUE),
 complete=list(Model=model,N.hat=N.hat,CI=c(Nminus,Nplus),pH.hat=pH.hat,AIC=-Inf,L.Failure=TRUE,log.lik=log.lik,N.range=Nrange,meaningful.mat.matrix=meaningful.mat,vec.cut=z.cut,mod.est=mod.est))
		}else{			
 out=switch(output,
 base=list(Model=model,N.hat=N.hat,CI=c(Nminus,Nplus),AIC=AIC,L.Failure=FALSE),
 complete=list(Model=model,N.hat=N.hat,CI=c(Nminus,Nplus),pH.hat=pH.hat,AIC=AIC,L.Failure=FALSE,log.lik=log.lik,N.range=Nrange,meaningful.mat.matrix=meaningful.mat,vec.cut=z.cut,mod.est=mod.est))
	}
}	 

if(logistic==0){
 if(all(diff(log.lik)>=0)  & failure.conditions==T &
(M*(t-1)-n0<=((M-1)*(t-1))/2-1) & n1==M){
     	cat("\n Likelihood Failure conditions are met! The likelihood is always increasing \n ===> N.hat --> infinity   \n\n")
	}
if(all(diff(log.lik)>=0) & 
((M*(t-1)-n0>((M-1)*(t-1))/2-1) | n1!=M)){
cat("\n Likelihood always increasing within the range of N values evaluated! \n Maximization failed \n Try with a larger upper bound (and use by.incr to skip some evaluations)\n\n")
	}		
}
else{
 if(all(diff(log.lik)>=0)){
     	cat("\n Likelihood always increasing within the range of N values evaluated! \n Maximization failed \n Try with a larger upper bound (and use by.incr to skip some evaluations)\n\n")
	}		
}

  return(out)
 }
