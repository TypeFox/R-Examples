Clest <-
function(d,maxK,alpha,B=15,B0=5,nstart=1000,L1=6,beta=0.1,pca=TRUE,silent=FALSE){
        n.rskc<-(maxK-1)*(B0*1*(1+1)+B*(1+1))
	cat("RSKC will be performed (maxK-1)*(B0*1*(1+1)+B*(1+1))=",n.rskc,"times\n")
        if (is.null(L1)){sparse<-FALSE} else{ sparse<-TRUE}
        ### Notes ###	        
        # In total, RSKC is run (MaxK-1)*(B*2+B0*1*2)      
        # dCER = obs CER - ref CER 
        call<-match.call()
        if (pca==TRUE){S<-var(d,na.rm=TRUE) # p by p
        	       P<-eigen(S)$vectors
                       Y<-d%*%P # n by p
                       }else{Y<-d}
        dmax<-apply(Y,2,max,na.rm=TRUE);
        dmin<-apply(Y,2,min,na.rm=TRUE)
        trans.scale=1/(dmax-dmin)
        trans.center=-dmin               
                       
        maxK1<-maxK-1      
        dCER<-rep(Inf,maxK1);n<-nrow(d);f<-ncol(d)
        Rf<-matrix(NA,B0,maxK1);colnames(Rf)<-c("k=2",3:maxK)
        result<-matrix(NA,maxK1,4);Pvalues<-rep(NA,maxK1);I<-2:maxK
        rownames(result)<-I;
        colnames(result)<-c("test.stat","obsCER","refCER","P-value")
        # 1) Clest.base for reference dataset
        # RSKC is performed B0*(maxK-1)*(1*(1+1))
        for ( i in 1 :B0){
              if (silent==FALSE)
              cat("\nAssessing a reference data",i, "out of",B0,"\n")
              NULLd<-reference.data(n,f,trans.center,trans.scale)
              if (pca==TRUE){NULLd<-NULLd%*%t(P)}
              # Clest.base returns vector of CER in [0,1]^(maxK-1) 
              # for reference data, B=1 is enough 
              # B = # partitinoing data into training and testing sets
              Rf[i,]<-Clest.base(NULLd,B=1,alpha,nstart,n,L1,maxK,sparse=sparse,silent=silent)$cer
              }   
         
        refCER<-apply(Rf,2,mean,na.rm=TRUE); 

        # 2) Clest.base for the observed data 
        # RSKC is performed (maxK-1)*(B*(1+1))
        cat("\nAssessing the observed data \n")
        cbase.obs<-Clest.base(d,B,alpha,nstart,n,L1,maxK,sparse=sparse,silent=silent)
        obsCER<-cbase.obs$cer
        observedCERs<-cbase.obs$bootstrap
        for ( ncl in 2:maxK){
            iK<-ncl-1
            krefs<-Rf[,iK];kobsCER<-obsCER[iK];krefCER<-refCER[iK];
            # in Clest (Dudoit and Fridlyand (2002)), the p-values(k) are computed as the proportion of
            # the reference agreements of two partitions that are at least as good as the observed one for each k
            # However, in this computation, the p-values are computed as the proportion of the refCERs that
            # are better than the obsCER. (No equality)
            Pvalues[iK]<-mean(krefs[!is.na(krefs)]<kobsCER)
            dCER[iK]<-(kobsCER-krefCER)
            result[iK,]<-c(dCER[iK],kobsCER,krefCER,Pvalues[iK])
            }
        Im1<-1:(maxK1)
        SigK<-Im1[Pvalues<beta]  
         # if SigK = {2,4} then K=3 and 5 are significant 
        if (length(SigK)==0){K<-1}else{K<-SigK[which.min(dCER[SigK])]+1}
        # take the number of clusters that 
        # corresponds to the largest signifance difference 
        return(list(call=call,K=K,result.table=result,referenceCERs=Rf,observedCERs=observedCERs))
}


reference.data <-
function(n,f,trans.center,trans.scale){ 
        # generate n*f sample pts from uniform(0,1) 
        nulld<-matrix(runif(n*f,0,1),nrow=n,ncol=f)
        # transform i^th column to unif(dmin_i,dmax_i)
        # X follows U(0,1)
        # P(X<x)=x => P(X<(x-min)/(max-min))=(x-min)/(max-min)
        # P((max-min)*X+min<x)=(x-min)/(max-min)
        # Y=(max-min)*X+min follow unif(min,max)
        # note that:
        # trans.scale=1/(dmax-dmin)
        # trans.center=-dmin
        nulld<-scale(nulld,center=FALSE,scale=trans.scale)
        nulld<-scale(nulld,center=trans.center,scale=FALSE)
        return(nulld)
        }

