est_multilogit <-
function(Y,Xdis,label=1:n,be=NULL,Pdis=NULL,dis=FALSE,fort=TRUE,ex=FALSE){
		
# fit multinomial logit model for reponses included in Y
# covariates in matrix X (two formats are allows for X: with binary responses it is n x ncov, with more response categories it is ncat x ncov x ndis)
# label, in case of redundacies are the labels say to which design matrix in Xdis each row of Y is referred		
# ex = TRUE if exit after only computing score and information		
		
time = proc.time()
# preliminaries
     n = dim(Y)[1]
 	 k = dim(Y)[2]
     nbe = dim(Xdis)[2]
     ndis = max(label)
# correct covariance matrix
     if(length(dim(Xdis))==2) Xdis = aperm(array(Xdis,c(ndis,nbe,1),c(3,2,1)))    
     if(dim(Xdis)[1]<k){
     	Xdis0 = Xdis
     	Xdis = array(0,c(k,nbe,n))
     	Xdis[2:k,,] = Xdis0
     } 
#print(c(4,proc.time()-time))
# starting values
	if(is.null(be)){
		be = rep(0,nbe)
		Pdis = NULL}
	if(is.null(Pdis)){
		Pdis = prob_multilogit(Xdis,be,label,fort)$Pdis
	}
    Ydis = matrix(0,ndis,k)	
	if(fort==F){
	    for(i in 1:ndis){
		    li = (label==i)
		    if(sum(li)==1) Ydis[i,] = Y[li,] else Ydis[i,] = colSums(Y[li,])
	    }
	}else{
        out = .Fortran("sum_Y",Ydis=Ydis,Y=Y,label=as.integer(label),ndis=as.integer(ndis),ns=as.integer(n),k=as.integer(k))
        Ydis = out$Ydis	
    }
    ny = rowSums(Ydis)
    lk = sum(Ydis*log(Pdis))
    if(dis) print(c(0,lk))
# iterate until convergence
    it = 0; lko = lk
    while((lk-lko>10^-6 & it<100) | it==0){
        it = it+1; lko = lk 
    	if(fort==F){
        	sc = 0; Fi = 0
	        for(i in 1:ndis){
		        pdis = Pdis[i,]
		        sc = sc+t(Xdis[,,i])%*%(Ydis[i,]-ny[i]*pdis)
		        Fi = Fi+ny[i]*t(Xdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%Xdis[,,i]
	        }
	    }else{
            out = .Fortran("nr_multilogit",Xdis=Xdis,be=be,Pdis=Pdis,Ydis=Ydis,ny=ny,k=as.integer(k),ndis=as.integer(ndis),ncov=as.integer(nbe),
            sc=rep(0,nbe),Fi=matrix(0,nbe,nbe))
            sc = out$sc; Fi = out$Fi
        }
        if(ex==FALSE){
        	dbe = as.vector(ginv(Fi)%*%sc)
	        mdbe = max(abs(dbe))
	        if(mdbe>0.5) dbe = dbe/mdbe*0.5
	        be0 = be
	        flag = TRUE
	        while(flag){
	        	be = be0+dbe
	        	Eta = matrix(0,ndis,k)
	        	for(i in 1:ndis) Eta[i,] = Xdis[,,i]%*%be
	        	if(max(abs(Eta))>100){
	        		dbe = dbe/2
	        		flag = TRUE	
	        	}else{
	        		flag = FALSE
	        	}	        	
	        }
	    }
# compute again probabilities
	        out = prob_multilogit(Xdis,be,label,fort)
	        P = out$P; Pdis = out$Pdis
            lk = sum(Ydis*log(Pdis))
            if(dis) print(c(it,lk,lk-lko))
	}
	out = list(be=be,P=P,Pdis=Pdis,sc=sc,Fi=Fi)
}
