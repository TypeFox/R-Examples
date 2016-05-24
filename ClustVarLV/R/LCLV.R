#' L-CLV for L-shaped data
#' 
#' Define clusters of X-variables aroud latent components. In each cluster, two latent components are extracted, 
#' the first one is a linear combination of the external information collected for the rows of X and the second one
#' is a linear combination of the external information associated with the columns of X.
#' 
#' @param X The matrix of variables to be clustered
#' @param Xr The external variables associated with the rows of X
#' @param Xu The external variables associated with the columns of X
#' @param ccX TRUE/FALSE : double centering of X (FALSE, by default) If FALSE this implies that cX = TRUE : column-centering of X
#' @param sX TRUE/FALSE : standardization or not of the columns X (TRUE by default) 
#' @param sXr TRUE/FALSE : standardization or not of the columns Xr (FALSE by default)\cr
#'        (predefined -> cXr    = TRUE : column-centering of Xr)
#' @param sXu TRUE/FALSE : standardization or not of the columns Xu (FALSE by default)\cr
#'        (predefined -> cXu= FALSE : no centering, Xu considered as a weight matrix)
#' @param nmax maximum number of partitions for which the consolidation will be done (by default nmax=20)
#' @return \item{tabres}{ Results of the clustering algorithm.
#'         In each line you find the results of one specific step of the hierarchical clustering.
#'         \itemize{
#'                \item {Columns 1 and 2}{ : The numbers of the two groups which are merged}
#'                \item {Column 3}{ : Name of the new cluster}
#'                \item {Column 4}{ : The value of the aggregation criterion for the Hierarchical Ascendant Clustering (HAC)}
#'                \item {Column 5}{ : The value of the clustering criterion for the HAC}
#'                \item {Column 6}{ : The percentage of the explained initial criterion value}
#'                \item {Column 7}{ : The value of the clustering criterion after consolidation}
#'                \item {Column 8}{ : The percentage of the explained initial criterion value after consolidation}        
#'                \item {Column 9}{ : number of iterations in the partitioning algorithm.\cr
#'                Remark: A zero in columns 7 to 9 indicates that no consolidation was done }
#'        }}                
#' @return \item{partition K}{ a list for each number of clusters of the partition, K=2 to nmax with
#'          \itemize{
#'                \item {clusters}{ :  in line 1, the groups membership before consolidation; in line 2 the groups membership after consolidation} 
#'                \item {compt}{ : The latent components of the clusters (after consolidation) defined according to the Xr variables}
#'                \item {compc}{ : The latent components of the clusters (after consolidation) defined according to the Xu variables}
#'                \item {loading_v}{ : loadings of the external Xr variables (after consolidation)}
#'                \item {loading_u}{ : loadings of the external Xu variables (after consolidation)}
#'          }}
#' @references Vigneau, E., Endrizzi, I.,& Qannari, E.M. (2011). Finding and explaining clusters of consumers using CLV approach. Food Quality and Preference, 22, 705-713.
#' @references Vigneau, E., Charles, M.,& Chen, M. (2014). External preference segmentation with additional information on consumers: A case study on apples. Food Quality and Preference, 32, 83-92. 
#' 
#' @export      
#' 
LCLV <-
function(X,Xr,Xu,ccX=FALSE,sX=TRUE,sXr=FALSE,sXu=FALSE,nmax=20)
{
 cX=TRUE
 cXr=TRUE
 cXu=FALSE
 # verification if some variables have constant values (standard deviation=0)
 who<-which(apply(X,2,sd)==0)
 if (length(who)>0) {
   listwho<-c(": ")
   for (r in 1:length(who)) {listwho=paste(listwho,colnames(X)[who[r]],",")}
   stop("The variables",listwho," have constant values (standard deviation=0). Please remove these variables from the X matrix.")
 }
 
 if (ccX==T) {
    X<- X-matrix(1,n,1)%*%colMeans(X)-rowMeans(X)%*%matrix(1,1,p)+mean(X)*matrix(1,n,p)
    cX=F
 }
 X<- scale(X, center=cX, scale=sX)
 p <- dim(X)[2] 
 n <- dim(X)[1] 
 if (is.null(Xu)) {stop("external variables Xu is missing")}                      
 EXTu<-1                   
 Xu<- scale(Xu, center=cXu, scale=sXu) 
 ptilde <- dim(Xu)[1]  
 m<-dim(Xu)[2] 
 if (p != ptilde)
 stop(" the number of consumers in X and Xu must are not the same") 

 M=diag(m)/(n-1) 
 if (is.null(Xr)) {stop("external variables Xr is missing")}                     
 EXTr<-1                                    
 Xr<- scale(Xr, center=cXr, scale=sXr)
 ntilde <- dim(Xr)[1]
 q<-dim(Xr)[2] 
 if (n != ntilde)
 stop("X and Xr must have the same number of observations")                              

 groupes <- 1:p                 
 fusions <- -groupes
 crit=c()
 for (i in 1:p){
   P= t(t(X[,i])) %*% Xu[i,] %*% M
   B=t(P)%*%Xr
   vp = eigen( B%*%t(B)) 
   critk= sqrt(vp$values[1])
   crit=c(crit,critk)
 }
 inertie <- sum(crit)
 sbegin <- sum(crit) 
 ncluster <- p        # nb clusters
 #print(paste('initial value of the criterion : ',round(inertie,2)))
 results = matrix(0,p-1,9)
 resultscc <- list()
 
 critstep <-  matrix(nrow=p, ncol=p)
 deltamin<-matrix(nrow=p,ncol=p)
 hmerge <-  matrix(0,nrow=p-1, ncol=2)
 delta <- matrix(0,nrow=p-1, ncol=1)
 ordr <- c()
#  matcov = t(X)%*%X
#  rescpp= critcpp(matcov,crit,dim(X)[1])
#  critstep = rescpp[[1]]
#  deltamin = rescpp[[2]]
 pg <- max(groupes) # nb groups
 for (i in 1:(pg - 1)) {
     for (j in (i + 1):pg) {
          P=as.matrix(X[,c(i,j)]) %*% as.matrix(Xu[c(i,j),]) %*%M
          B=t(P)%*%Xr
          vp = eigen(t(B)%*% B)
          critstep[i,j]= sqrt(vp$values[1])
          deltamin[i,j]<-crit[i]+crit[j]-critstep[i,j]
     }
 }  
                           
 for (level in 1:(p-1)) { 
 
# 		cmerge<-which(deltamin==min(deltamin,na.rm=TRUE),arr.ind=TRUE)
# 		if (nrow(cmerge)>1) cmerge<-cmerge[1,]  
		cmerge<-mincpp(deltamin)
		cmerge1 = cmerge[[1]]
		cmerge2 =  cmerge[[2]]
		ind1 <- which(groupes == cmerge1)
    ind2 <- which(groupes == cmerge2)
    ind<-c(fusions[ind1[1]],fusions[ind2[1]])
    hmerge[level,]<-ind
    listind<-c(ind1,ind2)  
    ordr <- order_var(ordr,listind)
    fusions[listind]<- level
    crit[listind[1]]<- critstep[cmerge1,cmerge2];
    autre <- listind[2:length(listind)]
    crit[autre]<-NA
   
    inertie <- rbind(inertie, (inertie[level] - deltamin[cmerge1,cmerge2]))    
	  delta[level]<-deltamin[cmerge1,cmerge2]
	  inertie_inv<- -inertie+sbegin
	  ncluster <- ncluster-1		                        
    results[level,1:6]<- c(ind,level,delta[level],inertie[level+1],100*inertie[level+1]/sbegin)     
    groupes[groupes==cmerge2]<-cmerge1
    groupes[groupes>cmerge2]<-groupes[groupes>cmerge2]-1
    deltamin<-deltamin[-cmerge2,-cmerge2]
    critstep<-critstep[-cmerge2,-cmerge2] 
    gr2<-which(groupes==cmerge1)
    if (cmerge1>1)
		{
			for (iter in 1:(cmerge1-1)){
				gr1<-which(groupes==iter) 
				P=as.matrix(X[,c(gr1,gr2)])%*%as.matrix(Xu[c(gr1,gr2),])%*%M
				B=t(P)%*%Xr
				vp = eigen(B %*% t(B))
				critstep[iter,cmerge1]= sqrt(vp$values[1])
        deltamin[iter,cmerge1]<-crit[gr1[1]]+crit[gr2[1]]-critstep[iter,cmerge1]	    
     }
    }
    if (cmerge1<max(groupes))
		 {
			for (iter in (cmerge1+1):max(groupes))
			{             
				gr1<-which(groupes==iter)
				P=as.matrix(X[,c(gr2, gr1)])%*%as.matrix(Xu[c(gr2, gr1),])%*%M
				B=t(P)%*%Xr
				vp = eigen(B %*% t(B))
				critstep[cmerge1,iter]= sqrt(vp$values[1])   
		    deltamin[cmerge1,iter]<-crit[gr2[1]]+crit[gr1[1]]-critstep[cmerge1,iter]		    
			}
    }
    if ((ncluster <= nmax) & (ncluster > 1) ) {
    cc_consol <- t(t(groupes))
    K <- ncluster
    T<-c()   
    maxiter=20
           
    for (i in 1:maxiter) {           
        critere <-rep(0,K)
        groupes_tmp <- cc_consol[,i]
        
        a <-matrix(0,nrow=q,ncol=K)
        rownames(a)<-colnames(Xr)
        colnames(a)<-paste("Comp",c(1:K),sep="")
        u<-matrix(0,nrow=m,ncol=K)
        rownames(u)<-colnames(Xu)
        colnames(u)<-paste("Comp",c(1:K),sep="")
        compc <-matrix(0,nrow=n,ncol=K)
        rownames(compc)<-rownames(X)
        colnames(compc)<-paste("Comp",c(1:K),sep="")
        compt <-matrix(0,nrow=n,ncol=K)
        rownames(compt)<-rownames(X)
        colnames(compt)<-paste("Comp",c(1:K),sep="")
        
        for (k in 1:K) {
            ind <- which(groupes_tmp == k)
            if (length(ind) > 0) {
              xgroupe<-X[,ind]
              xgroupe=as.matrix(xgroupe) 
              Xugroupe<-Xu[ind,]
              P=xgroupe%*%Xugroupe%*%M 
              B=t(P)%*%Xr    
              svd=svd(B)           
              critere[k]=svd$d[1]                            
              a[,k]=svd$v[,1]                
              compt[,k]=Xr%*%a[,k]
              u[,k]=svd$u[,1]
              compc[,k]=P%*%u[,k]           
            }
        }  
        T = cbind(T, sum(critere))
        
        for (j in 1:p) {
            critind=diag(t(u)%*%M%*%t(t(Xu[j,]))%*%t(X[,j])%*%compt)
            maxj=which(critind==max(critind))
            groupes_tmp[j] = maxj     
        }
 
      if (length(which((cc_consol[, i] == groupes_tmp) == FALSE, arr.ind = T)) == 0)    break
      cc_consol = cbind(cc_consol, groupes_tmp)
    }
    rownames(cc_consol) <- colnames(X)      
    names(cc_consol) = NULL
 
    initgroupes<-cc_consol[,1]
    lastgroupes<-cc_consol[,ncol(cc_consol)]
  
    listcc = list(clusters = rbind(initgroupes,lastgroupes),  compt=compt,compc=compc, loading_u=u,loading_v=a)
    results[level,7:9]<- c(sum(critere),(sum(critere)/sbegin)*100, i)  
    resultscc[[K]] <- listcc 
              	 	                                   
    }  
 }     

  results[p-1,7:9]<- c(results[p-1,5],results[p-1,6], 0)  
  group1<-matrix(1,nrow=2,ncol=p)
  colnames(group1) <- colnames(X)
  rownames(group1) <- c("initgroupes","lastgroupes")
 
  compc <-matrix(0,nrow=n,ncol=1)
  rownames(compc)<-rownames(X)
  colnames(compc)<-"Comp1"
  compt <-matrix(0,nrow=n,ncol=1)
  rownames(compt)<-rownames(X)
  colnames(compt)<-"Comp1"
        
  P=X%*%Xu%*%M
  B=t(P)%*%Xr
  svd = svd(B)
  a=svd$v[,1];                
  compt[,1]=Xr%*%a
  u=svd$u[,1]
  compc[,1]=P%*%u        
 
  listcc = list(clusters = group1,  compt=compt,compc=compc, loading_u=u,loading_v=a)
  resultscc[[1]] <- listcc 

  if (nmax == p) {
        groupes<-matrix(1:p,nrow=p,ncol=1)
        rownames(groupes) <- colnames(X)
        critere <-rep(0,p)
        compc <-matrix(0,nrow=n,ncol=p)
        rownames(compc)<-rownames(X)
        colnames(compc)<-paste("Comp",c(1:p),sep="")
        compt <-matrix(0,nrow=n,ncol=p)
        rownames(compt)<-rownames(X)
        colnames(compt)<-paste("Comp",c(1:p),sep="")
        a <-matrix(0,nrow=q,ncol=p) 
        u <-matrix(0,nrow=m,ncol=p)         
        for (k in 1:p) {
             ind <- which(groupes == k)
             xgroupe<-X[,ind]
             xgroupe=as.matrix(xgroupe) 
             Xugroupe<-Xu[ind,]
             P=xgroupe%*%Xugroupe%*%M 
             B=t(P)%*%Xr    
             svd=svd(B)           
             critere[k]=svd$d[1]
             a[,k]=svd$v[,1]                
             compt[,k]=Xr%*%a[,k]
             u[,k]=svd$u[,1]
             compc[,k]=P%*%u[,k]        
        }    
        groupes = cbind(groupes,groupes)  
        colnames(groupes) <- c("initgroupes","lastgroupes")
        listcc = list(clusters = t(groupes),compt=compt,compc=compc, loading_u=u,loading_v=a)
        resultscc[[p]] <- listcc      
  }  
 
 colnames(results)= c("merg1","merg2","new.clust","agg.crit.hac","clust.crit.hac",
                      "%S0expl.hac","clust.crit.cc","%S0expl.cc","iter")
 names(resultscc) = paste("partition",1:nmax,sep="")
 
 resultscc$tabres=results
 resultscc$param<-list(X=X,n = n, p = p,nmax = nmax,ccX=ccX,sX=sX,sXr=sXr,cXu=cXu,sXu=sXu,strategy="none")
 
 resultscah=list(labels=colnames(X),inertie=inertie, height=delta, merge=hmerge,order=ordr )    
 mytot<-resultscah  
 class(mytot)="hclust"
 mydendC=as.dendrogram(mytot)
  
 clvclt= c(resultscc, list(mydendC = mydendC))
 

 class(clvclt) = "lclv"
 return(clvclt) 
}
