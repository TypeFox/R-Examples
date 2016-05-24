# K-means of histogram data based on Wasserstein ----
#' K-means of a dataset of histogram-valued data
#' @description The function implements the k-means for a set of histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param k An integer, the number of groups.
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @return a list with the results of the k-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$centers A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references Irpino A., Verde R., Lechevallier Y. (2006). Dynamic clustering of histograms using Wasserstein 
#' metric. In: Rizzi A., Vichi M.. COMPSTAT 2006 - Advances in computational statistics. p. 869-876, 
#' Heidelberg:Physica-Verlag
#' @examples
#' results=WH_kmeans(x = BLOOD,k = 2, rep = 10,simplify = TRUE,qua = 10,standardize = TRUE)
#' @export
WH_kmeans =function (x,k, rep=5, 
                     simplify=FALSE,
                     qua=10,
                     standardize=FALSE){
  if ((k<2) || (k>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= k<",
              nrow(x@M), "(no. of individuals in x)\n")
    stop(str)
  }
  
  init="RPART"
  ind=nrow(x@M)
  vars=ncol(x@M)
  
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ## compute total sum of 
  TOTSSQ=0
  for (v in 1:vars){
    tmp=ComputeFastSSQ(MM[[v]])
    TOTSSQ=TOTSSQ+tmp$SSQ
    
  }
  solutions=list()
  criteria=numeric(0)
  for (repet in 1:rep){
    cat(paste("---------> rep  ",repet,"\n"))
    ## initialize clusters and prototypes
    proto=new("MatH",nrows=k,ncols=vars)
    cards=rep(0,k)
    SSQ=matrix(0,k,vars)
    GenCrit=Inf;
    if (init=="RPART"){
      rperm=sample(1:ind)
      tmp=rep(c(1:k),ceiling(ind/k))[1:ind]
      IDX=tmp[rperm]
      ##  compute prototypes
      for (clu in 1:k){
        sel=(1:ind)[IDX==clu]
        cards=length(sel)
        for (j in 1:vars){
          tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          proto@M[clu,j][[1]]=tmpH
          SSQ[clu,j]=tmp$SSQ
        }
      }
    }
    GenCrit=sum(SSQ)
    OK=1 
    maxit=100
    treshold=GenCrit*1e-10
    itcount=0
    
    while (OK==1){
      itcount=itcount+1
      #assign
      Dmat=array(0,c(ind,vars,k))
      MD=matrix(0,ind,k)
      for (i in 1:ind){
        for (cl in 1:k){
          for (j in 1:vars){
            tmp=ComputeFast_L2_SQ_WASS_D(cbind(MM[[j]][,i],proto@M[cl,j][[1]]@x,MM[[j]][,(ind+1)]))
            Dmat[i,j,cl]=tmp#WassSqDistH(x@M[i,j][[1]],proto@M[cl,j][[1]])
            MD[i,cl]=MD[i,cl]+Dmat[i,j,cl]
          }
        }
        IDX[i]=which.min(MD[i,])
        
      }
     
      #recompute
      OldCrit=GenCrit
      
      #check for empty clusters
      moved=numeric(0)
      for (i in 1:k){
        sel=(1:ind)[IDX==i]
        #show(sel)
        if (length(sel)==0){
          cat("empty cluster\n")
          which.max(apply(MD,1,max))
          tomove=which.max(apply(MD,1,max))
          IDX[tomove]=i
          moved=c(moved,tomove)
          MD[tomove,]=rep(0,k)
          #show(MD)
        }
      }
      for (i in 1:k){
        sel=(1:ind)[IDX==i]
        
        cards=length(sel)
        for (j in 1:vars){
          tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          proto@M[i,j][[1]]=tmpH
          SSQ[i,j]=tmp$SSQ
        }
      }
      GenCrit=sum(SSQ)
      cat(paste(itcount,GenCrit, "\n", sep="---->"))
      #check criterion
      if (abs(GenCrit-OldCrit)<treshold){OK=0}
    }
    cardinality=table(IDX)
    dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
    solutions=c(solutions,list(solution=list(IDX=IDX,cardinality=cardinality,centers=proto,Crit=GenCrit)))
    criteria=c(criteria,GenCrit)
    plot(proto,type="DENS")
  }
  return(best.solution=list(solution=solutions[[which.min(criteria)]],
                            quality=1-min(criteria)/TOTSSQ))
}


# Hierarchical clustering of histogram data based on Wasserstein ----
#' Hierarchical clustering of histogram data
#' @description The function implements a Hierarchical clustering 
#'  for a set of histogram-valued data, based on the L2 Wassertein distance.
#'  Extends the \code{hclust} function of the \pkg{stat} package. 
#' @param x A MatH object (a matrix of distributionH).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @param distance A string default "WDIST" the L2 Wasserstein distance (other distances will be implemented)
#' @param method A string, default="complete", is the the agglomeration method to be used.
#'  This should be (an unambiguous abbreviation of) one of "\code{ward.D}", "\code{ward.D2}",
#'   "\code{single}", "\code{complete}", "\code{average}" (= UPGMA), "\code{mcquitty}" 
#'   (= WPGMA), "\code{median}" (= WPGMC) or "\code{centroid}" (= UPGMC).
#' @seealso \code{\link{hclust}} of \pkg{stat} package for further details.
#' @return An object of class hclust which describes the tree produced by
#'  the clustering process. 
#' @references Irpino A., Verde R. (2006). A new Wasserstein based distance for the hierarchical clustering 
#' of histogram symbolic data. In: Batanjeli et al. Data Science and Classification, IFCS 2006. p. 185-192,
#'  BERLIN:Springer, ISBN: 3-540-34415-2
#' @examples
#' results=WH_hclust(x = BLOOD,simplify = TRUE, method="complete")
#' plot(results) # it plots the dendrogram
#' cutree(results,k = 5) # it returns the labels for 5 clusters
#' @importFrom stats hclust quantile as.dist
#' @export
WH_hclust =function (x, 
                     simplify=FALSE,
                     qua=10,
                     standardize=FALSE,
                     distance="WDIST",
                     method="complete"){
  ind=nrow(x@M)
  vars=ncol(x@M)
  
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ##compute distance matrix
  d=matrix(0,ind,ind)
  for (i in 1:(ind-1)){
    for (j in (i+1):ind){
      for (v in 1:vars){
        tmp=ComputeFast_L2_SQ_WASS_D(cbind(MM[[v]][,i],MM[[v]][,j],MM[[v]][,(ind+1)]))
        d[i,j]= d[i,j]+tmp#WassSqDistH(x@M[i,v][[1]],x@M[j,v][[1]])
      }
      if (distance=="WDIST"){d[i,j]=sqrt(d[i,j])}
      d[j,i]=d[i,j]
    }
  }
  rownames(d)=rownames(x@M)
  colnames(d)=rownames(x@M)
  hc=hclust(as.dist(d),method=method)
  return(hc)
}


# Adaptive-distance based k-means ----
#' K-means of a dataset of histogram-valued data using adaptive  Wasserstein distances
#' @description The function implements the k-means using adaptive distance for a set of histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param k An integer, the number of groups.
#' @param schema a number from 1 to 4 \cr
#' 1=A weight for each variable (default) \cr
#' 2=A weight for the average and the dispersion component of each variable\cr
#' 3=Same as 1 but a different set of weights for each cluster\cr
#' 4=Same as 2 but a different set of weights for each cluster 
#' @param init (optional, do not use) initialization for partitioning the data default is 'RPART', other strategies shoul be implemented.
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @param weight.sys a string. Weights may add to one ('SUM') or their product is equal to 1 ('PROD', default).
#' @param theta a number. A parameter if \code{weight.sys='SUM'}, default is 2.  
#' @param init.weights a string how to initialize weights: 'EQUAL' (default), all weights are the same, 
#' 'RANDOM', weights are initalised at random.
#' 
#' @return a list with the results of the k-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$centers A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references Irpino A., Rosanna V., De Carvalho F.A.T. (2014). Dynamic clustering of histogram data 
#' based on adaptive squared Wasserstein distances. EXPERT SYSTEMS WITH APPLICATIONS, vol. 41, p. 3351-3366, 
#' ISSN: 0957-4174, doi: http://dx.doi.org/10.1016/j.eswa.2013.12.001
#' @examples
#' results=WH_adaptive.kmeans(x = BLOOD,k = 2, rep = 10,simplify = TRUE,qua = 10,standardize = TRUE)
#' @importFrom stats runif
#' @export
WH_adaptive.kmeans =function (x,k,
                              schema=1, #1=VariableGLOBAL 2=componentGLOBAL 3=Variable x Cluster 4=Components x cluster 
                              init, rep,  
                              simplify=FALSE,
                              qua=10,
                              #empty.action="singleton",
                              standardize=FALSE,
                              weight.sys='PROD',
                              theta=2,
                              init.weights='EQUAL'){
  if ((k<2) || (k>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= k<",
              nrow(x@M), "(no. of individuals in x)\n")
    stop(str)
  }
  if (missing(init)){
    init="RPART"
  }
  if (missing(rep)){
    rep=5
  }
  ind=nrow(x@M)
  vars=ncol(x@M)
  
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ## compute total sum of 
  TOTSSQ=0
  for (v in 1:vars){
    tmp=ComputeFastSSQ(MM[[v]])
    TOTSSQ=TOTSSQ+tmp$SSQ
    }
  
  solutions=list()
  criteria=numeric(0)
  repet=0
  while(repet<rep){
    repet=repet+1
    cat(paste("---------> rep  ",repet,"\n"))
    #initialize matrix of variables' weights
    if (init.weights=='EQUAL'){
      cat("Weights initialization === EQUAL  \n")
      if (weight.sys=='PROD'){
        lambdas=matrix(1,2*vars,k)}
      else{lambdas=matrix(1/vars,2*vars,k)}
    }
    else{#Random initialization
      cat("Weights initialization === RANDOM  \n")
      m1=matrix(runif((vars*k),0.01,0.99),vars,k)
      m2=matrix(runif((vars*k),0.01,0.99),vars,k)
      m1=m1/matrix(rep(apply(m1,2,sum)),vars,k,byrow = TRUE)
      m2=m2/matrix(rep(apply(m2,2,sum)),vars,k,byrow = TRUE)
      if (weight.sys=='PROD'){
        m1=exp(m1*vars-1)
        m2=exp(m2*vars-1)
        
      }
      
      if (schema==1){
        m1=matrix(m1[,1],nrow = vars,ncol = k)
        m2=m1
        
      }
      if (schema==2){
        m1=matrix(rep(m1[,1],k),vars,k)
        m2=matrix(rep(m2[,1],k),vars,k)
      }
      if (schema==3){m2=m1}
      
      lambdas=matrix(0,2*vars,k)
      colnames(lambdas)=paste('Clust',c(1:k),sep="_")
      
      n1=paste('M_Var.',c(1:vars),sep="_")
      n2=paste('C_Var.',c(1:vars),sep="_")
      nr=list()
      for (nn in 1:vars){
        nr[[nn*2-1]]=n1[[nn]]
        nr[[nn*2]]=n2[[nn]]
      }
      rownames(lambdas)=nr
      lambdas[(c(1:vars)*2-1),]=m1
      lambdas[(c(1:vars)*2),]=m2
    }
    ## initialize clusters and prototypes
    proto=new("MatH",nrows=k,ncols=vars)
    cards=rep(0,k)
    SSQ=matrix(0,k,vars)
    GenCrit=Inf;
    if (init=="RPART"){
      rperm=sample(1:ind)
      tmp=rep(c(1:k),ceiling(ind/k))[1:ind]
      IDX=tmp[rperm]
      
    }
    memb=matrix(0,ind,k)
    for(indiv in 1:ind){memb[indiv,IDX[indiv]]=1}
    GenCrit=Inf
    OK=1 
    maxit=100
    treshold=1e-5
    itcount=0
    
    while (OK==1){
      itcount=itcount+1
      #STEP 1 ##  compute prototypes
      for (cluster in 1:k){
        sel=(1:ind)[IDX==cluster]
        #show(sel)
        cards=length(sel)
        for (variables in 1:vars){
          tmp=ComputeFastSSQ(MM[[variables]][,c(sel,(ind+1))])
          tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          proto@M[cluster,variables][[1]]=tmpH
          SSQ[cluster,variables]=tmp$SSQ
         
        }
      }
      #       TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto)
      #       cat('\n after prototypes SSQ:', TMP_SSQ,'\n')
      #STEP 2 ##  compute weights (Lambda) Fixed Prototypes and partition
      #######################################################################################
      distances=array(0,dim=c(vars,k,2))
      diINDtoPROT=array(0,dim=c(ind,vars,k,2))
      for (cluster in 1:k){
        for (variables in (1:vars)){
          for (indiv in 1:ind){
            tmpD=ComputeFast_L2_SQ_WASS_D(cbind(MM[[variables]][,indiv],proto@M[cluster,variables][[1]]@x,
                                                MM[[variables]][,(ind+1)]))
            #tmpD1=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
            
            #if (memb[indiv,cluster]>0){
            if (schema==1){#one weigth for one variable
              #tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*memb[indiv,cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD
            }
            if (schema==2){#two weigths for the two components for each variable
              #tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
              tmpD_mean=(x@M[indiv,variables][[1]]@m-proto@M[cluster,variables][[1]]@m)^2
              tmpD_centered=tmpD-tmpD_mean
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*memb[indiv,cluster]
              distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*memb[indiv,cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
              diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
            }
            if (schema==3){#a weigth for one variable and for each cluster
              #tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*memb[indiv,cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD
            }
            if (schema==4){#two weigths for the two components for each variable and each cluster
              #tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
              tmpD_mean=(x@M[indiv,variables][[1]]@m-proto@M[cluster,variables][[1]]@m)^2
              tmpD_centered=tmpD-tmpD_mean
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*memb[indiv,cluster]
              distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*memb[indiv,cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
              diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
            }
          }
          #}
        }
      }
      ##  compute weights
      # S2.2) Weights computation
      
      for (variables in (1:vars)){
        for (cluster in 1:k){
          #product
          if (weight.sys=='PROD'){
            if (schema==1){#one weigth for one variable
              if (cluster<=1){
                num=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/vars)
                denom=max(sum(distances[variables,,1]),1e-10)
                lambdas[(variables*2-1),cluster]=num/denom
                lambdas[(variables*2),cluster]=num/denom
              }else{lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
                    lambdas[(variables*2),cluster]=lambdas[(variables*2),1]}
            }
            if (schema==2){#two weigths for the two components for each variable
              if (cluster<=1){
                num=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/vars)
                denom=max(sum(distances[variables,,1]),1e-10)
                lambdas[(variables*2-1),cluster]=num/denom
                num=(prod(apply(distances[,,2],MARGIN=c(1),sum)))^(1/vars)
                denom=max(sum(distances[variables,,2]),1e-10)
                lambdas[(variables*2),cluster]=num/denom}
              else{
                lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
                lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
              }
            }
            if (schema==3){#a weigth for one variable and for each cluster
              num=(prod(distances[,cluster,1]))^(1/vars)
              denom=max(distances[variables,cluster,1],1e-10)
              lambdas[(variables*2-1),cluster]=num/denom
              lambdas[(variables*2),cluster]=num/denom
            }
            if (schema==4){#two weigths for the two components for each variable and each cluster
              num=max((prod(distances[,cluster,1]))^(1/vars),1e-10)
              denom=max(distances[variables,cluster,1],1e-10)
              
              lambdas[(variables*2-1),cluster]=num/denom
              num=max((prod(distances[,cluster,2]))^(1/vars),1e-10)
              denom=max(distances[variables,cluster,2],1e-10)
              lambdas[(variables*2),cluster]=num/denom
            }
          }else{
            #sum ugual to 1
            if (schema==1){#one weigth for one variable 1W GS Global-sum
              if (cluster<=1){
                num=rep(sum(distances[variables,,1]),vars)
                den=apply(distances[,,1],1,sum)
                lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
                lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
              }else{lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
                    lambdas[(variables*2),cluster]=lambdas[(variables*2),1]}
            }
            if (schema==2){#two weigths for the two components for each variable 2W GS Global-sum
              if (cluster<=1){
                num=rep(sum(distances[variables,,1]),vars)
                den=apply(distances[,,1],1,sum)
                lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
                num=rep(sum(distances[variables,,2]),vars)
                den=apply(distances[,,2],1,sum)
                lambdas[(variables*2),cluster]=1/sum((num/den)^(1/(theta-1)))
              }
              else{
                lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
                lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
              }
            }
            if (schema==3){#a weigth for one variable and for each cluster 1W LS Local-sum
              num=rep(distances[variables,cluster,1],vars)
              den=distances[,cluster,1]
              
              lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
              lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
            }
            if (schema==4){#two weigths for the two components for each variable and each cluster 2W LS Local-sum
              num=rep(distances[variables,cluster,1],vars)
              den=distances[,cluster,1]
              
              lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
              
              num=rep(distances[variables,cluster,2],vars)
              den=distances[,cluster,2]
              lambdas[(variables*2),cluster]=1/sum((num/den)^(1/(theta-1)))
              
            }
          }
        }
      }
      ############################################################################
      #       TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto)
      #       cat('\n after weights SSQ:', TMP_SSQ,'\n')      
      #recompute
      OldCrit=GenCrit
      # S3) affectation (prototype, weights are fixed)
      DiToClu=matrix(0,ind,k)
      for (indiv in 1:ind){
        for (cluster in 1:k){
          for (variables in 1:vars){
            DiToClu[indiv,cluster]=DiToClu[indiv,cluster]+
              lambdas[(variables*2-1),cluster]*diINDtoPROT[indiv,variables,cluster,1]+
              lambdas[(variables*2),cluster]*diINDtoPROT[indiv,variables,cluster,2]
          }
        }
        IDX[indiv]=which.min(DiToClu[indiv,])
      }
      
      #check for empty clusters and restart
      empty=FALSE;
      for (i in 1:k){
        sel=(1:ind)[IDX==i]
        #show(sel)
        if (length(sel)==0){
          
          cat("empty cluster\n")
          empty=TRUE
          OK=0
        }
      }
      if (empty){repet=repet-1}
      memb=matrix(0,ind,k)
      for(indiv in 1:ind){memb[indiv,IDX[indiv]]=1}
      TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto)
      GenCrit=TMP_SSQ
      if(is.na(GenCrit)){
        cat('isNAN')
      }
      
      cat(paste(itcount,GenCrit, "\n", sep="---->"))
      #check criterion
      if (abs(GenCrit-OldCrit)<treshold){OK=0}
    }
    if(!empty){  
      cardinality=table(IDX)
      dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
      solutions=c(solutions,list(solution=list(IDX=IDX,cardinality=cardinality,proto=proto, weights=lambdas,Crit=GenCrit)))
      criteria=c(criteria,GenCrit)
    }
  }
  best.solution=c(solutions[[which.min(criteria)]])
  plot(best.solution$proto,type="DENS")
   memb=matrix(0,ind,k)
   for (i in 1:ind){
     memb[i,best.solution$IDX[i]]=1;
   }
   resTOTSQ=WH.ADPT.KMEANS.TOTALSSQ(x,memb,m=1,
                                     lambdas=best.solution$weights,proto=best.solution$proto)
   GEN_proto=resTOTSQ$protogen
   DET_TSQ=resTOTSQ$TSQ
   TOTFSSQ=sum(resTOTSQ$TSQ)
  DET_TSQ2=resTOTSQ$TSQ2
  TOTFSSQ2=sum(resTOTSQ$TSQ2)
  BSQ=resTOTSQ$BSQ
   quality1=1-min(criteria)/TOTFSSQ
  quality2=1-min(criteria)/TOTFSSQ2
  best.solution=c(best.solution,TOTSSQ=TOTFSSQ2,BSQ=BSQ,WSQ=best.solution$Crit,quality=BSQ/TOTFSSQ2)
  return(best.solution)
#  return(0)
}


