# Batch SOM Kohonen Maps of HD using adaptive distances -----
#' Batch Kohonen self-organizing 2d maps using adaptive distances for  histogram-valued data
#' @description The function implements a Batch Kohonen self-organizing 2d maps algorithm for  histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param net a list describing the topology of the net \code{list(xdim=number of rows,
#' ydim=numbers of columns,topo=c('rectangular' or 'hexagonal'))}, see \code{somgrid} sintax in package\pkg{class}
#' default \code{net=list(xdim=4,ydim=3,topo=c('rectangular'))}
#' @param kern.param (default =2) the kernel parameter for the RBF kernel used in the algorithm
#' @param TMAX a parameter useful for the iterations (default=2)
#' @param Tmin a parameter useful for the iterations (default=0.2)
#' @param niter maximum number of iterations (default=30)
#' @param rep number of repetion of the algorithm (default=5), beacuase each launch may generate a local optimum
#' @param simplify a logical parameter for speeding up computations (default=FALSE). If true data are recoded in order to have fast computations
#' @param qua if \code{simplify=TRUE} number of equally spaced quantiles for recodify the histograms (default=10)
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one. 
#' @param schema a number from 1 to 4 \cr
#' 1=A weight for each variable (default) \cr
#' 2=A weight for the average and the dispersion component of each variable\cr
#' 3=Same as 1 but a different set of weights for each cluster\cr
#' 4=Same as 2 but a different set of weights for each cluster 
#' @param init.weights a string how to initialize weights: 'EQUAL' (default), all weights are the same, 
#' @param weight.sys a string. Weights may add to one ('SUM') or their product is equal to 1 ('PROD', default).
#' @param theta a number. A parameter if \code{weight.sys='SUM'}, default is 2.  
#' @param Wfix a logical parameter (default=FALSE). If TRUE the algorithm does not use adaptive distances
#'         
#' @return a list with the results of the Batch Kohonen map
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$MAP The map topology.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$proto A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot solution$Weights.comp the final weights assigned to each component of the histogram variables
#' @slot solution$Weight.sys a string the type of weighting system ('SUM' or 'PRODUCT')
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references 
#' Irpino A, Verde R, De Carvalho FAT (2012). Batch self organizing maps for interval and histogram data.
#' In: Proceedings of COMPSTAT 2012. p. 143-154, ISI/IASC, ISBN: 978-90-73592-32-2
#' @details An extension of Batch Self Organised Map (BSOM) is here proposed for  histogram data.
#'  These kind of data have been defined in the context of symbolic data analysis.
#'   The BSOM cost function is then based on a distance 
#'   function: the L2 Wasserstein distance. This distance has been widely proposed in several
#'    techniques of analysis (clustering, regression) when input data are expressed by distributions 
#'    (empirical by histograms or theoretical by probability distributions).
#'    The peculiarity of such distance is to be an Euclidean distance between quantile functions so 
#'    that all the properties proved for L2 distances are verified again. An adaptative versions of 
#'    BSOM is also introduced considering an automatic system of weights in the cost function in 
#'    order to take into account the different effect of the several variables in the Self-Organised Map 
#'    grid. 
#' 
#' @importFrom class somgrid
#' @importFrom stats lm dist runif
#' @examples
#' \dontrun{
#' results=WH_2d_Adaptive_Kohonen_maps(x = BLOOD,k = 2,
#'                                    net=list(xdim=2,ydim=3,topo=c('rectangular')), 
#'                                    rep = 2,simplify = TRUE,
#'                                    qua = 10,standardize = TRUE)
#'                                    }
#' @export
WH_2d_Adaptive_Kohonen_maps =function (x,net=list(xdim=4,ydim=3,topo=c('rectangular')), kern.param=2, TMAX=2, Tmin=0.2, 
                              niter=30,rep ,
                      simplify=FALSE,
                      qua=10,
                      standardize=FALSE, schema=4,
                      init.weights='EQUAL',weight.sys='PROD',theta=2,Wfix=FALSE){
  tol=1e-10
  #require(class)#for somgrid function
  # Check initial conditions and passed parameters
  
  ind=nrow(x@M)
  vars=ncol(x@M)
  ## we homogeneize data for speeding up the code if required
  if (simplify){
    p=(0:qua)/qua
    for (j in 1:vars){
      for (i in 1:ind){
        dom=numeric(0)
        for (q in 0:qua){
          dom=c(dom, compQ(x@M[i,j][[1]],q/qua))
        }
        tmp=new("distributionH",dom,p)
        tmp=(tmp-tmp@m)*(x@M[i,j][[1]]@s/tmp@s)+x@M[i,j][[1]]@m #transformation with invariance with respect mean and std
        x@M[i,j][[1]]=tmp
      }
    } 
  }
  else{
    for (j in 1:vars){
      tmp=registerMH(x[,j])
      for (i in 1:ind){
        x@M[i,j][[1]]=tmp@M[i,1][[1]]
      }
    }
  }
  ## standardize data if required
  if (standardize){
    cat("Standardizing data...\n")
    STAND=rep(0,vars)
    Mc=rep(0,vars)
    # compute varianaces
    for (v in 1:vars){
      STAND[v]=sqrt(WH.var.covar(x[,v]))
      Mc[v]=(WH.vec.mean(x[,v]))@m
      for (i in 1:ind){
        if (STAND[v]>0){
          x@M[i,v][[1]]=new("distributionH",x=(x@M[i,v][[1]]@x-Mc[v])/STAND[v],p=x@M[i,v][[1]]@p)
        }
      }
    }
  }
  # end of the preprocessing step
  
  TOTSSQ=0
  for (v in 1:vars){
    TOTSSQ=TOTSSQ+(WH.SSQ(x[,v]))
  }
  ##------------------------------------------------batchKOHONENMAP-----------------------------------------------------
  solutions=list()
  criteria=numeric(0)
  MAP=somgrid(net$xdim,net$ydim,net$topo)
  k=nrow(MAP$pts)
  if (missing(rep)){
    rep=5
  }
  for (repet in (1:rep)){
    data=x@M
    nd=nrow(data)
    #random selection of prototypes of neurons
    init=data[sample(1L:nd,k, replace=FALSE), ,drop=FALSE]
    proto=new("MatH")
    proto@M=init
    nhbrdist=as.matrix(dist(MAP$pts))
    TT=TMAX
    KT=exp(-(nhbrdist^2)/(2*TT^2))
    #initialize weights
    if (Wfix){lambdas=matrix(1,2*vars,k)}else{
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
    }}
    #assign objects to closest neuron
    #compute adaptive distances to prototypes
    distances=array(0,dim=c(vars,k,2))
    diINDtoPROT=array(0,dim=c(ind,vars,k,2))
    
    for (cluster in 1:k){
      for (variables in (1:vars)){
        for (indiv in 1:ind){
          if (schema==1){#one weight for one variable
            tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
            distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*lambdas[(variables*2-1),cluster]
            diINDtoPROT[indiv,variables,cluster,1]=tmpD
          }
          if (schema==2){#two weights for the two components for each variable
            tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
            tmpD_mean=as.numeric(tmpD[2])
            tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
            distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*lambdas[(variables*2-1),cluster]
            distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*lambdas[(variables*2),cluster]
            diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
            diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
          }
          if (schema==3){#a weight for one variable and for each cluster
            tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
            distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*lambdas[(variables*2-1),cluster]
            diINDtoPROT[indiv,variables,cluster,1]=tmpD
          }
          if (schema==4){#two weights for the two components for each variable and each cluster
            tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
            tmpD_mean=as.numeric(tmpD[2])
            tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
            distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*lambdas[(variables*2-1),cluster]
            distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*lambdas[(variables*2),cluster]
            diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
            diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
          }
        }
      }
    }
    Dind2Neu=matrix(Inf,ind,k)
    for (indiv in 1:ind){
      for (cluster in 1:k){#s
        if (Dind2Neu[indiv,cluster]==Inf) Dind2Neu[indiv,cluster]=0
        for (otherclust in 1:k){#l
        Dind2Neu[indiv,cluster]=Dind2Neu[indiv,cluster]+KT[otherclust,cluster]*sum(diINDtoPROT[indiv,,otherclust,])
        }
      }
    }
    #initialize matrix of memberships
    IDX=apply(Dind2Neu,1,which.min)
    memb=matrix(0,ind,k)
    for (i in (1:ind)){
      memb[i,IDX[i]]=1
    }
    cat(paste("---------> rep  ",repet,"\n"))
    ## initialize clusters and prototypes
    proto=new("MatH",list(new("distributionH")),nrows=k,ncols=vars)
    SSQ1=matrix(0,k,vars)
    GenCrit=Inf
    ##  compute initial criterion
    
    SSQ=0
    for (cluster in 1:k){
      for (variables in (1:vars)){
        for (indiv in 1:ind){
          if (memb[indiv,cluster]>0)
          SSQ=SSQ+Dind2Neu[indiv,IDX[indiv]]
        }
      }
    }
    GenCrit=SSQ
    OK=1 
    maxit=100

   t=0
    while (TT>=Tmin){
      t=t+1
      TT=TMAX*(Tmin/TMAX)^(t/(niter-1))
      KT=exp(-(nhbrdist^2)/(2*TT^2))
      #computing prototypes
      CardinalityOfNeuron=apply(memb,2,sum)
      
      for (cluster in 1:k){
        kern=rep(0,ind)
        for (individuals in 1:ind){
          kern[individuals]=KT[IDX[individuals],cluster]
        }
        for (variables in 1:vars){
          proto@M[cluster,variables][[1]]=WH.vec.mean(x[,variables],kern)
        }
      }
      
      #compute weights
      if (Wfix==FALSE){
      #first compute distances using kernel
      distances=array(0,dim=c(vars,k,2))
      diINDtoPROT=array(0,dim=c(ind,vars,k,2))
      
      for (cluster in 1:k){
        #if (length(proto@M[cluster,1]@M)>0){# Check for empty clusters i.e. null prototypes
          for (variables in (1:vars)){
            for (indiv in 1:ind){
              if (schema==1){#one weight for one variable
                tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
                distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*KT[IDX[indiv],cluster]
                diINDtoPROT[indiv,variables,cluster,1]=tmpD
              }
              if (schema==2){#two weights for the two components for each variable
                tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
                tmpD_mean=as.numeric(tmpD[2])
                tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
                distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*KT[IDX[indiv],cluster]
                distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*KT[IDX[indiv],cluster]
                diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
                diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
              }
              if (schema==3){#a weight for one variable and for each cluster
                tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
                distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*KT[IDX[indiv],cluster]
                diINDtoPROT[indiv,variables,cluster,1]=tmpD
              }
              if (schema==4){#two weights for the two components for each variable and each cluster
                tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
                tmpD_mean=as.numeric(tmpD[2])
                tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
                distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*KT[IDX[indiv],cluster]
                distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*KT[IDX[indiv],cluster]
                diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
                diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
              }
            }
          }
       # }
      }
      #Weights computation
      
      for (variables in (1:vars)){
        for (cluster in 1:k){
          #product
          if (weight.sys=='PROD'){
            if (schema==1){#one weight for one variable
              if (cluster<=1){
                num=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/vars)
                denom=max(sum(distances[variables,,1]),1e-10)
                lambdas[(variables*2-1),cluster]=num/denom
                lambdas[(variables*2),cluster]=num/denom
              }else{lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
                    lambdas[(variables*2),cluster]=lambdas[(variables*2),1]}
            }
            if (schema==2){#two weights for the two components for each variable
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
            if (schema==3){#a weight for one variable and for each cluster
              num=(prod(distances[,cluster,1]))^(1/vars)
              denom=max(distances[variables,cluster,1],1e-10)
              lambdas[(variables*2-1),cluster]=num/denom
              lambdas[(variables*2),cluster]=num/denom
            }
            if (schema==4){#two weights for the two components for each variable and each cluster
              num=(prod(distances[,cluster,1]))^(1/vars)
              denom=max(distances[variables,cluster,1],1e-10)
              lambdas[(variables*2-1),cluster]=num/denom
              num=(prod(distances[,cluster,2]))^(1/vars)
              denom=max(distances[variables,cluster,2],1e-10)
              lambdas[(variables*2),cluster]=num/denom
            }
          }else{
            #sum ugual to 1
            if (schema==1){#one weight for one variable 1W GS Global-sum
              if (cluster<=1){
                num=rep(sum(distances[variables,,1]),vars)
                den=apply(distances[,,1],1,sum)
                lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
                lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
              }else{lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
                    lambdas[(variables*2),cluster]=lambdas[(variables*2),1]}
            }
            if (schema==2){#two weights for the two components for each variable 2W GS Global-sum
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
            if (schema==3){#a weight for one variable and for each cluster 1W LS Local-sum
              num=rep(distances[variables,cluster,1],vars)
              den=distances[,cluster,1]
              
              lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
              lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
            }
            if (schema==4){#two weights for the two components for each variable and each cluster 2W LS Local-sum
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
      }
      #assign data to neurons
      #assign objects to closest neuron
      #compute adaptive distances to prototypes
      distances=array(0,dim=c(vars,k,2))
      diINDtoPROT=array(0,dim=c(ind,vars,k,2))
      
      for (cluster in 1:k){
        for (variables in (1:vars)){
          for (indiv in 1:ind){
            if (schema==1){#one weight for one variable
              tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*lambdas[(variables*2-1),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD
            }
            if (schema==2){#two weights for the two components for each variable
              tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
              tmpD_mean=as.numeric(tmpD[2])
              tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*lambdas[(variables*2-1),cluster]
              distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*lambdas[(variables*2),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
              diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
            }
            if (schema==3){#a weight for one variable and for each cluster
              tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]])
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*lambdas[(variables*2-1),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD
            }
            if (schema==4){#two weights for the two components for each variable and each cluster
              tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
              tmpD_mean=as.numeric(tmpD[2])
              tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
              distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*lambdas[(variables*2-1),cluster]
              distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*lambdas[(variables*2),cluster]
              diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
              diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
            }
          }
        }
      }
      Dind2Neu=matrix(Inf,ind,k)
      for (indiv in 1:ind){
        for (cluster in 1:k){#s
          if (Dind2Neu[indiv,cluster]==Inf) Dind2Neu[indiv,cluster]=0
          for (otherclust in 1:k){#l
            Dind2Neu[indiv,cluster]=Dind2Neu[indiv,cluster]+KT[otherclust,cluster]*sum(diINDtoPROT[indiv,,otherclust,])
          }
        }
      }
      #recompute criterion
      IDX=apply(Dind2Neu,1,which.min)
      OldCrit=GenCrit
      GenCrit=0
      for (individual in 1:ind){
      GenCrit=GenCrit+sum(Dind2Neu[individual, IDX[individual]])
      }
      cat(paste(t,GenCrit, "\n", sep="---->"))
      
    }
    #crisp assignment
    
    cardinality=table(IDX)
    #dimnames(cardinality)$IDX=paste("Cl",1:k,sep=".")
    solutions=c(solutions,list(solution=list(MAP=MAP, IDX=IDX,cardinality=cardinality,proto=proto,
                                             Crit=GenCrit,weights.comp=lambdas,Weight.sys=weight.sys)))
    criteria=c(criteria,GenCrit)
    plot(proto,type="DENS")
  }
  return(best.solution=list(solution=solutions[[which.min(criteria)]],
                            rep=which.min(criteria),
                            quality=1-min(criteria)/TOTSSQ))
}

# Batch SOM of HD -----
#' Batch Kohonen self-organizing 2d maps for  histogram-valued data
#' @description The function implements a Batch Kohonen self-organizing 2d maps algorithm for  histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param net a list describing the topology of the net \code{list(xdim=number of rows,
#' ydim=numbers of columns,topo=c('rectangular' or 'hexagonal'))}, see \code{somgrid} sintax in package\pkg{class} 
#' default \code{net=list(xdim=4,ydim=3,topo=c('rectangular'))}
#' @param kern.param (default =2) the kernel parameter for the RBF kernel used in the algorithm
#' @param TMAX a parameter useful for the iterations (default=2)
#' @param Tmin a parameter useful for the iterations (default=0.2)
#' @param niter maximum number of iterations (default=30)
#' @param rep number of repetion of the algorithm (default=5), beacuase each launch may generate a local optimum
#' @param simplify a logical parameter for speeding up computations (default=FALSE). If true data are recoded in order to have fast computations
#' @param qua if \code{simplify=TRUE} number of equally spaced quantiles for recodify the histograms (default=10)
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @return a list with the results of the Batch Kohonen map
#'@slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$MAP The map topology.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$proto A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references 
#' Irpino A, Verde R, De Carvalho FAT (2012). Batch self organizing maps for interval and histogram data.
#' In: Proceedings of COMPSTAT 2012. p. 143-154, ISI/IASC, ISBN: 978-90-73592-32-2
#' @details An extension of Batch Self Organised Map (BSOM) is here proposed for  histogram data.
#'  These kind of data have been defined in the context of symbolic data analysis.
#'   The BSOM cost function is then based on a distance 
#'   function: the L2 Wasserstein distance. This distance has been widely proposed in several
#'    techniques of analysis (clustering, regression) when input data are expressed by distributions 
#'    (empirical by histograms or theoretical by probability distributions).
#'    The peculiarity of such distance is to be an Euclidean distance between quantile functions so 
#'    that all the properties proved for L2 distances are verified again. An adaptative versions of 
#'    BSOM is also introduced considering an automatic system of weights in the cost function in 
#'    order to take into account the different effect of the several variables in the Self-Organised Map 
#'    grid. 
#' 
#' @importFrom class somgrid
#' @examples
#' \dontrun{
#' results=WH_2d_Kohonen_maps(x = BLOOD,k = 2,
#'                                    net=list(xdim=2,ydim=3,topo=c('rectangular')), 
#'                                    rep = 2,simplify = TRUE,
#'                                    qua = 10,standardize = TRUE)
#'                                    }
#' @export
WH_2d_Kohonen_maps =function (x,net=list(xdim=4,ydim=3,topo=c('rectangular')),
                              kern.param=2, TMAX=2, Tmin=0.2, 
                              niter=30,rep=5 ,      
                              simplify=FALSE,
                                    qua=10,
                                    standardize=FALSE){
  SOL=WH_2d_Adaptive_Kohonen_maps(x,net, kern.param, TMAX, Tmin, 
                                         niter,rep ,
                                         simplify,
                                         qua,
                                         standardize, schema=1,
                                         init.weights='EQUAL',weight.sys='PROD',theta=2,Wfix=TRUE)
  
  return(SOL)
}

WH.ADPT.KOHONEN.SSQ=function(x,memb,m,lambdas,proto){
  vars=ncol(x@M)
  ind=nrow(x@M)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  SSQ=0
  for (indiv in 1:ind){
  for (cluster in 1:k){
    for (variables in (1:vars)){
      tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
      tmpD_mean=as.numeric(tmpD[2])
      tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
      SSQ=SSQ+((memb[indiv,cluster])^m)*(lambdas[(variables*2-1),cluster]*tmpD_mean+lambdas[(variables*2),cluster]*tmpD_centered)
      
      }
    }
  }
  return(SSQ)
}
