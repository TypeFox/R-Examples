
# res=Prepare(x,simplify,qua,standardize)
Prepare =function (x, simplify, qua, standardize){
  ind=nrow(x@M)
  vars=ncol(x@M)
  MM=list()
  if (simplify){
    p=(0:qua)/qua
    for (j in 1:vars){
      MM[[j]]=matrix(0,(qua+1),(ind+1))
      MM[[j]][,(ind+1)]=p
      for (i in 1:ind){
        dom=numeric(0)
        for (q in 0:qua){
          dom=c(dom, compQ(x@M[i,j][[1]],q/qua))
        }
        
        tmp=new("distributionH",dom,p) 
        tmp=(tmp-tmp@m)*(x@M[i,j][[1]]@s/tmp@s)+x@M[i,j][[1]]@m #transformation with invariance with respect mean and std
        x@M[i,j][[1]]=tmp
        MM[[j]][,i]=tmp@x
      }
    } 
  }
  else{
    
    for (j in 1:vars){
      tmp=registerMH(x[,j])
      MM[[j]]=matrix(0,length(tmp@M[1,1][[1]]@x),ind+1)
      MM[[j]][,ind+1]=tmp@M[1,1][[1]]@p
      for (i in 1:ind){
        x@M[i,j][[1]]=tmp@M[i,1][[1]]
        MM[[j]][,i]=tmp@M[i,1][[1]]@x
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
          MM[[v]][,i]=(MM[[v]][,i]-Mc[v])/STAND[v]
          #x@M[i,v][[1]]@x=(x@M[i,v][[1]]@x-Mc[v])/STAND[v]
        }
      }
    }
    
  }
  return(list(MM=MM,x=x))
}

# compute fast SSQ
ComputeFastSSQ=function(subMM){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  SSQ=0
  m1=apply(as.matrix(subMM[,1:ind]),1,mean)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  cm=(m1[2:rr]+m1[1:(rr-1)])/2
  rm=(m1[2:rr]-m1[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
    ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
    SSQ=SSQ+sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2))
  }
  return(list(SSQ=SSQ,mx=m1, mp=subMM[,ind+1]))
}

# compute fast DIST
ComputeFast_L2_SQ_WASS_D=function(subMM){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  
  c1=(subMM[2:rr,1]+subMM[1:(rr-1),1])/2
  r1=(subMM[2:rr,1]-subMM[1:(rr-1),1])/2
  c2=(subMM[2:rr,2]+subMM[1:(rr-1),2])/2
  r2=(subMM[2:rr,2]-subMM[1:(rr-1),2])/2
  
  Dist=sum(p1*((c1-c2)^2)+ p1/3*((r1-r2)^2))
  
  return(Dist)
}

WH.ADPT.KMEANS.TOTALSSQ=function(x,memb,m,lambdas,proto){
  vars=ncol(x@M)
  ind=nrow(x@M)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  #Compute general weights for the global PROTOTYPE
  mu=apply(memb^m,2,sum)
  protoGEN=x[1,]
  for (variables in (1:vars)){
    protoGEN@M[1,variables][[1]]=WH.vec.mean(proto[,variables],mu)
  }
  #Compute general weights for the global PROTOTYPE2
  Mat_of_means=matrix(0,ind,vars)
  CentredD=x
  for (i in 1:ind){
    for (j in 1:vars){
      Mat_of_means[i,j]=x@M[i,j][[1]]@m
      CentredD@M[i,j][[1]]=x@M[i,j][[1]]-x@M[i,j][[1]]@m
    }
  }
  W_centers=matrix(0,ind,vars)
  W_dist_cent=matrix(0,ind,vars)
  for (i in 1:ind){
    for (j in 1:vars){
      W_centers[i,j]=lambdas[(j*2-1),which.max(memb[i,])]
      W_dist_cent[i,j]=lambdas[(j*2),which.max(memb[i,])]
    }
  }
  protoGEN2=x[1,]
  mprot=matrix(0,1,vars)
  for (variables in (1:vars)){
    mprot[1,variables]=sum(Mat_of_means[,variables]*W_centers[,variables])/sum(W_centers[,variables])
    protoGEN2@M[1,variables][[1]]=WH.vec.mean(CentredD[,variables],W_dist_cent[,variables])+mprot[1,variables]
  }
  #Compute BETWEEN
  BSQ_clu=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (variables in (1:vars)){
    GENPROTm=protoGEN2@M[1,variables][[1]]@m
    GENPROTcent=protoGEN2@M[1,variables][[1]]-GENPROTm
    for (clu in 1:k){
      LOCPROTm=proto@M[clu,variables][[1]]@m
      LOCPROTcent=proto@M[clu,variables][[1]]-LOCPROTm
      BSQ_clu[1,variables,clu]=sum(memb[,clu])*lambdas[(variables*2-1),clu]*(GENPROTm-LOCPROTm)^2
      BSQ_clu[2,variables,clu]=sum(memb[,clu])*lambdas[(variables*2),clu]*WassSqDistH(GENPROTcent,LOCPROTcent)
    }
  }
  #Compute the total fuzzy sum of SQUARES
  TSQ_clu=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],protoGEN@M[1,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        TSQ_clu[1,variables,cluster]=TSQ_clu[1,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        TSQ_clu[2,variables,cluster]=TSQ_clu[2,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
      }
    }
  }
  TSQ_clu2=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],protoGEN2@M[1,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        TSQ_clu2[1,variables,cluster]=TSQ_clu2[1,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        TSQ_clu2[2,variables,cluster]=TSQ_clu2[2,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
      }
    }
  }
  return(RES=list(protoGEN=proto,TSQdetailed=TSQ_clu, TSQ=sum(TSQ_clu),TSQdetailed2=TSQ_clu2, TSQ2=sum(TSQ_clu2),
                  BSQdetailed=BSQ_clu, BSQ=sum(BSQ_clu)))
}

ComputeFastSSQ_Fuzzy=function(subMM,memb,m){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  SSQ=0
  W=memb^m
  W2=W/sum(W)
  m1=subMM[,1:ind]%*%as.matrix(W2)
  
  m1=as.vector(m1)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  cm=(m1[2:rr]+m1[1:(rr-1)])/2
  rm=(m1[2:rr]-m1[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
    ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
    SSQ=SSQ+(sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2)))*W[indiv]
  }
  return(list(SSQ=SSQ,mx=m1, mp=subMM[,ind+1]))
}

ComputeFast_Fuzzy_TOT_SSQ=function(MM,proto,memb,m){
  ind=ncol(MM[[1]])-1
  var=length(MM)
  k=ncol(memb)
#computeGenProt
mus=apply(memb^m,2,sum)
mus=mus/sum(mus)
GP=new('MatH',1,var)
for (variable in 1:var){
  dG=MM[[variable]][,1]*0
  for (clu in 1:k){
    dG=dG+proto@M[clu,variable][[1]]@x*mus[clu]
  }
  tmp=new('distributionH',x=dG,p=MM[[variable]][,(ind+1)])
  GP@M[1,variable][[1]]=tmp
}
#compute TOTALSSQ
SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
for (variable in 1:var){
  rr=nrow(MM[[variable]])
  p1=MM[[variable]][2:rr,ind+1]-MM[[variable]][1:(rr-1),ind+1]
  cm=(GP@M[1,variable][[1]]@x[2:rr]+GP@M[1,variable][[1]]@x[1:(rr-1)])/2
  rm=(GP@M[1,variable][[1]]@p[2:rr]-GP@M[1,variable][[1]]@p[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(MM[[variable]][2:rr,indiv]+MM[[variable]][1:(rr-1),indiv])/2
    ri=(MM[[variable]][2:rr,indiv]-MM[[variable]][1:(rr-1),indiv])/2
    dist=sum(p1*((ci-cm)^2+1/3*(ri-rm)^2))
    dc=(sum(p1*ci)-GP@M[1,variable][[1]]@m)^2
      dv=dist-dc
    for (clu in 1:k){
      SSQ_det[1,variable,clu]=SSQ_det[1,variable,clu]+memb[indiv,clu]^m*dc
      SSQ_det[2,variable,clu]=SSQ_det[2,variable,clu]+memb[indiv,clu]^m*dv
    }
  }
}
return(list(SSQ=sum(SSQ_det),SSQ_det=SSQ_det,ProtoGEN=GP))
}

ComputeFast_Fuzzy_Adaptive_TOT_SSQ=function(MM,proto,memb,m,lambdas){
  ind=ncol(MM[[1]])-1
  var=length(MM)
  k=ncol(memb)
  #computeGenProt
  
  GP=new('MatH',1,var)
  for (variable in 1:var){
    rr=nrow(MM[[variable]])
    ci=(MM[[variable]][2:rr,1:ind]+MM[[variable]][1:(rr-1),1:ind])/2
    ri=(MM[[variable]][2:rr,1:ind]-MM[[variable]][1:(rr-1),1:ind])/2
    wi=(MM[[variable]][2:rr,(ind+1)]-MM[[variable]][1:(rr-1),(ind+1)])
    LMc=memb^m*matrix(lambdas[(variable*2-1),],ind,k, byrow=TRUE)
    LMc=LMc/sum(LMc)
    LMv=memb^m*matrix(lambdas[(variable*2),],ind,k, byrow=TRUE)
    LMv=LMv/sum(LMv)
    Cen=MM[[variable]][,1]*0
    MG=0
    McG=ci[,1]*0
    RcG=ri[,1]*0
    for (clu in 1:k){
      for (i in 1:ind){
      mui=sum(ci[,i]*wi)
      MG=MG+LMc[i,clu]*mui
      Cen=Cen+(MM[[variable]][,i]-mui)*LMv[i,clu]
      }
    }
    Cen=Cen+MG
    tmp=new('distributionH',x=Cen,p=MM[[variable]][,(ind+1)])
    GP@M[1,variable][[1]]=tmp
  }
  #compute TOTALSSQ
  SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
  for (variable in 1:var){
    rr=nrow(MM[[variable]])
    p1=MM[[variable]][2:rr,ind+1]-MM[[variable]][1:(rr-1),ind+1]
    cm=(GP@M[1,variable][[1]]@x[2:rr]+GP@M[1,variable][[1]]@x[1:(rr-1)])/2
    rm=(GP@M[1,variable][[1]]@p[2:rr]-GP@M[1,variable][[1]]@p[1:(rr-1)])/2
    for (indiv in 1:ind){
      ci=(MM[[variable]][2:rr,indiv]+MM[[variable]][1:(rr-1),indiv])/2
      ri=(MM[[variable]][2:rr,indiv]-MM[[variable]][1:(rr-1),indiv])/2
      dist=sum(p1*((ci-cm)^2+1/3*(ri-rm)^2))
      dc=(sum(p1*ci)-GP@M[1,variable][[1]]@m)^2
      dv=dist-dc
      for (clu in 1:k){
        SSQ_det[1,variable,clu]=SSQ_det[1,variable,clu]+lambdas[(variable*2-1),clu]*(memb[indiv,clu]^m)*dc
        SSQ_det[2,variable,clu]=SSQ_det[2,variable,clu]+lambdas[(variable*2),clu]*(memb[indiv,clu]^m)*dv
      }
    }
  }
  return(list(SSQ=sum(SSQ_det),SSQ_det=SSQ_det,ProtoGEN=GP))
}

#'From real data to distributionH.
#' 
#' @param data a set of numeric values.
#' @param algo (optional) a string. Default is "histogram", i.e. the function "histogram"
#' defined in the \code{\link[histogram]{histogram}}  package. \cr If "base" 
#' the \code{\link[graphics]{hist}} function is used. \cr
#' "FixedQuantiles" computes the histogram using as breaks a fixed number of quantiles.\cr
#' "ManualBreaks" computes a histogram where braks are provided as a vector of values.\cr
#' "PolyLine" computes a histogram using a piecewise linear approximation of the empirical
#' cumulative distribution function using the "Ramer-Douglas-Peucker algorithm", 
#'  \url{http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm}. 
#'  An \code{epsilon} parameter is required.
#'  The data are scaled in order to have a standard deviation equal to one.
#' @param type (optional) a string. Default is "combined" and generates 
#' a histogram having regularly spaced breaks (i.e., equi-width bins) and 
#' irregularly spaced ones. The choice is done accordingly with the penalization method described in 
#' \code{\link[histogram]{histogram}}. "regular" returns equi-width binned histograms, "irregular" returns
#' a histogram without equi-width histograms. 
#' @param qua a positive integer to provide if \code{algo="FixedQuantiles"} is chosen. Default=10.
#' @param breaks a vector of values to provide if  \code{algo="ManualBreaks"} is chosen.
#' @param epsilon a number between 0 and 1 to provide if \code{algo="PolyLine"} is chosen. Default=0.01.
#' @return A \code{distributionH} object, i.e. a distribution.
#' @importFrom histogram histogram
#' @importFrom stats quantile sd
#' @export
#' @examples
#' data=rnorm(n = 1000,mean = 2,sd = 3)
#' mydist=data2hist(data)
#' plot(mydist)
#' @seealso \code{\link[histogram]{histogram}} function
data2hist<-function(data, 
                    algo="histogram",
                    type="combined",
                    qua=10,
                    breaks=numeric(0),
                    epsilon=0.01){
  a=switch(algo,"base"=1, "histogram"=2, "FixedQuantiles"=3,"ManualBreaks"=4 ,
           "PolyLine"=5,2)
  t=switch(type,"regular"=1, "irregular"=2, "combined"=3, 3)
  if (a==1){
    h <- hist(data)
    x=h$breaks
    counts=h$counts
    counts[which(counts==0)]=0.001
    p=c(0,cumsum(counts/sum(counts)))
  }
  if (a==2){
   
    if (t==1){
      h<-histogram(data,type="regular", verbose=FALSE,plot=FALSE)
    }
    if (t==2){
      h<-histogram(data,type="irregular",verbose=FALSE,plot=FALSE)
    }
    if (t==3){
      h<-histogram(data,type="combined",verbose=FALSE,plot=FALSE)
    }
    x=h$breaks
    counts=h$counts
    counts[which(counts==0)]=0.001
    p=c(0,cumsum(counts/sum(counts)))
  }
  if (a==3){
    p=c(0:qua)/qua
    x=rep(0,qua+1)
    for (i in 0:qua){      
      x[i+1]=as.numeric(quantile(data,probs = p[i+1]))
      if (i>0){
        if(x[i+1]<=x[i]) x[i+1]=x[i]+(1e-10)
      }
    }    
  }
  if (a==4){
    #checking limits
    
    ini=min(data)
    end=max(data)
    rr=end-ini
    tol=min(1e-10,(rr*(1e-10)))
    end=end+tol
    breaks=sort(unique(c(breaks,ini,end)))
    breaks=breaks[which(breaks>=ini)]
    breaks=breaks[which(breaks<=end)]
    #end check
    data.cut = cut(data, breaks, right=FALSE)
    counts=as.vector(table(data.cut))
    x=breaks
    p=c(0,cumsum(counts/sum(counts)))
  }
  if (a==5){
    data=sort(data)
    stdev=sd(data)
    
    cums=c(0:(length(data)-1))/(length(data)-1)
    points=cbind(data/stdev,cums)
    resu=DouglasPeucker(as.matrix(points),epsilon=epsilon)
    x=as.vector(resu[,1])*stdev
    p=as.vector(resu[,2])
  }
  p[1]=0
  p[length(p)]=1
  mydist<- distributionH(x,p)
  return(mydist)
}
#'Ramer-Douglas-Peucker algorithm for curve fitting with a PolyLine
#' 
#' @param points a 2D matrix with the coordinates of 2D points
#' @param epsilon an number between 0 and 1. Recomended 0.01.
#' @return A matrix with the points of segments of a Poly Line.
#' @export
#' @seealso \code{\link[HistDAWass]{data2hist}} function
#' @export
DouglasPeucker=function(points,epsilon){
  dmax=0
  index=0
  end=nrow(points)
  ResultList=numeric(0)
  if (end<3) return (ResultList=rbind(ResultList,points))
  for (i in 2:(end-1)){
    d=ShortestDistance(points[i,], line=rbind(points[1,],points[end,]))
    if (d>dmax){
      index=i
      dmax=d
    }
  }
  #if dmax is greater than epsilon recursively apply
  if (dmax>epsilon){
   # print(dmax)
    recResults1=DouglasPeucker(points[1:index,],epsilon)
    recResults2=DouglasPeucker(points[index:end,],epsilon)
    ResultList=rbind(ResultList,recResults1,recResults2)
   
  }
  else
  {
    ResultList=rbind(ResultList,points[1,],points[end,])
  }
  ResultList=as.matrix(ResultList[!duplicated(ResultList),])
  colnames(ResultList)=c("x","p")
  return(ResultList)
}
#' Shortes distance from a point o a 2d segment
#' 
#' @param p coordinates of a point
#' @param line a 2x2 matrix with the coordinates of two points defining a line
#' @return A numeric value, the Euclidean distance of point \code{p} to the \code{line}.
#' @export
#' @seealso \code{\link[HistDAWass]{data2hist}} function and \code{\link[HistDAWass]{DouglasPeucker}} function
#' @export
ShortestDistance=function(p, line){
  x1=line[1,1]
  y1=line[1,2]
  x2=line[2,1]
  y2=line[2,2]
  x0=p[1]
  y0=p[2]
  d=abs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)
  return(as.numeric(d))
}
