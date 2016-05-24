dcl.boot<-function(dcl.par,sigma2,Ntriangle,boot.type=2,B=999,Tail=TRUE,
                   summ.by="diag",Tables=TRUE,num.dec=2,n.cal=NA)
{
  ## 1. Parameters in the model
  if (missing(sigma2)) sigma2<-dcl.par$sigma2
  
  if (is.na(sigma2) | sigma2<=0) stop('The variance estimate is null or negative. See the documentation for suggestions.')
  
  Ntriangle<-as.matrix(Ntriangle)
  mu.adj<-dcl.par$mu.adj
  Vy<-dcl.par$Vy
  Ey<-dcl.par$Ey
  pj<-dcl.par$pj
  adj<-dcl.par$adj
  m<-length(pj);d<-m-1
  Nhat<-dcl.par$Nhat
  alpha.N<-dcl.par$alpha.N
  beta.N<-dcl.par$beta.N
  if (missing(Ntriangle)) Ntriangle<-Nhat
  
  # 2. Pointwise predictions
  preds<-dcl.predict(dcl.par,Ntriangle,Model=2,Tail,Tables=FALSE)
  Drbns<-preds$Drbns
  Rrbns<-preds$Rrbns
  Dibnr<-preds$Dibnr
  Ribnr<-preds$Ribnr
  Dtotal<-preds$Dtotal
  Rtotal<-preds$Rtotal
  
  # 3. Bootrapping RBNS and IBNR and total reserves
  # Generate the bootstrap distribution (b=1,...,B)
  
  # Only the variance process (boot I)
  BootI<-function(b,Ntriangle,Nhat,m,d,pj,Ey,Vy)
  {
    # b is the number of the bootstrap sample
    # here we do not estimate the parameters again form the bootstrapped triangles
    delta<- Ey^2 / Vy
    nu<- Ey/Vy
    
    ## Similar to BootII below but simpler since we have not to bootstrappint the parameters
    delta_star<- delta
    nu_star<- nu
    Npred_star<-Nhat
    pj_star<-pj
    
    #  IBNR bootstrapped IBNR_star 
    Npred_star[row(Npred_star)+col(Npred_star)<=m+1]<-NA
    Nest_ext<-cbind(Npred_star,matrix(NA,nrow=m,ncol=d))
    
    Nlist_star<-apply(Nest_ext, MARGIN=1:2, function(v) if (!is.na(v)) 
      as.vector(rmultinom(n=1,size=v,prob=pj_star)) else NA) 
    
    Narray_star<-array(NA,dim=c(m,m+d,d+1))
    for (i in 2:m)
    {  
      for (j in (m-i+2):m)
      {
        Narray_star[i,j,]<-rmultinom(n=1,size=Npred_star[i,j],prob=pj_star)
      }
    }
    Nibnr_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 2:m) 
    {
      for (j in (m-i+2):(m+d)) 
      {
        lim2<- j:(max(m-i+2,j-d))
        for (ll in lim2)
        {
          Nibnr_star[i,j]<- sum(c(Nibnr_star[i,j],Narray_star[i,ll,j-ll+1]),na.rm=T) 
        }
      }
    } 
    rm(Narray_star)
    
    Xibnr_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 2:m) 
    {
      for (j in (m-i+2):(m+d))  
      {
        Xibnr_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nibnr_star[i,j]*delta_star[i]) 
      }
    }
    ## RBNS bootstrapped (Xrbns_star) from Nrbns_star 
    
    # Generate the delay and get Nrbns_star  from observed incurred counts (Ntriangle) 
    
    Nlist_star<-apply(Ntriangle, MARGIN=1:2, FUN=function(v) { 
      if (!is.na(v)) as.vector(rmultinom(n=1,size=v,prob=pj_star)) else NA} )
    # Nlist_star is a matrix m*m with lists at each cell
    Narray_star<-array(NA,dim=c(m,m,d+1))
    for (i in 1:m) 
    {	
      for (j in 1:(m-i+1))
      {
        Narray_star[i,j,]<-rmultinom(n=1,size=Ntriangle[i,j],prob=pj_star)
      }
    }
    
    Nrbns_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 1:m) 
    {
      for (j in (m-i+2):(m+d-i+1)) 
      {
        lim2<- (max(1,j-d)):(max(m-i+1,j-d))
        for (ll in lim2)
        {
          Nrbns_star[i,j]<- sum(c(Nrbns_star[i,j],Narray_star[i,ll,j-ll+1]),na.rm=T) 
        }
      }
    }
    rm(Narray_star)
    Xrbns_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 1:m) 
    {
      for (j in (m-i+2):(m+d-i+1)) 
      {
        Xrbns_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nrbns_star[i,j]*delta_star[i]) 
      }
    }
    
    ## Return the RBNS and IBNR forecasts into a list
    return(list(Xibnr_star=Xibnr_star,Xrbns_star=Xrbns_star))
  }
  
  #Taking into account the uncertainty of the parameters
  BootII<-function(b,Ntriangle,alpha.N,beta.N,m,d,pj,Ey,Vy)
  {
    # b is the number of the bootstrap sample
    # Ntriangle are the original N data to generate Ntriangle_star 
    # which we use to get (new) bootstrapped counts parameters (only IBRN boot.)
    set.triangle<-expand.grid(1:m,1:m) 
    
    ## The next is only for IBNR
    # Generate Ntriangle_star from Pois(Ntriangle) 
    # Ntriangle_star<-matrix(NA,nrow=m,ncol=m) 
    # we complete only the upper triangle to estimate again the N's parameters
    
    Npois<-apply(set.triangle,MARGIN=1,
                 FUN= function(v) { 
                   j<-as.numeric(v[1])
                   i<-as.numeric(v[2])
                   if (j<(m-i+2)) 
                   { 
                     if (Ntriangle[i,j]>0) v.e<-rpois(n=1,lambda=round(Ntriangle[i,j])) else v.e<-0
                   } else v.e<-NA
                   return(v.e) 
                 })
    Ntriangle_star<-matrix(Npois,m,m,byrow=T)
    # Estimate the counts parameters from Ntriangle_star
    #    then we get bootstrapped count parameters
    clm.N.star<-clm(Ntriangle_star,n.cal=n.cal) 
    Npred_star<-clm.N.star$triangle.hat
    
    
    ## The next is common for IBNR and RBNS
    
    ## 1. Generate the delay and get Npaid*  from upper triangle in Ntriangle
    
    Nlist_star<-apply(Ntriangle, MARGIN=1:2, 
                      function(v) {
                        if (!is.na(v)) { 
                          if (v>0) return(as.vector(rmultinom(n=1,size=v,prob=pj))) else return(rep(0,d)) 
                        }})
    #Nlist_star is a matrix m*m with lists at each cell
    
    v.Npaid<-apply(set.triangle,MARGIN=1, FUN= function(v){
      j<-as.numeric(v[1]);i<-as.numeric(v[2]);
      if (j<(m-i+2)) {lim.m<-0:(min(d,j-1))
                      v.n<- sapply(lim.m, function(vv) unlist(Nlist_star[i,j-vv])[vv+1]); 
                      v.n<-sum(v.n,na.rm=T)} else v.n<-NA
      return(v.n) })
    
    Npaid_star<-matrix(v.Npaid,m,m,byrow=T)
    Npaid_star<-cbind(Npaid_star,matrix(NA,nrow=m,ncol=d))
    
    ## 2. Generate new Xs (Xtriangle*) from the estimated distribution
    #    consider only a distribution gamma with no zero claims
    Xtriangle_star<-matrix(NA,nrow=m,ncol= m)
    ## the parameters delta and nu are m-dimensional vectors
    delta<- Ey^2 / Vy
    nu<- Ey/Vy
    
    v.X<-apply(set.triangle,MARGIN=1,FUN= function(v) {
      j<-as.numeric(v[1]);i<-as.numeric(v[2]);
      if (j<(m-i+2)) v.X<-rgamma(n=1,scale=1/nu[i],shape=Npaid_star[i,j]*delta[i]) else v.X<-NA
      return(v.X)})
    Xtriangle_star<-matrix(v.X,m,m,byrow=T)
    
    ## 3. From Xtriangle_star and Ntriangle get the bootstrapped parameters: pj*,Ey*,Vy*
    
    par.star<-dcl.estimation(Xtriangle_star,Ntriangle,adj,Tables=FALSE,n.cal=n.cal)
    pj_star<-par.star$pj
    Ey_star<-as.vector(par.star$Ey)
    mu_star<-Ey_star[1]
    Vy_star<-as.vector(par.star$Vy)
    
    
    ## 4. Generate the paymments X* in rbns and ibnr (gamma distribution)
    if (any(Vy_star<=0)) Vy_star<-Vy
    delta_star<- Ey_star^2 / Vy_star
    nu_star<- Ey_star/Vy_star
    
    ##  IBNR bootstrapped (IBNR_star) from lower triangle in Npred_star and pj_star
    
    Npred_star[row(Npred_star)+col(Npred_star)<=m+1]<-NA
    Nest_ext<-cbind(Npred_star,matrix(NA,nrow=m,ncol=d))
    
    Nlist_star<-apply(Nest_ext, MARGIN=1:2, function(v) if (!is.na(v)) 
      as.vector(rmultinom(n=1,size=v,prob=pj_star)) else NA) 
    
    Narray_star<-array(NA,dim=c(m,m+d,d+1))
    for (i in 2:m)
    {	
      for (j in (m-i+2):m)
      {
        Narray_star[i,j,]<-rmultinom(n=1,size=Npred_star[i,j],prob=pj_star)
      }
    }
    Nibnr_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 2:m) 
    {
      for (j in (m-i+2):(m+d)) 
      {
        lim2<- j:(max(m-i+2,j-d))
        for (ll in lim2)
        {
          Nibnr_star[i,j]<- sum(c(Nibnr_star[i,j],Narray_star[i,ll,j-ll+1]),na.rm=T) 
          # it is j-m+1 because my delay has index which start in l=1 intead of l=0
        }
      }
    } 
    rm(Narray_star)
    
    Xibnr_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 2:m) 
    {
      for (j in (m-i+2):(m+d))	
      {
        Xibnr_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nibnr_star[i,j]*delta_star[i]) 
      }
    }
    
    ## RBNS bootstrapped (Xrbns_star) from Nrbns_star and bootstrapped parameters
    
    # Generate the delay and get Nrbns_star  from observed incurred counts (Ntriangle) 
    # but with bootstrapped parameters pj_star
    
    Nlist_star<-apply(Ntriangle, MARGIN=1:2, FUN=function(v) {
      if (!is.na(v)) as.vector(rmultinom(n=1,size=v,prob=pj_star)) else NA} )
    # Nlist_star is a matrix m*m with lists at each cell
    Narray_star<-array(NA,dim=c(m,m,d+1))
    for (i in 1:m) 
    {	
      for (j in 1:(m-i+1))
      {
        Narray_star[i,j,]<-rmultinom(n=1,size=Ntriangle[i,j],prob=pj_star)
      }
    }
    
    Nrbns_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 1:m) 
    {
      for (j in (m-i+2):(m+d-i+1)) 
      {
        lim2<- (max(1,j-d)):(max(m-i+1,j-d))
        for (ll in lim2)
        {
          Nrbns_star[i,j]<- sum(c(Nrbns_star[i,j],Narray_star[i,ll,j-ll+1]),na.rm=T) 
          # put j-m+1 becuase the index  0 in the delay is the first element
        }
      }
    }
    rm(Narray_star) 
    Xrbns_star<-matrix(NA,nrow=m,ncol= m+d)
    for (i in 1:m) 
    {
      for (j in (m-i+2):(m+d-i+1)) 
      {
        Xrbns_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nrbns_star[i,j]*delta_star[i]) 
      }
    }
    
    
    ## Return the RBNS and IBNR forecasts into a list
    return(list(Xibnr_star=Xibnr_star,Xrbns_star=Xrbns_star))
  }
  
  
  app<-logical(0)
  array.rbns.boot<-array.ibnr.boot<-array(0,dim=c(m,m+d,B))
  for (b in 1:B)
  {
    ## overwrite each simulation
    if (b==1) {app<-F;
               print('Please wait, simulating the distribution...')
    } else app<-T
    
    if (boot.type==1)  res.boot<-BootI(b,Ntriangle,Nhat,m,d,pj,Ey,Vy)
    
    if (boot.type==2)  res.boot<-BootII(b,Ntriangle,alpha.N,beta.N,m,d,pj,Ey,Vy)
    
    array.rbns.boot[,,b]<-as.matrix(res.boot$Xrbns_star)    
    array.ibnr.boot[,,b]<-as.matrix(res.boot$Xibnr_star)
    if (b==B) print('Done!')
  }
  
  ## Bootstrapped RBNS and IBNR claims
  if (summ.by=='diag')
  {  
    # By diagonals
    dd<-m+d-1 +1
    Mat_rbns<-Mat_ibnr<-matrix(0,B,dd)
    for (b in 1:B)
    {
      Xrbns_star<-array.rbns.boot[,,b]
      Xibnr_star<-array.ibnr.boot[,,b]
      if (Tail==FALSE) {Xrbns_star[,(m+1):( 2*m-1)]<-0;Xibnr_star[,(m+1):( 2*m-1)]<-0}
      Drbns_star<-sapply(split(Xrbns_star, row(Xrbns_star)+col(Xrbns_star)), sum, na.rm=T)
      Drbns_star<-as.vector(Drbns_star[-(1:m)])
      Mat_rbns[b,]<-c(Drbns_star,sum(Drbns_star,na.rm =T))
      Dibnr_star<-sapply(split(Xibnr_star, row(Xibnr_star)+col(Xibnr_star)), sum, na.rm=T)
      Dibnr_star<-as.vector(Dibnr_star[-(1:m)])
      Mat_ibnr[b,]<-c(Dibnr_star,sum(Dibnr_star,na.rm =T))
    }
  } else {
    dd<-m+1
    Mat_rbns<-Mat_ibnr<-matrix(0,B,dd)
    for (b in 1:B)
    {
      Xrbns_star<-array.rbns.boot[,,b]
      Rrbns_star<-as.vector(rowSums(Xrbns_star))
      Mat_rbns[b,]<-c(Rrbns_star,sum(Rrbns_star,na.rm =T))
      
      Xibnr_star<-array.ibnr.boot[,,b]
      Ribnr_star<-as.vector(rowSums(Xibnr_star))
      Mat_ibnr[b,]<-c(Ribnr_star,sum(Ribnr_star,na.rm =T))
    }
  }  
  # Totals (ibnr+rbns)
  Mat_total<-Mat_rbns+Mat_ibnr ## has dd rows
  
  ####  BOOTSTRAP SUMMARY: MEAN,SD,QUANTILES
  # RBNS
  mean.Drbns<-colMeans(Mat_rbns) 
  mean.Drbns<-round(mean.Drbns,num.dec)
  sd.Drbns<-apply(Mat_rbns,2,sd)
  sd.Drbns<-round(sd.Drbns,num.dec)
  quantiles.Drbns<-apply(Mat_rbns,2,quantile,c(0.01,0.05,0.5,0.95,0.99))
  quantiles.Drbns<-round(quantiles.Drbns,num.dec)
  # IBNR
  mean.Dibnr<-colMeans(Mat_ibnr) 
  mean.Dibnr<-round(mean.Dibnr,num.dec)
  sd.Dibnr<-apply(Mat_ibnr,2,sd)
  sd.Dibnr<-round(sd.Dibnr,num.dec)
  quantiles.Dibnr<-apply(Mat_ibnr,2,quantile,c(0.01,0.05,0.5,0.95,0.99))
  quantiles.Dibnr<-round(quantiles.Dibnr,num.dec)
  # TOTAL= ibnr+rbns
  mean.Dtotal<-colMeans(Mat_total) 
  mean.Dtotal<-round(mean.Dtotal,num.dec)
  sd.Dtotal<-apply(Mat_total,2,sd)
  sd.Dtotal<<-round(sd.Dtotal,num.dec)
  
  quantiles.Dtotal<-apply(Mat_total,2,quantile,c(0.01,0.05,0.5,0.95,0.99))
  quantiles.Dtotal<-round(quantiles.Dtotal,num.dec)
  
  if (Tables==TRUE)
  {
    period<-c(1:(dd-1),'Tot.')
    resD_rbns<-data.frame(period,rbns=round(Drbns,num.dec),
                          mean.rbns=round(mean.Drbns,num.dec),sd.rbns=round(sd.Drbns,num.dec),
                          Q1.rbns=round(quantiles.Drbns[1,],num.dec),
                          Q5.rbns=round(quantiles.Drbns[2,],num.dec)
                          ,Q50.rbns=round(quantiles.Drbns[3,],num.dec),
                          Q95.rbns=round(quantiles.Drbns[4,],num.dec),
                          Q99.rbns=round(quantiles.Drbns[5,],num.dec))
    resD_ibnr<-data.frame(period,ibnr=round(Dibnr,num.dec),
                          mean.ibnr=round(mean.Dibnr,num.dec),sd.ibnr=round(sd.Dibnr,num.dec),
                          Q1.ibnr=round(quantiles.Dibnr[1,],num.dec),
                          Q5.ibnr=round(quantiles.Dibnr[2,],num.dec),
                          Q50.ibnr=round(quantiles.Dibnr[3,],num.dec),
                          Q95.ibnr=round(quantiles.Dibnr[4,],num.dec),
                          Q99.ibnr=round(quantiles.Dibnr[5,],num.dec))
    resD_total<-data.frame(period,total=round(Dtotal,num.dec),
                           mean.total=round(mean.Dtotal,num.dec),sd.total=round(sd.Dtotal,num.dec),
                           Q1.total=round(quantiles.Dtotal[1,],num.dec),
                           Q5.total=round(quantiles.Dtotal[2,],num.dec),
                           Q50.total=round(quantiles.Dtotal[3,],num.dec),
                           Q95.total=round(quantiles.Dtotal[4,],num.dec),
                           Q99.total=round(quantiles.Dtotal[5,],num.dec))
    
    print(resD_rbns)
    print(resD_ibnr)
    print(resD_total)
    res<-list(array.rbns.boot=array.rbns.boot,array.ibnr.boot=array.ibnr.boot,
              Mat.rbns=Mat_rbns,Mat.ibnr=Mat_ibnr,Mat.total=Mat_total,
              summ.rbns=resD_rbns,summ.ibnr=resD_ibnr,summ.total=resD_total)  
  } else {
    res<-list(array.rbns.boot=array.rbns.boot,array.ibnr.boot=array.ibnr.boot,
                    Mat.rbns=Mat_rbns,Mat.ibnr=Mat_ibnr,Mat.total=Mat_total)  
  }
  return(res)
  
}
