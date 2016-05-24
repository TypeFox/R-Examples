dcl.boot.prior<-function(Xtriangle,Ntriangle,sigma2,mu,
                         inflat.i,inflat.j,Qi,Model=2,adj=1,
                         boot.type=2,B=999,Tail=TRUE,summ.by="diag",
                         Tables=TRUE,num.dec=2,n.cal=NA,Fj.X=NA,Fj.N=NA)
  
{
  m<-nrow(Ntriangle)
  Xtriangle.adj<-Xtriangle
  
  if (missing(inflat.j)) {
    if (Tail==TRUE) inflat.j<-rep(1,2*m-1) else inflat.j<-c(rep(1,m),rep(NA,m-1))
  } else {
    if (length(inflat.j)<m) {
      message("Not valid development inflation (not used here). The length should be 
              the number of columns in the triangle (m) or 2*m-1 if 
              the tail is required. ")
      if (Tail==TRUE) inflat.j<-rep(1,2*m-1) else inflat.j<-c(rep(1,m),rep(NA,m-1))
    } else {
      if (length(inflat.j)<2*m-1){
        if (Tail==TRUE) message("The length of the dev. inflation should be 2*m-1 if 
                                the tail is required. You are not getting the tail. ")
        inflat.j<-c(inflat.j,rep(NA,m-1))
      }
      inflat.j.correct<-inflat.j
      inflat.j.correct[inflat.j==0]<-1 
      ##  Remove the inflat.j effect 
      Xtriangle.adj<-t( t(Xtriangle.adj)/inflat.j.correct[1:m] )
    }
    }
  
  
  if (missing(inflat.i)) fix.inflat<-FALSE else {
    if (length(inflat.i)<m) {
      message("Not valid underwriting inflation, it will be not used")
      fix.inflat<-FALSE
    } else {
      fix.inflat<-TRUE
    }
  }
  
  if (missing(Qi)) {
    Qi<-rep(0,m)
    zero<-FALSE
  } else {
    if (length(Qi)<m | any(Qi<0) | any(Qi>1)) {
      message("Not valid zero-claims probabilities, they will be not used")
      Qi<-rep(0,m)
      #   Xtriangle.adj<-Xtriangle
      zero<-FALSE
    } else {
      Qi.correct<-Qi
      Qi.correct[Qi==1]<-0 
      #  Remove the zero-claims inflation
      Xtriangle.adj<-Xtriangle.adj/(1-Qi.correct)
      zero<-TRUE
    }  
  }
  inflat.zero<-1-Qi
  
  ## 2. Estimating the model parameters from the adjusted triangle (Xtriangle.adj)
  par.dcl<-dcl.estimation(Xtriangle.adj,Ntriangle,adj,Tables=FALSE,n.cal=n.cal,Fj.X=Fj.X,Fj.N=Fj.N)  
  pj<-par.dcl$pj
  d<-m-1
  pi.delay<-par.dcl$pi.delay
  mu.adj<-par.dcl$mu.adj
  Xhat<-as.matrix(par.dcl$Xhat)
  alpha.X<-par.dcl$alpha.X
  beta.X<-par.dcl$beta.X
  alpha.N<-par.dcl$alpha.N
  beta.N<-par.dcl$beta.N
  Nhat<-as.matrix(par.dcl$Nhat)
  if (missing(mu)) {
    mu.adj<-par.dcl$mu.adj 
    fix.mu<-FALSE
  } else { if (is.na(mu)) {    
    mu.adj<-par.dcl$mu.adj 
    fix.mu<-FALSE
  } else {
    mu.adj<-mu
    fix.mu<-TRUE
  }
  }
  
  if (missing(sigma2)){
    sigma2<-par.dcl$sigma2
    fix.sigma2<-FALSE
  } else {
    if (is.na(sigma2)) {sigma2<-par.dcl$sigma2
                        fix.sigma2<-FALSE} else  fix.sigma2<-TRUE
  }  
  
  ## if inflat.i (fix.inflat==TRUE) is provided then we do not estimate it (as in the BDCL paper)
  if (fix.inflat==TRUE){
    inflat<-inflat.i/inflat.zero
    inflat.correct<-inflat
    inflat.correct[inflat.correct==0]<-1
  } else{
    inflat<-par.dcl$inflat
    inflat.correct<-inflat  ## to multiply the variance and avoid it is zero
  }
  
  if (is.na(sigma2) | sigma2<=0){
    stop('There is not enough information to estimate the variance of the 
         individual size claims. See the documentation for suggestions.')
  } else {
    Ey<-mu.adj*inflat #inflat.correct
    # and Vy= inflat^2 * (sigma2*(1-Qi) - Qi* mu^2) the variance
    Vy<-inflat.correct^2 * (sigma2*(1-Qi) - Qi* mu.adj^2) #if Qi=0 then is just DCL-variance
    
  }
  
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
        if (zero==F)  Xibnr_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nibnr_star[i,j]*delta_star[i]) 
        if (zero==T)
        {
          numpay<-rbinom(n=1,size=Nibnr_star[i,j],prob=1-Qi[i])
          Xibnr_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=numpay*delta_star[i])
        }
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
        if (zero==F)	Xrbns_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nrbns_star[i,j]*delta_star[i]) 
        if (zero==T)
        {
          numpay<-rbinom(n=1,size=Nrbns_star[i,j],prob=1-Qi[i])
          Xrbns_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],
                                  shape=numpay*delta_star[i])
        }
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
    clm.N.star<-clm(Ntriangle_star,n.cal) 
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
    
    
    ## 2. Generate a new-bootstrapped X triangle (Xtriangle_star with dimension m) 
    # We use a the gamma distribution with bootstrapped parameters
    # note thatthe INFLATION makes that the gamma is different for each row
    # then the parameters are vectors with dimension m (rows) as happend with Ey, Vy (both vectors)
    delta<- Ey^2 / Vy
    nu<- Ey/Vy
    
    if (zero==F){
      ### original DCL
      v.X<-apply(set.triangle,MARGIN=1,
                 FUN= function(v)
                 { j<-as.numeric(v[1]);i<-as.numeric(v[2]);
                   if (j<(m-i+2)) v.X<-rgamma(n=1,scale=1/nu[i],shape=Npaid_star[i,j]*delta[i])  else v.X<-NA
                   ## we do not need to include a filter when size=0, in this case it returns a zeros vector
                   return(v.X) } )
      # the X's values are into v.X but we need to put properly into a m-dimensional matrix                  
      Xtriangle_star<-matrix(v.X,m,m,byrow=T)
    } else {
      ## nuevo geramos la distribucion mixta
      ##if (zero==T)
      
      ## Indiv Payments puede ser 0 con probabilidad Q o la gamma con probabilidad 1-Q
      ## Entonces genero para cada una de las Npaid_star[i,j] claims un par (Z,Y) donde Z->B(1,1-Q) y Y->Gamma(delta,nu) and define
      ## the IndPay[i,j,m]=Z*Y for m=1,...,Nrbns[i,j]
      ## Esto es lo mismo que generar el total de pagos no zero: numpay->B(Nrbns[i,j],1-Q)
      ## y luego calcular igual que antes: Xtriangle_star[i,j]=sum( rgamma(n=numpay,delta, nu))
      v.X<-apply(set.triangle,MARGIN=1,
                 FUN= function(v)
                 { j<-as.numeric(v[1]);i<-as.numeric(v[2]);
                   if (j<(m-i+2)) {
                     numpay<-rbinom(n=1,size=Npaid_star[i,j],prob=1-Qi[i]);
                     v.X<-rgamma(n=1,scale=1/nu[i],shape=numpay*delta[i])
                   }else v.X<-NA
                   ## we do not need to include a filter when size=0, in this case it returns a zeros vector
                   return(v.X) } )
      # the X's values are into v.X but we need to put properly into a m-dimensional matrix                  
      Xtriangle_star<-matrix(v.X,m,m,byrow=T)
    }  
    
    ## We again adjust the bootstrapped triangle to estimate the parameters using DCL
    Xtriangle_star_adj<-Xtriangle_star/(1-Qi)  
    
    
    
    ## 3. From Xtriangle_star and Ntriangle get the bootstrapped parameters: pj*,Ey*,Vy*
    
    par.star<-dcl.estimation(Xtriangle_star,Ntriangle,adj,Tables=FALSE,n.cal=n.cal) 
    pj_star<-par.star$pj
    mu.adj.star<-par.star$mu.adj
    sigma2.star<-par.star$sigma2
    inflat.star<-par.star$inflat
    
    if (fix.mu==TRUE) mu.adj.star<-Ey[1]
    if (fix.inflat==TRUE) inflat.star<-Ey/Ey[1]  
    if (fix.sigma2==TRUE) sigma2.star<-Vy[1]
    
    
    Ey_star<-mu.adj.star*inflat.star
    Vy_star<-inflat.star^2 * (sigma2.star*(1-Qi) - Qi* mu.adj.star^2)
    
    
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
        if (zero==F)  Xibnr_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nibnr_star[i,j]*delta_star[i]) 
        if (zero==T)
        {
          numpay<-rbinom(n=1,size=Nibnr_star[i,j],prob=1-Qi[i])
          Xibnr_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=numpay*delta_star[i])
        }
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
        if (zero==F)	Xrbns_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=Nrbns_star[i,j]*delta_star[i]) 
        if (zero==T)
        {
          numpay<-rbinom(n=1,size=Nrbns_star[i,j],prob=1-Qi[i])
          Xrbns_star[i,j]<-rgamma(n=1,scale=1/nu_star[i],shape=numpay*delta_star[i])
        }
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
    
    xrbns<-as.matrix(res.boot$Xrbns_star)
    xibnr<-as.matrix(res.boot$Xibnr_star)
    ## Calculate the predictions by diagonals after including inflat.j 
    array.rbns.boot[,,b]<-t( t(xrbns)*inflat.j )
    array.ibnr.boot[,,b]<-t( t(xibnr)*inflat.j )
    
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
    resD_rbns<-data.frame(period,
                          mean.rbns=round(mean.Drbns,num.dec),sd.rbns=round(sd.Drbns,num.dec),
                          Q1.rbns=round(quantiles.Drbns[1,],num.dec),
                          Q5.rbns=round(quantiles.Drbns[2,],num.dec)
                          ,Q50.rbns=round(quantiles.Drbns[3,],num.dec),
                          Q95.rbns=round(quantiles.Drbns[4,],num.dec),
                          Q99.rbns=round(quantiles.Drbns[5,],num.dec))
    resD_ibnr<-data.frame(period,
                          mean.ibnr=round(mean.Dibnr,num.dec),sd.ibnr=round(sd.Dibnr,num.dec),
                          Q1.ibnr=round(quantiles.Dibnr[1,],num.dec),
                          Q5.ibnr=round(quantiles.Dibnr[2,],num.dec),
                          Q50.ibnr=round(quantiles.Dibnr[3,],num.dec),
                          Q95.ibnr=round(quantiles.Dibnr[4,],num.dec),
                          Q99.ibnr=round(quantiles.Dibnr[5,],num.dec))
    resD_total<-data.frame(period,
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
