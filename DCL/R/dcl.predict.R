dcl.predict<-function(dcl.par,Ntriangle,Model=2,Tail=TRUE,Tables=TRUE,
                      summ.by="diag",num.dec=2)
{
  inflat<-dcl.par$inflat
  mu<-dcl.par$mu
  mu.adj<-dcl.par$mu.adj
  pi.delay<-dcl.par$pi.delay
  ## warning if the DCL model is not working
  if (max(abs(pi.delay))>1) {
    message("Warning: Check if the data follow the DCL model, the estimated delay function is not valid.")
  }
  
  pj<-dcl.par$pj
  m<-length(pj);d<-m-1
  Nhat<-dcl.par$Nhat
  Nhat<-as.matrix(Nhat)
  Xhat<-dcl.par$Xhat
  Xhat<-as.matrix(Xhat)
  if (missing(Ntriangle) | Model==0) Ntriangle<-Nhat 
  Ntriangle<-as.matrix(Ntriangle)
  
  two.models<-function(pj,mu)
  {
    ## The functions Ec.rbns and Ec.ibnr calculate 
    # the conditional expectation in the RBNS and IBNR prediction sets
    Ec.rbns<-function(Ntriangle,pj,mu,inflat,m,d) 
    {
      # Expected RBNS:  i=1,...,m ; j=m-i+2,...,m+d-i+1  TRIANGLE J1+J_2 
      # X_{ij}^{rbns}=mu * sum_{l=(i-m+j-1)^(min((j-1),d))  p_l N_{i,j-l}
      set.rect<-expand.grid(1:(m+d),1:m) #the whole rectangle dim=(m,m+d)
      v.expect<-apply(set.rect,MARGIN=1,FUN= function(v) { j<-v[1];i<-v[2];
                                                           if ((j>(m-i+1))& (j<(m+d-i+2))) {limq<-(i-m+j-1):(min((j-1),d));limN<-j-limq;
                                                                                            v.e<-sum(Ntriangle[i,limN]* pj[limq+1] * mu*inflat[i],na.rm=T)} else v.e<-NA
                                                           return(v.e) })
      Xrbns<-matrix(v.expect,m,m+d,byrow=T)
      return(Xrbns)
    }
    
    Ec.ibnr<-function(Nhat,pj,mu,inflat,m,d) 
    {
      # Expected IBNR:  i=2,...,m ; j= m-i+2,..., m+d TRIANGLE J1+(J_2-row1)+J3
      # X_{ij}^{ibnr}=sum_{l=0}^min{d,j+i-m} p_l Nhat_{i,j-l}
      set.rect<-expand.grid(1:(m+d),1:m)  #the whole rectangle dim=(m,m+d)
      Nhat_ext<-cbind(Nhat,matrix(NA,nrow=m,ncol=d))
      v.expect<-apply(set.rect,MARGIN=1, FUN=function(v) {
        j<-as.numeric(v[1]);i<-as.numeric(v[2]);
        if ((j>(m-i+1))& (j<(m+d+1))& (i>1)) {
          limq<-0:(min((i-m+j-2),d));limN<-j-limq;
          v.e<-sum(Nhat_ext[i,limN]*pj[limq+1]*mu*inflat[i],na.rm=T)} else v.e<-NA
        return(v.e)})
      Xibnr<-matrix(v.expect,m,m+d,byrow=T)
      return(Xibnr)
    }
    
    Xrbns<-Ec.rbns(Ntriangle,pj,mu,inflat,m,d)
    
    Xibnr<-Ec.ibnr(Nhat,pj,mu,inflat,m,d) 
    if (Tail==FALSE) {Xrbns[,(m+1):( 2*m-1)]<-0 ; Xibnr[,(m+1):( 2*m-1)]<-0}
    Xrbns.zero<-Xrbns
    Xrbns.zero[is.na(Xrbns)]<-0
    Xibnr.zero<-Xibnr
    Xibnr.zero[is.na(Xibnr)]<-0
    Xtotal<-Xrbns.zero+Xibnr.zero
    
    if (summ.by=="diag")
    {
      # Estimated reserve for RBNS claims,Drbns, for each diagonal (l=i+j=m+2,...,m+d)
      dd<-m+d-1
      Drbns<-sapply(split(Xrbns, row(Xrbns)+col(Xrbns)), sum, na.rm=T)
      Drbns<-as.vector(Drbns[-(1:m)])
      Drbns<-c(Drbns,sum(Drbns,na.rm =T))
      # Estimated reserve for IBNR claims: Dibnr
      Dibnr<-sapply(split(Xibnr, row(Xibnr)+col(Xibnr)), sum, na.rm=T)
      Dibnr<-as.vector(Dibnr[-(1:m)])
      Dibnr<-c(Dibnr,sum(Dibnr,na.rm =T))
      
      ## Total reserves (ibnr+rbns) by diagonals
      Dtotal<-Drbns+Dibnr
      
      if (Tables==TRUE)
      {
        # compare with the simple CL predictions from paid data (Xhat)
        ## By diagonals (calendar year): we have m-1 values
        Diag.CL.paid<-sapply(split(Xhat, row(Xhat)+col(Xhat)), sum, na.rm=T)
        Dclm<-c(Diag.CL.paid[-(1:m)])
        Total.CL.paid<- sum(Dclm,na.rm=TRUE)
        Dclm<-c(Dclm,rep(NA,dd-length(Dclm)),Total.CL.paid)
        
        ## Output: Table with predicted reserve by diagonals 
        
        table.diags<-data.frame(Future.years=c(1:dd,'Tot.'),rbns=round(Drbns,num.dec),
                                ibnr=round(Dibnr,num.dec),total=round(Dtotal,num.dec),
                                clm=round(Dclm,num.dec))
        print(table.diags)
        
      } 
      return(list(Drbns=Drbns,Xrbns=Xrbns,Dibnr=Dibnr,Xibnr=Xibnr,
                  Dtotal=Dtotal,Xtotal=Xtotal))
    } 
    
    if (summ.by=="row")
    {
      Rrbns<-rowSums(Xrbns,na.rm=TRUE)
      Rrbns<-c(Rrbns,sum(Rrbns,na.rm=TRUE))
      Ribnr<-rowSums(Xibnr,na.rm=TRUE)
      Ribnr<-c(Ribnr,sum(Ribnr,na.rm=TRUE))
      ## Total reserves (ibnr+rbns) by diagonals
      Rtotal<-Rrbns+Ribnr
      
      if (Tables==TRUE)
      {
        # compare with the simple CL predictions from paid data (Xhat)
        ## By rows (accident year): we have m values
        #Ri.CL.paid<-rowSums(Xhat)-rowSums(Xtriangle,na.rm=T)
        Xhat.low<-Xhat
        Xhat.low[row(Xhat)+col(Xhat)<=(m+1)]<-NA
        Ri.CL.paid<-rowSums(Xhat.low,na.rm=T)
        Total.CL.paid<- sum(Ri.CL.paid,na.rm=TRUE)
        Rclm<-c(Ri.CL.paid,Total.CL.paid)
        
        ## Output: Table with predicted reserves 
        table.rows<-data.frame(Rows=c(1:m,'Tot.'),rbns=round(Rrbns,num.dec),
                               ibnr=round(Ribnr,num.dec),total=round(Rtotal,num.dec),
                               clm=round(Rclm,num.dec))
        print(table.rows)
      } 
      return(list(Rrbns=Rrbns,Xrbns=Xrbns,Ribnr=Ribnr,Xibnr=Xibnr,
                  Rtotal=Rtotal,Xtotal=Xtotal))
    } 
    
    if (summ.by=="cell") 
    {
      return(list(Xrbns=Xrbns,Xibnr=Xibnr,Xtotal=Xtotal))
    }
  }
  if (Model==1 | Model==0) res.model<-two.models(pi.delay,mu) else res.model<-two.models(pj,mu.adj)
  return(res.model)
  
}
