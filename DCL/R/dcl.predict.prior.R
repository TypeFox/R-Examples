dcl.predict.prior<-function(Ntriangle,Xtriangle,
                            inflat.i,inflat.j,Qi,Model=2,adj=2,Tail=FALSE,
                            Tables=TRUE,summ.by="diag",num.dec=2)
{
  m<-nrow(Ntriangle)
  Xtriangle.adj<-Xtriangle
  
  if (missing(inflat.j)) {
    if (Tail==TRUE) inflat.j<-rep(1,2*m-1) else inflat.j<-rep(1,m)
  } else {
    if (length(inflat.j)<m) {
      message("Not valid development inflation, it will be not used")
      if (Tail==TRUE) inflat.j<-rep(1,2*m-1) else inflat.j<-rep(1,m)
    } else {
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
      inflat.i<-inflat.i/inflat.i[1] # normalize as our inflation
      inflat.i.correct<-inflat.i
      inflat.i.correct[inflat.i==0]<-1 
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
  par.dcl<-dcl.estimation(Xtriangle.adj,Ntriangle,adj,Tables=FALSE)  
  pj<-par.dcl$pj
  d<-m-1
  pi.delay<-par.dcl$pi.delay
  ## warning if the DCL model is not working
  if (max(abs(pi.delay))>1) {
    message("Warning: Check if the data follow the DCL model, the estimated delay function is not valid.")
  }
  mu<-par.dcl$mu
  mu.adj<-par.dcl$mu.adj
  Xhat<-as.matrix(par.dcl$Xhat)
  alpha.X<-par.dcl$alpha.X
  beta.X<-par.dcl$beta.X
  alpha.N<-par.dcl$alpha.N
  beta.N<-par.dcl$beta.N
  Nhat<-as.matrix(par.dcl$Nhat)
  
  ## if inflat.i (fix.inflat==TRUE) is provided then we do not estimate it (as in the BDCL paper)
  # update dcl.par if neccesary
  par.dcl.1<-par.dcl
  if (fix.inflat==TRUE){
    inflat<-inflat.i
    inflat.adj<-inflat/inflat.zero
    par.dcl.1$inflat<-inflat.adj
  }
  inflat.adj<-par.dcl.1$inflat
  
  # 2. Pointwise predictions
  
  preds<-dcl.predict(par.dcl.1,Ntriangle,Model,Tail=Tail,Tables=FALSE,summ.by="cell")
  Xrbns<-preds$Xrbns
  Xibnr<-preds$Xibnr
  Xtotal<-preds$Xtotal
  ## CLM as a benchmark (no prior)
  Xhat<-par.dcl$Xhat
  Xhat<-cbind(Xhat,matrix(NA,m,m-1))
  for (i in 1:m){
    for (j in 1:(2*m-1)){
      Xrbns[i,j]<-Xrbns[i,j]*inflat.j[j]*inflat.zero[i]
      Xibnr[i,j]<-Xibnr[i,j]*inflat.j[j]*inflat.zero[i]
      Xtotal[i,j]<-Xtotal[i,j]*inflat.j[j]*inflat.zero[i]
      Xhat[i,j]<-Xhat[i,j]*inflat.zero[i]*inflat.j[j]
    }
  }
  
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
      #Dclm<-c(Dclm,rep(NA,dd-length(Dclm)),Total.CL.paid)
      Dclm<-c(Dclm,Total.CL.paid)
      ## Output: Table with predicted reserve by diagonals 
      
      table.diags<-data.frame(Future.years=c(1:dd,'Tot.'),rbns=round(Drbns,num.dec),
                              ibnr=round(Dibnr,num.dec),total=round(Dtotal,num.dec),
                              clm=round(Dclm,num.dec))
      print(table.diags)
      
    } 
    return(list(par.dcl.upd=par.dcl.1,Drbns=Drbns,Xrbns=Xrbns,Dibnr=Dibnr,Xibnr=Xibnr,
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
      Ri.CL.paid<-rowSums(Xhat)-rowSums(Xtriangle,na.rm=T)
      Total.CL.paid<- sum(Ri.CL.paid,na.rm=TRUE)
      Rclm<-c(Ri.CL.paid,Total.CL.paid)
      
      ## Output: Table with predicted reserves 
      table.rows<-data.frame(Rows=c(1:m,'Tot.'),rbns=round(Rrbns,num.dec),
                             ibnr=round(Ribnr,num.dec),total=round(Rtotal,num.dec),
                             clm=round(Rclm,num.dec))
      print(table.rows)
    } 
    return(list(par.dcl.upd=par.dcl.1,Rrbns=Rrbns,Xrbns=Xrbns,Ribnr=Ribnr,Xibnr=Xibnr,
                Rtotal=Rtotal,Xtotal=Xtotal))
  } 
  
  if (summ.by=="cell") 
  {
    return(list(par.dcl.upd=par.dcl.1,Xrbns=Xrbns,Xibnr=Xibnr,Xtotal=Xtotal))
  }
  
  ## the function ends
}