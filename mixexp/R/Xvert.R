Xvert = function(nfac=3,uc=c(0,0),lc=c(0,0),nlc=0,lb=c(0,0),ub=c(0,0),coef,ndm=0,plot=TRUE,
                 cornerlabs = c("x1","x2","x3"), axislabs = c("x1","x2","x3"),pseudo=TRUE) 
{
  # checks for the number of factors
  if (nfac>12){
    stop(" The maximum number of factors allowed is 12")
  }
  # checks to see if mxaimum number of constraints is less than 12
  n.uc<-length(uc)
  n.lc<-length(lc)
  if (max(n.uc,n.lc)>12) {
    stop(" The maximum number of mixture components is 12")
  }
  if(max(n.uc,n.lc)==0){
    stop(" No constraints given")
  }
  if (n.uc != n.lc) {
    stop(" the number of upper constraints supplied is different than the number of lower constraints")
  }
  # Create the constraints matrix for crvtave
  #  ck<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
  ck<-c(lc[1],uc[1])
  for (i in 2:n.uc){
    ck<-cbind(ck,c(lc[i],uc[i]))
  }
  for (i in (n.uc+1):12) {
    ck<-cbind(ck,c(0,1))
  }
  
  for (i in 1:12){
    cks<-ck[1,i]+(1-ck[2,i])
    if (cks!=0) {
      nfacc=i 
    } else {break}
  }
  if (nfacc>nfac) {
    stop(" The number of upper and lower limits supplied exceeds the number of factors")
  }
  
  v<-c(-ck[1,1],ck[2,1])
  for (i in 2:nfac) {
    v<-c(v,-ck[1,i],ck[2,i])
  }
  
  
  #Creates conmx corresponding to upper and lower constraints on components
  Ip<-diag(nfac)
  In<--1*Ip
  conmx<-interleave(Ip,In)
  conmx<-cbind(conmx,v)
  
  # Create constraint matrix for linear constraints
  if (nlc>0) {
    loc<-nrow(coef)
    if(loc != nlc){
      stop(" The number of rows of the coefficient matrix must equal the number of linear constraints")
    }
    loc2<-ncol(coef)
    if (loc2 != nfac){
      stop (" The number of columns of the coefficient matrix must equal the number of mixture components")
    }
    lolb<-length(lb)
    if (lolb!=nlc) { 
      stop(" The number of lower bounds for linear constraints is not equal to the number of linear constraints")
    }
    loub<-length(ub)
    if (loub!=nlc) { 
      stop(" The number of upper bounds for linear constraints is not equal to the number of linear constraints")
    }
    #    lincon<-matrix(coef,byrow=T,nrow=nlc)
    lincon<-coef
    
    nlcon1<-nrow(lincon)
    
    
    #set negative of coef on top of coef in a matrix
    
    nlinc<- -1*lincon
    lincon<-rbind(nlinc,lincon)
    
    # set upper bounds on top of negative of lower bounds in a vector
    v<- ub
    v<-c(v,-lb)
    
    # add vector of bounds to the right of lincon
    lincon<-cbind(lincon,v)
    
    # append lincon to the bottom of conmx
    conmx<-rbind(conmx,lincon)
  }
  
  # delete rows where contraint is zero
  conmx<-conmx[abs(conmx[,nfac+1])>0, ]
  
  # calls crvtave to create exteme vertices design plus centroid
  des<-crvtave(ndm,conmx)
  des<-data.frame(des)
  
  
  if (nfac==3 & plot) {
    DesignPoints(des[ ,1:3],x1lower=lc[1],x1upper=uc[1],x2lower=lc[2],x2upper=uc[2],x3lower=lc[3],x3upper=uc[3],cornerlabs=cornerlabs,axislabs=axislabs,pseudo=pseudo) }
  #  }
  return(des)
}
####End Function #################################