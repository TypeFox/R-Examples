##
##conditional maximisation on 2 dim
##
cm2<- function(rtot=NULL,ctot=NULL,m=matrix(1,length(rtot),length(ctot)),tol=1e-05,maxit=500,verbose=TRUE)
{
  if(round(sum(rtot))!=round(sum(ctot))) 
    stop("row and column totals are not equal, ensure sum(rtot)==sum(ctot)")
  i<-dim(m)[1];  j<-dim(m)[2]
  alpha <- rep(1,i)
  beta <- rep(1,j)
  if(verbose==TRUE){
    rd<-paste("%.",nchar(format(tol,scientific=FALSE))-2,"f",sep="")
    cat(sprintf(rd,c(alpha,beta)), fill = TRUE)
  }
  alpha.old <- alpha+1; beta.old <- beta+1
  it<-1;  max.diff<-tol*2
  while(max.diff>tol & it<maxit ){
    beta.old <- beta
    for(j in 1:j) {
      beta[j] <- ctot[j]/sum(alpha * m[, j])
    }
    alpha.old <- alpha
    for(i in 1:i) {
      alpha[i] <- rtot[i]/sum(beta * m[i,  ])
    }
    it<-it+1
    max.diff<-max(abs(alpha-alpha.old), abs(beta-beta.old))
    if(verbose==TRUE)
      cat(sprintf(rd,c(alpha,beta)), fill = TRUE)
  }
  return(list(N=alpha%*%t(beta)*m,
              theta=c(mu=1,alpha=alpha,beta=beta)))
}

##
##conditional maximization in 3 dimensions
##
cm3<-function(rtot=NULL,ctot=NULL,m,tol=1e-05,maxit=500,verbose=FALSE)
{
  if(round(sum(rtot))!=round(sum(ctot)))
    stop("row and column totals are not equal, ensure sum(rtot)==sum(ctot)")
  i<-dim(m)[1];  j<-dim(m)[2]
  alpha <- rep(1,i)
  beta <- rep(1,j)
  if(verbose==TRUE){
    rd<-paste("%.",nchar(format(tol,scientific=FALSE))-2,"f",sep="")
    cat(sprintf(rd,c(alpha,beta)), fill = TRUE)
  }
  alpha.old <- alpha+1; beta.old <- beta+1
  it<-1;  max.diff<-tol*2
  while(max.diff>tol & it<maxit){
    beta.old <- beta
    for(j in 1:j) {
      beta[j] <- ctot[j]/colSums(alpha*apply(m,c(1,2),sum))[j]
    }
    alpha.old <- alpha
    for(i in 1:i) {
      alpha[i] <- rtot[i]/colSums(beta*apply(m,c(2,1),sum))[i]
    }
    it<-it+1
    max.diff<-max(abs(alpha-alpha.old), abs(beta-beta.old))
    if(verbose==TRUE)
      cat(sprintf(rd,c(alpha,beta)), fill = TRUE)
  }
  return(list(N=c(alpha%*%t(beta))*m,
              theta=c(mu=1,alpha=alpha,beta=beta)))
}


##
##ipf (rather than cm)...gives different estimates
##
ipf2<-function(rtot=NULL,ctot=NULL,m=matrix(1,length(rtot),length(ctot)),tol=1e-05,maxit=500,verbose=FALSE){
  if(!is.null(rtot) & !is.null(ctot))
    if(round(sum(rtot))!=round(sum(ctot))) 
      stop("row and column totals are not equal, ensure sum(rtot)==sum(ctot)")
  n<-list(i=rtot,
          j=ctot)
  mu<-m
  mu.marg<-n
  m.fact<-n
  if(verbose==TRUE)
    rd<-paste("%.",nchar(format(tol,scientific=FALSE))-2,"f",sep="")
  it<-0; max.diff<-tol*2
  while(it==0 | max.diff>tol & it<maxit ){
    if(!is.null(ctot)){
      mu.marg$j <- apply(mu,2,sum)
      m.fact$j <- n$j/mu.marg$j
      m.fact$j[is.nan(m.fact$j)]<-0
      m.fact$j[is.infinite(m.fact$j)]<-0
      mu <- sweep(mu, 2, m.fact$j, "*")
    }
    if(!is.null(rtot)){
      mu.marg$i <- apply(mu,1,sum)
      m.fact$i <- n$i/mu.marg$i
      m.fact$i[is.nan(m.fact$i)]<-0
      m.fact$i[is.infinite(m.fact$i)]<-0
      mu <- sweep(mu, 1, m.fact$i, "*")
    }
    it<-it+1
    #max.diff<-max(abs(unlist(n)-unlist(mu.marg)))
    #speeds up a lot if get rid of unlist (new to v1.7)
    max.diff<-max(abs(c(n$i-mu.marg$i, n$j-mu.marg$j)))
    
    if(verbose==TRUE)
      cat(sprintf(rd,unlist(m.fact)), fill = TRUE)
  }
  return(list(mu=mu,it=it,tol=max(abs(unlist(n)-unlist(mu.marg)))))
}

##
##ipf3
##
ipf3<-function(rtot=NULL,ctot=NULL,m=NULL,tol=1e-05,maxit=500,verbose=TRUE){
  if(any(round(colSums(rtot))!=round(rowSums(ctot))))
    stop("row and column totals are not equal for one or more sub-tables, ensure colSums(rtot)==rowSums(ctot)")
  n<-list(ik=rtot,
          jk=t(ctot))
  R<-dim(rtot)[1]
  
  #set up offset
  if(is.null(m)){
    m<-array(1,c(dim(rtot),dim(rtot)[1]))
    dimnames(m)<-list(orig=dimnames(rtot)[[1]],dest=dimnames(ctot)[[1]],pob=dimnames(rtot)[[2]])
  }
  
  mu<-m
  mu.marg<-n
  m.fact<-n
  it<-0; max.diff<-tol*2
  while(max.diff>tol & it<maxit){
    mu.marg$ik <- apply(mu,c(1,3),sum)
    m.fact$ik <- n$ik/mu.marg$ik
    m.fact$ik[is.nan(m.fact$ik)]<-0
    m.fact$ik[is.infinite(m.fact$ik)]<-0
    mu <- sweep(mu, c(1,3), m.fact$ik, "*")
    
    mu.marg$jk <- apply(mu, c(2,3), sum)
    m.fact$jk <- n$jk/mu.marg$jk
    m.fact$jk[is.nan(m.fact$jk)]<-0
    m.fact$jk[is.infinite(m.fact$jk)]<-0
    mu <- sweep(mu, c(2,3), m.fact$jk, "*")
    
    it<-it+1
    #max.diff<-max(abs(unlist(n)-unlist(mu.marg)))
    #speeds up a lot if get rid of unlist (new to v1.7)
    max.diff<-max(abs(c(n$ik-mu.marg$ik, n$jk-mu.marg$jk)))
    
    if(verbose==TRUE)
      cat(c(it, max.diff), "\n")
  }
  return(list(mu=mu,it=it,tol=max.diff))
}

##
##flow matrix
##
fm<-function(y){
  R<-dim(y)[3]
  dg<-diag(apply(y,c(1,2),sum))
  od<-apply(y,c(1,2),sum)
  if(R==dim(y)[1])  fl<-od-diag(dg,R,R)
  if(R!=dim(y)[1]){
    y.adj<-y
    for(i in 1:R){
      y.adj[i,,i]<-y[i,,i]+y[R+1,,i]
    }
    od<-apply(y.adj,c(1,2),sum)
    fl<-od[1:R,1:R]-diag(dg[1:R],R,R)
  }
  diag(fl)<-0
  fl
}


##
##ipf3.qi
##
# P1.adj=P1;P2.adj=P2
#rtot=t(P1.adj);ctot=P2.adj;dtot=NULL;verbose=TRUE;tol=1e-05;maxit=500;speed=TRUE;m=NULL
ipf3.qi<-function(rtot=NULL,ctot=NULL,dtot=NULL,m=NULL,speed=TRUE,tol=1e-05,maxit=500,verbose=TRUE){
  if(any(round(colSums(rtot))!=round(rowSums(ctot))))
    stop("row and column totals are not equal for one or more sub-tables, ensure colSums(rtot)==rowSums(ctot)")
  
  n<-list(ik=rtot,
          jk=t(ctot),
          ijk=dtot)
  R<-dim(rtot)[1]
  
  #set up diagonals
  df1<-expand.grid(a = 1:R, b = 1:R)
  if(is.null(dtot)){
    dtot<-array(1,c(R,R,R))
    dtot<-with(df1, replace(dtot, cbind(a, a, b),  apply(cbind(c(n$ik),c(n$jk)),1,min) ))
    n$ijk<-dtot
  }
  
  #set up offset
  if(is.null(m)){
    m<-array(1,c(dim(rtot),dim(rtot)[1]))
  }
  if(is.null(dimnames(m))){
    dimnames(m)<-list(orig=dimnames(rtot)[[1]],dest=dimnames(ctot)[[1]],pob=dimnames(rtot)[[2]])
  }
  
  #alter ss (to speed up)
  if(speed==TRUE){
    n$ik<-n$ik-(apply(n$ijk,c(1,3),sum)-(R-1))
    n$jk<-n$jk-(apply(n$ijk,c(2,3),sum)-(R-1))
    n$ijk<-with(df1, replace(n$ijk, cbind(a, a, b),  0 ))
  }
  
  mu<-m
  mu.marg<-n
  m.fact<-n
  it<-0; max.diff<-tol*2
  while(max.diff>tol & it<maxit){
    mu.marg$ik <- apply(mu,c(1,3),sum)
    m.fact$ik <- n$ik/mu.marg$ik
    m.fact$ik[is.nan(m.fact$ik)]<-0
    m.fact$ik[is.infinite(m.fact$ik)]<-0
    mu <- sweep(mu, c(1,3), m.fact$ik, "*")
    
    mu.marg$jk <- apply(mu, c(2,3), sum)
    m.fact$jk <- n$jk/mu.marg$jk
    m.fact$jk[is.nan(m.fact$jk)]<-0
    m.fact$jk[is.infinite(m.fact$jk)]<-0
    mu <- sweep(mu, c(2,3), m.fact$jk, "*")
    
    mu.marg$ijk <- with(df1, replace(n$ijk, cbind(a, a, b),  c(apply(mu,3,diag)) ))
    m.fact$ijk <- n$ijk/mu.marg$ijk
    m.fact$ijk[is.nan(m.fact$ijk)]<-0
    m.fact$ijk[is.infinite(m.fact$ijk)]<-0
    mu <- mu*m.fact$ijk
    
    it<-it+1
    #max.diff<-max(abs(unlist(n)-unlist(mu.marg))) 
    #speeds up a lot if get rid of unlist (new to v1.6)
    max.diff<-max(abs(c(n$ik-mu.marg$ik, n$jk-mu.marg$jk, n$ijk-mu.marg$ijk)))
    if(verbose==TRUE)
      cat(c(it, max.diff), "\n")
  }
  if(speed==TRUE){
    mu<-with(df1, replace(mu, cbind(a, a, b), c(sapply(1:R, function(i) diag(dtot[,,i]))) ))
  }
  return(list(mu=mu,it=it,tol=max.diff))
}
#rm(n,mu,mu.marg,m.fact)
#ipf3.qi(rtot=t(P1.adj),ctot=P2.adj,m=m)#

##
##flows from stock
##
#P1=P1;P2=P2;d=d;b=b.adj;m=m
#method="stocks"; d.mat=NULL;b.mat=NULL;b.deduct="native.gt0"
#m=NULL;
ffs<-function(P1,P2,d,b,m=NULL,method="stocks",b.mat=NULL,d.mat=NULL,b.deduct="native.gt0",...){
  if(!(method %in% c("outside","stocks","deaths")) | length(method)!=1)
    stop("method must be one of outside, stocks or deaths")
  if(!(b.deduct %in% c("native.only","native.gt0")) | length(b.deduct)!=1)
    stop("method must be one of outside, stocks or deaths")
  R<-nrow(P1)
  #set up offset
  if(is.null(m)){
    m<-array(1,c(dim(P1),dim(P1)[1]))
  }
  if(is.null(dimnames(m))){
    dimnames(m)<-list(orig=dimnames(P1)[[1]],dest=dimnames(P2)[[1]],pob=dimnames(P1)[[2]])
  }
  if(method=="stocks"){
    if(round(sum(P2-P1),1)!=round(sum(b-d),1)){
      message("sum(P2-P1): ", sum(P2-P1))
      message("sum(b-d):   ", sum(b-d))
      stop("difference in stock tables must be equal to sum of all births - sum of all deaths")
    }
  } 
  
  #set up migration matrix with births/deaths and others rows and columns
  y<-m
  y<-array(0,dim(m)+c(2,2,0))
  dimnames(y)<-list(o=c(dimnames(m)[[1]],"b","O"),d=c(dimnames(m)[[1]],"d","O"),b=dimnames(m)[[3]])
  y[,,]<-0
  
  #step 1-2a take off deaths
  message("Adjust for Deaths...")
  if(is.null(d.mat))
    d.mat<-ipf2(ctot=d,m=P1)$mu
  if(method=="deaths"){
    d.mat<-ipf2(ctot=d,rtot=rowSums(P1) + b - rowSums(P2),m=P1)$mu  
  }
  y[1:R,R+1,]<-t(d.mat)
  P1.adj<-P1-d.mat
  
  #step 1-2b take off births
  message("Adjust for Births...")
  if(is.null(b.mat))
    b.mat<-diag(b)
  y[R+1,1:R,]<-b.mat
  P2.adj<-P2-b.mat
  #adjust for negative nb sums after birth deduction...spread births to all population where needed (new to v1.6)
  if(b.deduct=="native.gt0"){
    ii<-diag(P2.adj<0)
    if(sum(ii)>0){
      b.mat[,ii]<-ipf2(ctot=b,m=P2)$mu[,ii]
      P2.adj<-P2-b.mat
    }
  }
  
  #step 3-4a take off moves in from external or adjust P1.adj rows
  message("Remaining Adjustments (P1)...")
  dif<-rowSums(P1.adj) - rowSums(P2.adj)
  if(method=="outside" | method=="deaths"){
    #following is in versions <1.3. is wrong. those leaving contolled for in P1, like those who die
    #this (in the #) is labelled in.mat but should be out.mat (where the dif>0). should have offset P1.adj, where they leave from, not P2.adj
    #in.mat<-t(ipf2(ctot=pmax(dif,0),m=t(P2.adj))$mu)
    #P1.adj<-P1.adj-in.mat
    #y[R+2,1:R,]<-t(in.mat)
    out.mat<-t(ipf2(ctot=pmax(dif,0),m=t(P1.adj))$mu)
    P1.adj<-P1.adj-out.mat
    y[1:R,R+2,]<-t(out.mat)
  }
  if(method=="stocks"){
    #was not reaching convergence of near zero difference in maxdiff in <v1.5, i.e. either row or col totals in estimates were not matching arguments.
    P1.adj<-ipf2(rtot=rowSums(P1.adj)-dif/2,ctot=colSums(P1.adj), m=P1.adj, maxit=100000, tol=0.1)$mu
  }
  
  #step 3-4b take off moves out from external or adjust P2.adj rows
  message("Remaining Adjustments (P2)...")
  if(method=="outside" | method=="deaths"){
    #following is in versions <1.3. is wrong. those arriving contolled for in P2, like those who are born
    #this (in the #) is labelled out.mat but should be in.mat (where the dif<0). should have offset P2.adj, where they arrive too, not P1.adj
    #out.mat<-t(ipf2(ctot=pmax(-dif,0),m=t(P1.adj))$mu)
    #P2.adj<-P2.adj-out.mat
    #y[1:R,R+2,]<-t(out.mat)
    in.mat<-t(ipf2(ctot=pmax(-dif,0),m=t(P2.adj))$mu)
    P2.adj<-P2.adj-in.mat
    y[R+2,1:R,]<-t(in.mat)
  }
  if(method=="stocks"){
    #was not reaching convergence of near zero difference in maxdiff in <v1.5, i.e. either row or col totals in estimates were not matching arguments.
    P2.adj<-ipf2(rtot=rowSums(P2.adj)+dif/2,ctot=colSums(P2.adj),m=P2.adj, maxit=100000, tol=0.1)$mu
  }
  
  #step 5 calculate
  message("Calculate Flows...")
  #ipf<-ipf3.qi(rtot=t(P1.adj),ctot=P2.adj,m=m)
  ipf<-ipf3.qi(rtot=t(P1.adj),ctot=P2.adj,m=m,...)
  y[1:R,1:R,]<-ipf$mu
  return(c(ipf,list(y=y)))
}


##
##rogers castro curve
##
rc9.fund<-list(a1=0.02, alpha1=0.1, a2=0.06, alpha2=0.1, mu2=20, lambda2=0.4, c=0.003)

rc9<-function(x, param=rc9.fund, scaled=TRUE){
  if(!is.list(param))
    stop("param must be a list")
  if(sum(is.na(match(names(param),names(rc9.fund))))!=0)
    stop("param must be a list with correct names, see for example rc9.fund")
  m<-param$a1*exp(-param$alpha1*x)+param$a2*exp(-param$alpha2*(x-param$mu2)-exp(-param$lambda*(x-param$mu2)))+param$c
  if(scaled==TRUE) m<-m/sum(m)
  m
}

##
##weighted funtion for b.mat and d.mat
##
