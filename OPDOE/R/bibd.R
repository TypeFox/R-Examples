# BIBD-petr.R, still not finished

iterate.factor.comb <- function (p, n, i) 
{
  ret <- NULL
  for(col in 1:n){
    ret <- c(ret, ((i-1)%/%(p^(n-col)))%%p)
  }
  ret
}


matrix.search<-function(A, l)
{
  for (i in 1:nrow(A))
    if (all(A[i,]==l)) {return(i); break;}
  return(0)
}

multi.modp<-function(a,b,p){
  la <- length(a)
  lb <- length(b)
  
  dummy <- numeric(la+lb-1)        
  
  for (i in 1:(la+lb-1)){
    dummy[i] <- sum(a[max(i-lb+1,1):min(i,la)]*b[min(lb,i):max(1,i-la+1)]) %%p
  } 
  return(dummy)
}


myGF<-function(p, n)
{
  x<-GF(p,n)

  obj<-x[[1]]
  pol<-x[[2]][1,]

  if (n>1) {

    addition<-matrix(0,p^n,p^n)

    for(i in 1:(p^n))
      for(j in i:(p^n))
        {
          newobj<-(obj[i,]+obj[j,]) %%p
          addition[i,j]<-addition[j,i]<-matrix.search(obj,newobj)
        }

    multiply<-matrix(0,p^n,p^n)
    for(i in 1:(p^n))
      for(j in i:(p^n))
        {
          newobj<-redu.modp(multi.modp(obj[i,],obj[j,],p),pol,p)
          multiply[i,j]<-multiply[j,i]<-matrix.search(obj,newobj)
        }
  }
  else
    {

      addition<-matrix(0,p^n,p^n)

      for(i in 0:(p-1))
        for(j in 0:(p-1))
          {
            addition[i+1,j+1]<-addition[j+1,i+1]<-(i+j)%%p+1
          }

      multiply<-matrix(0,p^n,p^n)
      for(i in 0:(p-1))
        for(j in 0:(p-1))
          {
            multiply[i+1,j+1]<-multiply[j+1,i+1]<-(i*j)%%p+1
          }

    }


  return(list(objects=obj,add=addition,multiply=multiply))

}

elements.of.matrix<-function(A, i, j)
{
  apply(cbind(i,j),1,function(x) A[x[1],x[2]])
}

################################################################
# PRINT, TEST, ETC...
################################################################

#calls all other methods
bibd<-function(v, k, method=0, ...)
{
  if(v<k) stop("v < k ")
  # automatically choose method
  if(method==0){
    if(k>v/2) method <- 2
    if(v<6 || k<3) method <- 1
    if(v==6 && k==3) method <- 3 # 9
    if(v==7 && k==3) method <- 3 # 9, 16
    if(v==8 && k==3) method <- 6 # 9
    if(v==8 && k==4) method <- 4 # 5, 6
    if(v==9 && k==3) method <- 4 # 6, 9, 18
    if(v==9 && k==4) method <- 6 # 7, 19
    if(v==10 && k==3) method <- 9 # 10
    if(v==10 && k==4) method <- 7 # 18
    if(v==10 && k==5) method <- 18 #
    if(v==11 && k==3) method <- 6 # 7, 9
    if(v==11 && k==4) method <- 6 # 7
    if(v==11 && k==5) method <- 6 # 7, 15
    if(v==12 && k==3) method <- 9 #
    if(v==12 && k==4) method <- 8 #
    if(v==12 && k==5) method <- 15 # 
    if(v==12 && k==6) method <- 5 # 8, 18
    if(v%%4==3) method <- 16
  }
    
  # call method:
  if(method==1) ret <- bibd.method1(v,k)
  if(method==2) ret <- bibd.method2(v,k)
  if(method==3) ret <- bibd.method3(v,k)
  if(method==4) ret <- bibd.method4(v,k)
#  if(method==5) ret <- bibd.method5(v,k)
  if(method==6) ret <- bibd.method6(v,k)
  if(method==7) ret <- bibd.method7(v,k)
  if(method==8) ret <- bibd.method8(v,k)
  if(method==9) ret <- bibd.method9(v,k)
  if(method==10) ret <- bibd.method10(v,k)
  if(method==11) ret <- bibd.method11(v,k)
  if(method==12) ret <- bibd.method12(v,k)
  if(method==13) ret <- bibd.method13(v,k)
  if(method==14) ret <- bibd.method14(v,k)
  if(method==15) ret <- bibd.method15(v,k)
  if(method==16) ret <- bibd.method16(v,k)
  if(method==17) ret <- bibd.method17(v,k)
  if(method==18) ret <- bibd.method18(v,k)
  if(method==19) ret <- bibd.method18(v,k)
    
  class(ret)<-"bibd"
  return(ret)
}

#print bibd
print.bibd<-function(x, width=getOption("width"), ...)
{
  
  # sort by columns 
  x<-x[do.call(order,lapply(1:ncol(x), function(i) x[,i])),]

  # get characteristics
  t<-test.bibd(x)

  vdgts<-ceiling(log(t$v+1)/log(10))
  bdgts<-ceiling(log(t$b+1)/log(10))
#  cat(to.print," block  treatments",sep="")
  cat(" block  treatments\n")
  for (i in 1:nrow(x))
    {
      to.print<-paste(i)
      
      while (nchar(to.print)<bdgts)
        to.print<-paste(" ",to.print,sep="")
      
      cat(to.print,"   (",sep="")
      for (j in 1:ncol(x))  
        {
          to.print<-paste(x[i,j])
          while (nchar(to.print)<vdgts)
            to.print<-paste(" ",to.print,sep="")
          cat(to.print)
          if (j!=ncol(x)) cat(", ") else cat(")\n")
        }


    }
  cat(paste("\n v = ",t$v,
            "   k = ",t$k,
            "   b = ",t$b,
            "   r = ",t$r,
            "   lambda = ", t$lambda,
            "\n\n", sep=""))

}


# compact form -> incidence matrix
compact.2.imatrix<-function(A)
{
  v<-length(unique(as.vector(A)))
  b<-nrow(A)

  B<-matrix(0,v,b)
  for (i in 1:nrow(A))
    B[A[i,],i]<-1
  
  class(B)<-"bibdmatrix"
  return(B)
}

# printing method
print.bibdmatrix<-function(B)
{

  cat("\nIncidence matrix of bibd:\n\n")

  B<-unclass(B)
  NextMethod(B)
}

# incidence matrix form -> compact
imatrix.2.compact<-function(B)
{
  A<-matrix(0,ncol(B),sum(B[,1]))

  for (i in 1:ncol(B))
    A[i,]<-which(B[,i]==1)

  class(A)<-"bibd"
  return(A)
}


test.bibd<-function(A, fast=TRUE)
{
  flag<-TRUE
  abc<-unique(as.vector(A))
  B<-compact.2.imatrix(A)

  v<-length(abc);
  k<-ncol(A);
  b<-nrow(A);
  r<-sum(B[1,]==1)
  lambda<-sum(B[1,]==1&B[2,]==1)

  if(!fast)
    {
      if (any(apply(B,1,sum)!=r)) warning("Desing is not equi-replicated.")
      flag<-TRUE
      for (i in 1:(nrow(B)-1))
        for (j in (i+1):nrow(B))
          if (flag&sum(B[i,]==1&B[j,]==1)!=lambda)
            {flag<-FALSE; warning("Desing is not balanced.")}
    }
  
  return(list(v=v,k=k,b=b,r=r,lambda=lambda)) 
}

# check necessary conditions
necessary.conditions<-function(v, k, b, r, lambda)
{
  if (v*r!=b*k) return(1); #equi and proper
  if (lambda*(v-1)!=r*(k-1)) return(2); #balanc
  if (b<v) return(3) #fisher  
  if (2*lambda*(v*r-b-(v-1))*(b-3)/r/(r-1)/(v*r-b-v+3)<1 & v<6) return(4) #calinski  
  return(0);
}



################################################################
# METHOD 3 - PROJECTIVE GEOMETRY
################################################################

bibd.method3<-function(v, k, verbose=T)  #n,s,h=1,m=n-1)
{
  primen100 <- c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
      61,67,71,73,79,83,89,97)
  s <- 0; n <- 0
  for(p in primen100){

    for(N in 1:100){
      v1 <- (as.bigz(p)^as.bigz(N+1)-1)/as.bigz(p-1)
#      if(v1<MAXv)
      if(verbose)cat(paste("v:",as.numeric(as.character(v1)),"\n"))  
      for(M in (N-1):1){
        k1 <- (as.bigz(p)^as.bigz(M+1)-1)/as.bigz(p-1)
        if(verbose)cat(paste("k:",as.numeric(as.character(k1)),"\n"))
        if(v==v1 && k==k1){
          s=p; n=N; m=M;
          break
        }
      }
      if(s!=0) break
    }
    if(s!=0) break
  }
  
  if(s==0 && n==0) stop("no s and n found")
  h=1#;m=n-1;
  is.null.comb<-function(q,w)
    {
      e<-rep(0,length(q))
      for (i in 1:length(q)) e[i]<-gf$multiply[q[i],w[i]]
      dummy<-1
      for (i in 1:length(q)) dummy<-gf$add[dummy,e[i]]
      return(dummy==1)
    }

  p<-s
  gf<-myGF(p,h)
  if((p^h)^(n+1)>2^32)stop("p and h too large")
  poi<-coef<-factor.comb(p^h,n+1)[-1,]+1

  just.this<-rep(TRUE,nrow(poi))

  for (i in 1:(nrow(poi)-1))
    for (j in 2:(p^h))
      for (k in (i+1):nrow(poi))
        {
          dummy<-elements.of.matrix(gf$multiply,poi[i,],rep(j,ncol(poi)))
          if (all(dummy==poi[k,])) just.this[k]<-FALSE
        }
  poi<-poi[just.this,]

  results<-NULL
  for (i in 1:nrow(coef))
    {
      res<-c()
      for (j in 1:nrow(poi))
        if (is.null.comb(coef[i,],poi[j,])) res<-c(res,j)
      results[[i]]<-res+1  
    }

  k<-max(sapply(results,length))
  for (i in 1:nrow(coef))
    if (length(results[[i]])==k-1) results[[i]]<-c(1,results[[i]])

  results<-unique(results)

  if (m!=n-1) 
    for (i in (n-2):m)
      {
        new.results<-NULL
        index<-0
        for (j in 1:(length(results)-1))
          for (k in (j+1):length(results))
            { index<-index+1
              new.results[[index]]<-intersect(results[[j]],results[[k]])
            }
        l<-max(sapply(new.results,length))
        results<-NULL
        new.index<-0
        for (j in 1:index)
          if (length(new.results[[j]])==l)
            {
              new.index<-new.index+1
              results[[new.index]]<-new.results[[j]]
            }
        results<-unique(results)
      }
  
  results<-matrix(as.numeric(as.factor(unlist(results))),byrow=TRUE,nrow=length(results))
  
  
  return(results)
}


bibd.method3.new<-function(n, s, h=1, m=n-1)
{

  is.null.comb<-function(q,w)
    {
      e<-rep(0,length(q))
      for (i in 1:length(q)) e[i]<-gf$multiply[q[i],w[i]]
      dummy<-1
      for (i in 1:length(q)) dummy<-gf$add[dummy,e[i]]
      return(dummy==1)
    }

  p<-s
  gf<-myGF(p,h)
  poi <- NULL
  npoi <- (p^h)^(n+1)
  mpoi <- n+1
  I <- 1
  while(I<=npoi){
    if(I%%1000==0)cat(paste("I:",I,"\n"))
    poiI <- iterate.factor.comb(p^h,n+1,I)[-1]+1
  # poi<-coef<-factor.comb(p^h,n+1)[-1,]+1

  # just.this<-rep(TRUE,nrow(poi))

#  for (i in 1:(nrow(poi)-1))
    i <- I
    j <- 2
    while(j<=p^h){
    #for (j in 2:(p^h))
      k <- i+1
      while(k<=npoi){
            if(k%%1000==0)cat(paste("k:",k,"\n"))
      #for (k in (i+1):npoi)
       
          poiK <- iterate.factor.comb(p^h,n+1,k)[-1]+1
          dummy<-elements.of.matrix(gf$multiply,poiI,rep(j,mpoi))
          if (!all(dummy==poiK)) # just.this[k]<-FALSE
            poi <- rbind(poi,poiI)
          k <- k+1
        }
      j <- j+1
    }
    I <- I+1
  }
  #poi<-poi[just.this,]
  
  results<-list()
  lresults <- 1
  i <- 1
  while(i<=npoi){
        if(i%%1000==0)cat(paste("i:",i,"\n"))
    poiI <- iterate.factor.comb(p^h,n+1,i)[-1]+1
#  for (i in 1:nrow(coef))
#    {
    res<-c()
    for (j in 1:nrow(poi))
      if (is.null.comb(poiI,poi[j,])) res<-c(res,j)
    results[[lresults]]<-res+1
    lresults <- lresults+1
    i <- i+1
  }
  
  k<-max(sapply(results,length))
  for (i in 1:length(results))
    if (length(results[[i]])==k-1) results[[i]]<-c(1,results[[i]])

  results<-unique(results)

  if (m!=n-1) 
    for (i in (n-2):m)
      {
        new.results<-NULL
        index<-0
        for (j in 1:(length(results)-1))
          for (k in (j+1):length(results))
            { index<-index+1
              new.results[[index]]<-intersect(results[[j]],results[[k]])
            }
        l<-max(sapply(new.results,length))
        results<-NULL
        new.index<-0
        for (j in 1:index)
          if (length(new.results[[j]])==l)
            {
              new.index<-new.index+1
              results[[new.index]]<-new.results[[j]]
            }
        results<-unique(results)
      }
  
  results<-matrix(as.numeric(as.factor(unlist(results))),byrow=TRUE,nrow=length(results))
  
  
  return(results)
}

################################################################
# METHOD 4 - EUCLIDEAN GEOMETRY
################################################################

bibd.method4<-function(n, s, h=1, m=n-1)
{

  is.null.comb<-function(q,w)
    {
      e<-rep(0,length(q))
      for (i in 1:length(q)) e[i]<-gf$multiply[q[i],w[i]]
      dummy<-w[length(q)+1]
      for (i in 1:length(q)) dummy<-gf$add[dummy,e[i]]
      return(dummy==1)
    }

  p<-s
  gf<-myGF(p,h)
  poi<-factor.comb(p^h,n)+1
  coef<-factor.comb(p^h,n+1)[-1,]+1

  results<-NULL
  for (i in 1:nrow(coef))
    {
      res<-c()
      for (j in 1:nrow(poi))
        if (is.null.comb(poi[j,],coef[i,])) res<-c(res,j)
      results[[i]]<-res 
    }

#k<-max(sapply(results,length))
#for (i in 1:nrow(coef))
#  if (length(results[[i]])==k-1) results[[i]]<-c(1,results[[i]])

  results<-unique(results)

  if (m!=n-1) 
    for (i in (n-2):m)
      {
        new.results<-NULL
        index<-0
        for (j in 1:(length(results)-1))
          for (k in (j+1):length(results))
            { index<-index+1
              new.results[[index]]<-intersect(results[[j]],results[[k]])
            }
        l<-max(sapply(new.results,length))
        results<-NULL
        new.index<-0
        for (j in 1:index)
          if (length(new.results[[j]])==l)
            {
              new.index<-new.index+1
              results[[new.index]]<-new.results[[j]]
            }
        results<-unique(results)
      }

  new.results<-results
  l<-max(sapply(new.results,length))
  results<-NULL
  new.index<-0
  for (j in 1:length(new.results))
    if (length(new.results[[j]])==l)
      {
        new.index<-new.index+1
        results[[new.index]]<-new.results[[j]]
      }
  results<-unique(results)

  results<-matrix(as.numeric(as.factor(unlist(results))),byrow=TRUE,nrow=length(results))

  return(results)
}

######################
# bibd repository
######################

#initial setting
#bibd.repo.table<-data.frame(v=NULL,k=NULL,b=NULL,r=NULL,lambda=NULL,id=NULL)
#bibd.repo<-NULL
#bibd.repo.method<-NULL

#load repository
#load("bibd.repo.save")

#save.repo<-function()
#{
#  save(bibd.repo,bibd.repo.method,bibd.repo.table,file="bibd.repo.save")
#}  

#add.to.repo<-function(bibd.to.add, method)
#{
#  temp<-test.bibd(bibd.to.add)
#  N<-length(bibd.repo)
#  founded<-FALSE
#  for (i in seq_along(bibd.repo.table[,1]))
#    if (bibd.repo.table[i,]$v==temp$v & bibd.repo.table[i,]$k==temp$k) 
#      {founded<-TRUE; break;}
#  if (!founded) {
#    bibd.repo[[N+1]]<<-bibd.to.add
#    bibd.repo.method[[N+1]]<<-method
#    bibd.repo.table<<-rbind(bibd.repo.table,
#                            data.frame(v=temp$v,k=temp$k,b=temp$b,r=temp$r,lambda=temp$lambda,id=N+1))
#    return(TRUE)
#  }
#  else if (temp$b<bibd.repo.table[i,]$b) {
#    bibd.repo[[bibd.repo.table[i,]$id]]<<-bibd.to.add
#    bibd.repo.method[[bibd.repo.table[i,]$id]]<-method
#    bibd.repo.table[i,-6]<<-data.frame(v=temp$v,k=temp$k,b=temp$b,r=temp$r,lambda=temp$lambda)
#  } else return(FALSE)
#}

#add.by.call<-function(m)
#{
#  bibd.to.add<-do.call(what=bibd,arg=m)
#  add.to.repo(bibd.to.add,m)
#}

#minimal.bibd<-function(v, k)
#{
#  result<-NULL
#  if (k==2|k==v-1) return(bibd(v=v,k=k,method=1))
#  if (k<=v/2) {
#    for (i in seq_along(bibd.repo.table[,1]))
#      if (bibd.repo.table[i,]$v==v & bibd.repo.table[i,]$k==k) 
#        {result<-bibd.repo[[bibd.repo.table[i,]$id]]; break;}}
#  else {
#    for (i in seq_along(bibd.repo.table[,1]))
#      if (bibd.repo.table[i,]$v==v & bibd.repo.table[i,]$k==v-k) 
#        {result<-bibd.repo[[bibd.repo.table[i,]$id]]; break;}
#    if (!is.null(result)) result<-bibd(result,method=2)    
#  }      
#  return(result)
#}

################################################################
# METHOD 1 - TRIVIAL
################################################################

bibd.method1<-function(v, k)
{
  return(t(combn(v,k)))
}

################################################################
# METHOD 2 - DUAL DESIGN
################################################################

bibd.method2<-function(v, k, method=0)
{
# todo: switch
  A <- bibd(v,v-k,method)
# or even bibd ... with method=0
  
  abc<-unique(as.vector(A))
  B<-NULL
  for (i in 1:nrow(A))
    B<-rbind(B,setdiff(abc,A[i,]))
  return(B)
}

################################################################
# METHOD 6 - LATIN SQUARES
################################################################

#require(crossdes)
bibd.method6<-function(v, k, Nsim=1, seed=NA)
{
  if (!is.na(seed)) set.seed(seed);
  # get prime factorization
  f <- as.numeric(as.character(factorize(v)))
  fu <- unique(f)
  if(length(fu)==1){ #is it a prime power?
    p <- as.integer(fu)
    m <- length(f)

  A<-MOLS(p,m)
  flag<-TRUE
  solution<-matrix(0,1,1)

  for (index in 1:Nsim)
    {
      B<-matrix(0,k,dim(A)[3]*dim(A)[2])
      elmnts<-sample(p^m,k)
      for (i in 1:dim(A)[3])
        for (j in 1:dim(A)[2])
          B[,(i-1)*dim(A)[2]+j]<-sort(match(elmnts,A[,j,i]))
      B<-unique(t(B))
      if (flag | nrow(B)<nrow(solution)) solution<-B
      flag<-FALSE
    }
  
  return(solution)
  } else {
    stop("method 6 not applicable")
  }
}

################################################################
# METHOD 5 - TWO DESIGNS TOGETHER
################################################################

bibd.method5<-function(A)
{
  l<-(nrow(A)-3)/4
  if (l!=round(l) | length(unique(as.vector(A)))!=4*l+3 | ncol(A)!=2*l+1) stop("Assumptions violated!")
  
  B<-bibd.method2(A) 
  A<-cbind(A,4*l+4)
  
  return(rbind(A,B))
}

################################################################
# METHOD 16 - HADAMARD MATRIX
################################################################

#require(survey)
bibd.method16<-function(v, k)
{
  if ((v+1) %% 4 != 0) stop("v+1 must be divisible by 4")
  if(k!=(v-1)/2) stop(" ... wrong k ....")
  A<-hadamard.matrix(v+1)[-1,-1]

  B<-NULL
  for (i in 1:nrow(A))
    B<-rbind(B,which(A[i,]==1))
  return(B)
}

#################################################################
# METHOD 7&8 - DIFFERENCE SETS
#################################################################

dset.table<-
  rbind(c(7,3,1,1,1,7),
        c(9,3,3,3,4,3),
        c(9,4,3,1,2,9),
        c(10,3,4,2,6,5),
        c(10,4,4,2,3,5),
        c(11,3,3,1,5,11),
        c(11,4,6,1,5,11),
        c(11,5,2,1,1,11),
        c(13,3,1,1,2,13),
        c(13,4,1,1,4,13),
        c(13,4,4,1,3,13),
        c(13,6,5,1,2,13),
        c(15,3,3,3,7,5),
        c(15,7,3,1,1,15),
        c(16,3,2,1,5,16),
        c(16,5,4,1,3,16),
        c(17,4,3,1,4,17),
        c(17,5,5,1,4,17),
        c(17,8,7,1,2,17),
        c(19,3,1,1,3,19),
        c(19,4,2,1,3,19),
        c(19,6,5,1,3,19),
        c(19,9,4,1,1,19),
        c(21,3,3,3,10,7),
        c(21,4,3,1,5,21),
        c(21,5,1,1,5,21),
        c(21,6,9,3,6,7),
        c(22,4,4,2,7,11),
        c(22,7,8,2,4,11),
        c(23,11,1,1,11,23),
        c(25,3,1,1,4,25))
colnames(dset.table)<-c("v","k","rho","s","t","m")
dset.table<-as.data.frame(dset.table)

dsets<-NULL
dsets[[1]]<-matrix(c(1,2,4),1,3,byrow=T)
dsets[[2]]<-matrix(c(1,2,0+1i,1+1i,2+1i,0+2i,1+2i,2+2i,0,0,0+1i,0+2i),4,3,byrow=T)
dsets[[3]]<-matrix(c(0,1,2,4,0,3,4,7),2,4,byrow=T)
dsets[[4]]<-matrix(c(0,3,1+1i,1,2,1+1i,1+1i,4+1i,4,0+1i,4+1i,1+1i,0,3,4+1i,1,2,4+1i),6,3,byrow=T)
dsets[[5]]<-matrix(c(0,3,0+1i,4+1i,0,0+1i,1+1i,3+1i,1,2,3,0+1i),3,4,byrow=T)            
dsets[[6]]<-matrix(c(0,1,3,0,1,5,0,2,7,0,1,8,0,3,5),5,3,byrow=T)
dsets[[7]]<-matrix(c(0,1,8,9,0,2,5,7,0,1,4,5,0,2,3,5,0,4,5,9),5,4,byrow=T)
dsets[[8]]<-matrix(c(1,3,4,5,9),1,5,byrow=T)
dsets[[9]]<-matrix(c(0,1,4,0,2,7),2,3,byrow=T)
dsets[[10]]<-matrix(c(1,3,9,13),1,4,byrow=T)
dsets[[11]]<-matrix(c(0,1,2,4,8,0,1,3,6,12,0,2,5,6,10),3,5,byrow=T)
dsets[[12]]<-matrix(c(0,1,3,6,7,11,0,1,2,3,7,11),2,6,byrow=T)
dsets[[13]]<-matrix(c(1,4,0+1i,2,3,0+1i,1,4+1i,0+2i,2+1i,3+1i,0+2i,1+2i,4+2i,0,2+2i,3+2i,0,0,0+1i,0+2i),7,3,byrow=T)                                 
dsets[[14]]<-matrix(c(1,2,4,5,8,10,15),1,7,byrow=T)
dsets[[15]]<-matrix(c(0,1,3,0,4,9,0,2,8,0,3,7,0,1,6),5,3,byrow=T)
dsets[[16]]<-matrix(c(0,1,2,4,7,0,1,5,8,10,0,1,3,7,11),3,5,byrow=T)
dsets[[17]]<-matrix(c(0,1,3,7,0,1,5,7,0,1,6,9,0,2,5,9),4,4,byrow=T)
dsets[[18]]<-matrix(c(0,1,4,13,16,0,3,5,12,14,0,2,8,9,15,0,6,7,10,11),4,5,byrow=T)
dsets[[19]]<-matrix(c(0,1,3,7,8,12,14,15,0,2,3,4,7,8,9,11),2,8,byrow=T)

dsets.table2<-rbind(c(6,3,2,1,1,1,5),
                    c(8,4,3,1,1,1,7),
                    c(10,5,4,1,1,1,7),
                    c(12,3,2,1,3,1,11),
                    c(12,6,5,1,1,1,11),
                    c(14,7,6,1,1,1,13),
                    c(15,5,4,1,2,1,14),
                    c(15,6,10,2,3,2,7),
                    c(16,8,7,1,1,1,15))

dsets2<-NULL
dsets2[[1]]<-matrix(c(0,1,3,NaN,0,1),2,3,byrow=T)
dsets2[[2]]<-matrix(c(0,1,2,4,NaN,1,2,4),2,4,byrow=T)
dsets2[[3]]<-matrix(c(0,1,2,4,8,NaN,0,1,4,6),2,5,byrow=T)
dsets2[[4]]<-matrix(c(0,1,3,0,1,5,0,2,5,NaN,0,4),4,3,byrow=T)
dsets2[[5]]<-matrix(c(0,1,3,7,0,1,3,5,NaN,0,1,5),3,4,byrow=T)
dsets2[[6]]<-matrix(c(0,1,3,4,5,9,NaN,0,4,5,6,8),2,6,byrow=T)
dsets2[[7]]<-matrix(c(0,1,3,4,9,10,12,NaN,1,3,4,9,10,12),2,7,byrow=T)
dsets2[[8]]<-matrix(c(0,1,4,9,11,0,1,4,10,12,NaN,0,1,2,7),3,5,byrow=T)
dsets2[[9]]<-matrix(c(1,2,4,0+1i,1+1i,3+1i,2,3,5,0+1i,1+1i,3+1i,0,4,5,0+1i,1+1i,3+1i,NaN,0,0+1i,1+1i,2+1i,4+1i,NaN,0,3,5,6,0+1i),5,6,byrow=T)
dsets2[[10]]<-matrix(c(3,4,5,6,8,10,11,14,NaN,0,1,2,7,9,12,13),2,8,byrow=T)


bibd.method7<-function(ds, m)
{
  result<-NULL
  for (i in 1:nrow(ds))
    for (j in 0:(m-1))
      result<-rbind(result,complex(real=(Re(ds[i,])+j)%%m,imaginary=Im(ds[i,])))  
  dr<-dim(result)
  result<-as.numeric(as.factor(result))
  dim(result)<-dr 
  for (i in 1:nrow(result)) result[i,]<-sort(result[i,])
  return(result) 
}

bibd.method8<-bibd.method7

admissible.triples <- function(v){
  U <- 0:(v-1)
  W <- c(0:(v-2), Inf)
  d <- matrix(0,nrow=v,ncol=v)
  for(u in U){
    for(w in W){
      if(is.infinite(w)){
        d[u+1,length(W)] <- Inf
      }
      else
        d[u+1,w+1] <- min((u-w)%%v,(w-u)%%v)
    }
  }
#        browser()
  # 
  ds <- unique(c(d))
  maxd <- max(ds,na.rm=T)
  choices.index <- factor.comb(length(ds),3)+1
  ret <- NULL
  collect <- NULL
  for(j in 1:dim(choices.index)[1]){
    triple <- c(ds[choices.index[j,1]],
                ds[choices.index[j,2]],
                ds[choices.index[j,3]])
    triple <- sort(triple)
    if(triple[1]+triple[2]==triple[3]){
      chr <- ""
      for(k in triple){
        chr <- paste(chr,k,sep="")
      }
      if(!(chr %in% collect)){
        ret <- rbind(ret,triple)
        collect <- cbind(collect,chr)
      }
    }
  }
  dimnames(ret) <- NULL
  ret
}

admissible.triples.repeated <- function(v, lambda){
  t <- admissible.triples(v)
  
}

not.yet.implemented <- function(){
  stop("not yet implemented")
}

bibd.method9 <- function(v, k){
  not.yet.implemented()
}
bibd.method10 <- function(v, k){
  not.yet.implemented()
}
bibd.method11 <- function(v, k){
  not.yet.implemented()
}
bibd.method12 <- function(v, k){
  not.yet.implemented()
}
bibd.method13 <- function(v, k){
  not.yet.implemented()
}
bibd.method14 <- function(v, k){
  not.yet.implemented()
}
bibd.method15 <- function(v, k){
  not.yet.implemented()
}
bibd.method17 <- function(v, k){
  not.yet.implemented()
}
bibd.method18 <- function(v, k){
  not.yet.implemented()
}



