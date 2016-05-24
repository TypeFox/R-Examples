ReducedKM <- function(data,nclus,ndim,nstart=100,smartStart=FALSE){
  #Reduced k-means, minimize f(A,U,Y)=|| X - UYA' ||^2  
  
  data = data.matrix(data)
  n = dim(data)[1]
  m = dim(data)[2]
  conv=1e-6  # convergence criterion
  func={}; index={}; AA = {}; FF = {}; YY = {}; UU={}
  for (run in c(1:nstart)) {
    
    # Starting method
    if(smartStart==TRUE){ # k-means starting solution
      outkstart=kmeans(data,nclus,nstart=100)
      v=as.factor(outkstart$cluster)
      U = diag(nlevels(v))[v,] #dummy cluster membership
      #update A
      D = pseudoinverse(diag(apply(U,2,sum)))
      s = t(data)%*%U%*%D%*%t(U)%*%data
      #Sort and normalize eigenvectors of s
      u = ed(0.5*s+0.5*t(s))$u
      A=u[,c(1:ndim)]
    } else { #random start
      v = as.factor(ceiling(nclus*runif(n)))
      U = matrix(0,n,nclus) #  random start for U
      for (i in 1:n) {
        U[i,v[i]] =1
      }
      #random A
      A = orth(matrix(rnorm(m*ndim),m,ndim)) # random start for A
    } 
    
    #update Y
    su = apply(U,2,sum)
    ind = which(su > 0, arr.ind=TRUE)
    su[ind] = (su[ind])^(-1)
    D = diag(su)
    DU = D%*%t(U)
    Y = DU%*%data%*%A%*%pseudoinverse(t(A)%*%A)  
    
    a= data - U%*%Y%*%t(A)
    f = ssq(a) 
    fold = f+2 * conv*f
    iter = 0
    
    #iterative part
    while (f<fold-conv*f) {
      fold=f
      iter=iter+1
      # update partitioning U  (given A and Y)
      P = matrix(0,n,nclus)
      for (j in (1:nclus)) {
        Pj = (data%*%A-data.matrix(rep(1,n))%*%Y[j,])^2
        if (ndim>1) 
          P[,j] = t(apply(t(Pj),2,sum))  # squared distances to cluster j
        if (ndim==1)
          P[,j] = Pj	# squared distances to cluster j
      }
      
      mm = apply(t(P),2,min)
      v = as.factor(apply(t(P),2,which.min))
      U = matrix(0,n,nclus)
      for (i in 1:n) {
        U[i,v[i]] =1
      }    
      
      # correction in case zero columns emerge (developed by Rocci)
      # split one cluster into two by rosplit
      su=apply(U,2,sum)
      F = data %*% A
      while (sum(su==0)>0) {
        dw=t(U)%*%((F-U%*%pseudoinverse(t(U)%*%U)%*%t(U)%*%F)^2)%*%data.matrix(rep(1,ndim))
        
        mm = apply(data.matrix(su),2,min)
        p1 = apply(data.matrix(su),2,which.min)
        
        mm = apply(dw,2,max)
        p2 = apply(dw,2,which.max)
        
        ind = which(U[,p2]!=0, arr.ind=TRUE)
        ind1=ind[1:floor(su[p2]/2)]
        U[ind1,p1] = 1
        U[ind1,p2] = 0
        U1=rosplit(F[ind,],U[ind,c(p1,p2)])$U
        U[ind,c(p1,p2)] = U1
        su=apply(U,2,sum)
      }
      
      # update A
      D = pseudoinverse(diag(apply(U,2,sum)))
      s = t(data)%*%U%*%D%*%t(U)%*%data
      #Sort and normalize eigenvectors of s
      u = ed(0.5*s+0.5*t(s))$u
      A=u[,c(1:ndim)]
      
      #update Y
      su = apply(U,2,sum)
      ind = which(su > 0, arr.ind=TRUE)
      su[ind] = (su[ind])^(-1)
      D = diag(su)
      DU = D%*%t(U)
      Y = DU%*%data%*%A
      
      # criterion
      a = data - U%*%Y%*%t(A)
      f = ssq(a) #sum(sum(a^2))
      #fpX=1-sum(sum(a^2)) /sum(sum((data.matrix(data))))
      
    }
  
  func[run] = f
  #fpXunc[run]=fpX
  FF[[run]] = F
  AA[[run]] = A
  YY[[run]] = Y
  UU[[run]] = U
}
  #assign output
  out=list()
  mi = which.min(func)
  out$obscoord=FF[[mi]]
  rownames(out$obscoord) = rownames(data)
  out$attcoord=as.matrix(AA[[mi]])
  rownames(out$attcoord) = colnames(data)
  out$centroid=YY[[mi]]
  U=UU[[mi]]
  out$cluID=apply(U,1,which.max)
  out$criterion=func[mi]
  out  
}

ssq = function(a) {
  t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
}

