rsm <-
function(X, sub, Klist, nbredo=5, disp=F, maxit=50){
  
  if(!is.matrix(X)) X= as.matrix(X)
  
  if(dim(X)[1] != dim(X)[2]) 
    stop("The adjacency matrix must be a square matrix.\n")
  
  if(length(sub) != dim(X)[1])
    stop("Each vertex must be in one subgraph.\n")
  
  if(length(Klist) == 1)
    if(Klist == 1)
      stop("Klist must be at least 2.\n")
  
  R = max(sub); N = nrow(X); C = max(X)
  a = matrix(as.numeric(X != 0),nrow=N,byrow=F)
  
  g=array(1,dim=c(R,R,2))
  g0=g
  
  totreg=matrix(0,R,R)
  areg=matrix(0,R,R)
  for (r in 1:R){
    for (s in 1:R){
      areg[r,s]=sum(a[which(sub==r),which(sub==s)])
      totreg[r,s]=length(which(sub==r))*length(which(sub==s))
    }
    
    totreg[r,r]=totreg[r,r]-length(which(sub==r))
  }
  
  gama = matrix(NA, ncol=R, nrow= R)
  for (r in 1:R){
    for (s in 1:R){
      g[r,s,1]=g0[r,s,1]+areg[r,s]
      g[r,s,2]=g0[r,s,2]+totreg[r,s]-areg[r,s]
      gama[r,s] = g[r,s,1]/sum(g[r,s,])
    }
  }  
  
  
  output <- list()

  for(K in Klist){
    if (disp) cat('* K=',K,':',sep='')
    
    lower.max <- -Inf
    
    output[[K]] <- list()
    
    chi=matrix(1,R,K)
    chi0=chi
    Tau=matrix(0,N,K)
    ppi=array(1,dim=c(K,K,C))
    ppi0=ppi      
    
    for(redo in 1:nbredo){
      init=.rsm.initialbis(X,sub,K)
      initbis=init
      
      savetab=sub
      for(r0 in 1:R){
        
        Tau=initbis[r0,,]
        zinit=Tau
        chi=matrix(1,R,K)
        chiold=chi+1
        ppi=array(1,dim=c(K,K,C))
        ppiold=ppi+1
        
        tmp = .rsm.Mstep(X,chi,ppi,Tau,sub,g)
        lower=tmp$lower
        chi=tmp$chi
        ppi=tmp$ppi
        
        
        former=lower-100
        
        itercnt=0
        
        condition=TRUE
        while(condition & itercnt<maxit){
          if (disp) cat('x')
          itercnt=itercnt+1
          
          if(former>lower+10^-10) cat("the likelihood decreases to",lower-former,"\n")
          
          chiold=chi
          ppiold=ppi
          former=lower
          Tau_temp=.rsm.Estep(X,chi,ppi,Tau,sub,g)
          tmp=.rsm.Mstep(X,chi,ppi,Tau,sub,g)
          chi=tmp$chi
          ppi=tmp$ppi
          lower=tmp$lower
          
          prec=max(c((sum(abs(chi-chiold))+sum(abs(ppi-ppiold))),lower-former))
          condition=(prec>10^-2)
          
          if (lower > lower.max){
            lower.max = lower
            chi.max = chi
            ppi.max = ppi                         
            z.max = Tau
          }
          
        }
        
        if (disp) cat('\n')
        if (disp) cat('Criterion value:',lower,'\n')
      }
    } 
    alpha_temp = alpha = matrix(NA, nrow= R, ncol=K)
    for(r in 1:R) alpha_temp[r, ] =  chi.max[r,]/ sum(chi.max[r,])

    Pi_temp = Pi = array(NA, dim=c(K, K, C))
    for(k in 1:K){
      for(l in 1:K)
        Pi_temp[k, l, ] = ppi.max[k, l, ]/ sum(ppi.max[k, l, ]) 
    }
    
    Zcol_temp = Z_col = seq(0,0,length=N)
    for (n in 1:N) Zcol_temp[n]=which.max(z.max[n,])
    
    output[[K]]$lower <- lower.max
    output[[K]]$alpha <- alpha_temp
    output[[K]]$Pi <- Pi_temp
    output[[K]]$Tau <- z.max  
    output[[K]]$Zcol <- Zcol_temp
  } 

  vec <- rep(-Inf, max(Klist))
  
  for(K in Klist)
    vec[K] <- output[[K]]$lower
  K_star = which.max(vec)
  
  result = list(N =N, R=R, C=C, Klist = Klist, X=X, K_star = K_star, gama = gama, output = output) 
  
  class(result) <- "rsm"
  
  return(result)
}
