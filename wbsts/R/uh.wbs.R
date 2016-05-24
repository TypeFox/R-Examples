uh.wbs <-
function(z,C_i, del=-1, epp, scale,M=0,cstar=0.75){
  
  l <- dim(z)[2]
  n<-nz<-dim(z)[1]
  estimates<-br.list<-NULL
  
  if(epp[1]==-1){
    epp<-c()
    for (j in 1:l) {
      ep = round(max(2*n/2^scale[j], ceiling(sqrt(n)/2)))
      epp = c(epp,ep)
    }
    epp=round(epp/2) ###Be careful here!!!!!!!
  }
  if(del<0){ 
    del<-round(max(2*epp[1]+2, ceiling(log(n)*sqrt(n)/4)))
  }
  ep=epp[1]
  if(nz<del) stop("Input vector too short")
  criterion <- C_i[1] * log(n)   
  breakpoints<-NULL
  f<-NULL
  tree<-list(matrix(0, 6, 1))
  tree[[1]][1,1]<-1
  s<-tree[[1]][4,1]<-1
  e<-tree[[1]][6,1]<-nz
  b=-1
  while(b < 0 | b > e) { 
    temp.r = cr.rand.max.inner.prod(XX=z[s:e,],Ts=n,C_i=C_i,epp=epp,M=M,cstar=cstar,Plot=0)
    b <- temp.r[[1]]
  }
  tree[[1]][5,1] <- b
  d<-temp.r[[2]]
  if(abs(d)>criterion){
    
    tree[[1]][2,1]<-d
    breakpoints <- c(breakpoints, b)
    f<-c(f, d)  
    j<-1
    while(length(tree)==j){
      if(sum(tree[[j]][6, ]-tree[[j]][4, ]-rep(2*del, dim(tree[[j]])[2]))>0){ 
        no.parent.coeffs<-dim(tree[[j]])[2]
        no.child.coeffs<-0
        for(i in 1:no.parent.coeffs){
          if(tree[[j]][5, i]-tree[[j]][4, i]>max(ep, del)+1){
            s<-tree[[j]][4, i]
            e<-tree[[j]][5, i]
            ind.max=-1
            while(ind.max < 0 | ind.max > e) { 
              temp.RIP = cr.rand.max.inner.prod(XX=z[s:e,],Ts=n,C_i=C_i,epp=epp,M=M,cstar=cstar)
              ind.max = temp.RIP[[1]]
            }
            b<-s+ind.max-1
            d<-temp.RIP[[2]]
            
            if(abs(d)>criterion){
              if(length(tree)==j) tree<-c(tree, list(matrix(0, 6, 0)))
              no.child.coeffs<-no.child.coeffs+1
              tree[[j+1]]<-matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, no.child.coeffs)
              tree[[j+1]][1, no.child.coeffs]<-2*tree[[j]][1, i]-1
              tree[[j+1]][2, no.child.coeffs]<-d
              tree[[j+1]][4, no.child.coeffs]<-s
              tree[[j+1]][6, no.child.coeffs]<-e
              tree[[j+1]][5, no.child.coeffs]<-b
              if (sum(abs(breakpoints-b)>del)==length(breakpoints))  breakpoints<-c(breakpoints,b);
              f<-c(f, d)
            }
          }
          if(tree[[j]][6, i]-tree[[j]][5, i]>max(ep, del)+1){
            s<-tree[[j]][5, i]+1
            e<-tree[[j]][6, i]
            ind.max=-1
            while(ind.max < 0 | ind.max > e) { 
              temp.RIP = cr.rand.max.inner.prod(XX=z[s:e,],Ts=n,C_i=C_i,epp=epp,M=M,cstar=cstar)
              ind.max = temp.RIP[[1]]
              
            }  
            b<-s+ind.max-1
            d<-temp.RIP[[2]] 
            
            if(abs(d)>criterion){
              if(length(tree)==j) tree<-c(tree, list(matrix(0, 6, 0)))
              no.child.coeffs<-no.child.coeffs+1
              tree[[j+1]]<-matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, no.child.coeffs)
              tree[[j+1]][1, no.child.coeffs]<-2*tree[[j]][1, i]
              tree[[j+1]][2, no.child.coeffs]<-d 
              tree[[j+1]][4, no.child.coeffs]<-s
              tree[[j+1]][6, no.child.coeffs]<-e
              tree[[j+1]][5, no.child.coeffs]<-b
              if (sum(abs(breakpoints-b)>del)==length(breakpoints))  breakpoints<-c(breakpoints,b)
              f<-c(f, d)
              
            }
          }
        }
        
      }
      j<-j+1
    }
  }
  list(tree=tree, breakpoints=breakpoints, f=f)  
}
