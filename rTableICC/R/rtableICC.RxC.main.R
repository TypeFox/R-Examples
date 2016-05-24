rtableICC.RxC.main <-
function(p=NULL,theta,M,row.margins=NULL,col.margins=NULL,sampling="Multinomial",N=1,lambda=NULL,zero.clusters=FALSE,print.regular,print.raw){

  
  T=length(theta)
   
  if (sampling=="Product"){ 
    if ((is.null(p)==TRUE) | (length(dim(p))!=2) | (min(p)<0) | (max(p)>1) | (round(sum(p))!=1)){
      stop("Dimension of matrix (p>0) must be RxC and it must sum up to one!")
    } else {
      R=nrow(p)
      C=ncol(p)  
      if (min(R,C)<=1){
        stop("Minimum colum or row length must be greater than 1!")
      } 
    }
    if ((is.null(col.margins)==TRUE) & (is.null(row.margins)==TRUE)){
      stop("Row or columun margins must be entered as a scalar vector!")
    } else if ((length(row.margins)>1) & (length(col.margins)>1)){
      stop("Number of either row or columun margins must be greater than one. Both cannot be greater than one!")
    } else if ((length(row.margins)<2) & (length(col.margins)<2)){      
      stop("Number of fixed row or columun margins must be greater than one under product multinomial sampling!")
    }
    rTable=array(0,dim=c(1,R*C)) 
    g.t=array(0,dim=c(R,C,(T-1)))
    g.tilde=array(0,dim=max((T-1)))
    
    if (length(row.margins)>1){   
      rTable.raw=array(0,dim=c(R,C,sum(row.margins))) 
      row.margins=abs(round(row.margins))
      N=sum(row.margins)
      prob.margins=row.margins/N
      if ((apply(p,1,sum)!=prob.margins)){
        stop("Mismatch between cell probabilities and margin probabilities!")
      }            
     
      if (length(M)==1){
        M=rep(M,times=length(row.margins))
      }
      
      say=1
      say2=1
      for (row.index in 1:length(row.margins)){
        cluster.size=0
        if (zero.clusters==TRUE){
          cluster.size=rmultinom(1, row.margins[row.index], rep(1/M[row.index],M[row.index]))
        }else if (zero.clusters==FALSE){
          if ((row.margins[row.index]-M[row.index])<0){
            stop("Because number of individuals is less than the total number of clusters, it is impossible to allocate an individual to each cluster! Set zero.clusters = TRUE and re-run the function.")
          }
          cluster.size=rmultinom(1, (row.margins[row.index]-M[row.index]), rep(1/M[row.index],M[row.index]))+1
        }
        gen=rtableICC.RxC.engine(1,C,T,M[row.index],p[row.index,]/(row.margins[row.index]/N),row.margins[row.index],cluster.size,theta)
        
        rTable[say2:(say2+C-1)]=gen$rTable
        say2=C+1
        rTable.raw[row.index,1:C,say:(say+row.margins[row.index]-1)]=gen$rTable.raw
        g.t[row.index,1:C,]=gen$g.t
        say=row.margins[row.index]+1
        g.tilde=g.tilde+gen$g.tilde
      }
      if (print.regular==TRUE){
        rTable.regular=array(0,dim=c(R,C))
        say=0
        for (i in 1:R){      
          for (j in 1:C){
            say=say+1
            rTable.regular[i,j]=rTable[say]      
          }
        }
        list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,rTable.regular=rTable.regular,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
             print.raw=print.raw,print.regular=print.regular)
      } else {
        list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
             print.raw=print.raw,print.regular=print.regular)
      }
    }else if (length(col.margins)>1){
      rTable.raw=array(0,dim=c(R,C,sum(col.margins))) 
      col.margins=abs(round(col.margins))
      N=sum(col.margins)
      prob.margins=col.margins/N
      if ((apply(p,2,sum)!=prob.margins)){
        stop("Mismatch between cell probabilities and margin probabilities!")
      }
      
      if (length(M)==1){
        M=rep(M,times=length(col.margins))
      }
      
      say=1
      say2=1
      for (col.index in 1:length(col.margins)){
        cluster.size=0
        if (zero.clusters==TRUE){
          cluster.size=rmultinom(1, col.margins[col.index], rep(1/M[col.index],M[col.index]))
        }else if (zero.clusters==FALSE){
          if ((col.margins[col.index]-M[col.index])<0){
            stop("Because number of individuals is less than the total number of clusters, it is impossible to allocate an individual to each cluster! Set zero.clusters = TRUE and re-run the function.")
          }          
          cluster.size=rmultinom(1, (col.margins[col.index]-M[col.index]), rep(1/M[col.index],M[col.index]))+1
        }
        gen=rtableICC.RxC.engine(R,1,T,M[col.index],p[,col.index]/(col.margins[col.index]/N),col.margins[col.index],cluster.size,theta)
        rTable[say2:(say2+R-1)]=gen$rTable
        say2=R+1
        rTable.raw[1:R,col.index,say:(say+col.margins[col.index]-1)]=gen$rTable.raw
        g.t[1:R,col.index,]=gen$g.t
        say=col.margins[col.index]+1
        g.tilde=g.tilde+gen$g.tilde
      }
      
      if (print.regular==TRUE){
        rTable.regular=array(0,dim=c(R,C))
        say=0
        for (i in 1:C){      
          for (j in 1:R){
            say=say+1
            rTable.regular[j,i]=rTable[say]      
          }
        }
        list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,rTable.regular=rTable.regular,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
             print.raw=print.raw,print.regular=print.regular)
      } else { 
        list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
           print.raw=print.raw,print.regular=print.regular)
      }
    }   
    
  } else if (sampling=="Multinomial"){
    if ((is.null(p)==TRUE) | (length(dim(p))!=2) | (min(p)<0) | (max(p)>1) | (round(sum(p))!=1)){
      stop("Dimension of matrix of cell probabilities (p>0) must be RxC and it must sum up to one!")
    } else {
      R=nrow(p)
      C=ncol(p)  
      if (min(R,C)<=1){
        stop("Minimum colum or row length must be greater than 1!")
      }
    }    
    if ((length(N)!=1) | (is.finite(N)==FALSE)){ 
      stop("Total number of observation should be entered as a scalar under multinomial samlping plan.")      
    }else{
      N=abs(round(N))      
      
      cluster.size=0
      if (zero.clusters==TRUE){
        cluster.size=rmultinom(1, N, rep(1/M,M))
      }else if (zero.clusters==FALSE){
        if ((N-M)<0){
          stop("Because number of individuals is less than the total number of clusters, it is impossible to allocate an individual to each cluster! Set zero.clusters = TRUE and re-run the function.")
        }
        cluster.size=rmultinom(1, (N-M), rep(1/M,M))+1
      }
      
      gen=rtableICC.RxC.engine(R,C,T,M,p,N,cluster.size,theta)
      rTable=gen$rTable
      rTable.raw=gen$rTable.raw
      g.t=gen$g.t
      g.tilde=gen$g.tilde
    }
    if (print.regular==TRUE){
      rTable.regular=array(0,dim=c(R,C))
      say=0
      for (i in 1:R){      
        for (j in 1:C){
          say=say+1
          rTable.regular[i,j]=rTable[say]      
        }
      }
      list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,rTable.regular=rTable.regular,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
           print.raw=print.raw,print.regular=print.regular)
    } else {
       list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
            print.raw=print.raw,print.regular=print.regular)
    }
  }else if (sampling=="Poisson"){
    R=nrow(lambda)
    C=ncol(lambda)  
    if (min(R,C)<=1){
      stop("Minimum colum or row length must be greater than 1!")
    }
    cell.counts=rpois(R*C,t(lambda))
    N=sum(cell.counts)
    p=cell.counts/N    
    p=matrix(p,nrow = R,ncol = C)
    cluster.size=0
    if (zero.clusters==TRUE){
       cluster.size=rmultinom(1, N, rep(1/M,M))
    }else if (zero.clusters==FALSE){
      if ((N-M)<0){
        stop("Because number of individuals is less than the total number of clusters, it is impossible to allocate an individual to each cluster! Set zero.clusters = TRUE and re-run the function.")
      }
      cluster.size=rmultinom(1, (N-M), rep(1/M,M))+1
    }

    gen=rtableICC.RxC.engine(R,C,T,M,p,N,cluster.size,theta)
    rTable=gen$rTable
    rTable.raw=gen$rTable.raw
    g.t=gen$g.t
    g.tilde=gen$g.tilde   
    if (print.regular==TRUE){
      rTable.regular=array(0,dim=c(R,C))
      say=0
      for (i in 1:R){      
        for (j in 1:C){
          say=say+1
          rTable.regular[i,j]=rTable[say]      
        }
      }
      list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,rTable.regular=rTable.regular,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
           print.raw=print.raw,print.regular=print.regular)
    } else{
        list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,N=N,cluster.size=cluster.size,sampling=sampling,M=M,R=R,C=C,T=T,ICC=TRUE,structure="RxC",
             print.raw=print.raw,print.regular=print.regular)    
    }  
  }  
  
  
  }
