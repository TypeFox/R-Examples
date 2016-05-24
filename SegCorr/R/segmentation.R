###############################################################################
#########################            segmentation                   ###########
################################################################################
### Summary: performing block diagonal segmentation and then the exact test to detect the highly 
###          correlated regions

### Input: CHR:             chromosome name
###        EXP:             gene expression
###        genes:           name of the genes


### Output: Results:       a matrix containing the highly correlated regions in descending order (the most significant on top).
###                        CHR: the chromosome of the region. Start/End: are the order of the gene within each chromosome.
###                        Rho: is the correlation.  length: number of genes in the region.
###                        first/last gene: the name of the first/last gene in the region
###                        pvalue: the pvalue as obtained from the test. genes: the names of the genes belonging to the region
###         rho0:          background correlation
###         likelihood:    the log-likelihood of the model
###         K:             number of segments
#options(warn=1)

segmentation  = function(CHR,EXP,genes,S){
  
  if(missing(CHR)){
    warning('Chromosome name missing') 
    CHR = NA
  }else{
    cat(paste('Performing block diagonal segmentation for chromosome =',CHR,sep=" "),sep="\n")
  }
  
  if(missing(EXP)){
    
    stop("segmentation: Expression matrix missing")
    
    
  }else{
    
    
    p = dim(EXP)[1]
    n = dim(EXP)[2]

    
    
    if(missing(genes)){
      warning('segmentation: Gene names missing. Using artificial names')
      genes = 1:p
    }else{
      if(length(genes)!=p){
        stop("segmentation: Dimension mismatch. Expression matrix and gene vector")
      }
    }
    
    EXP = t(EXP)
    EXP = sqrt(n/(n-1))*scale(EXP, center = T, scale = T)
    G = cor(EXP)
    
    if(missing(S)){
      S=0.7
      warning("segmentation: Using default value for S")
    }
      
    #######################################
    ####### Define Kmax ###################
    ########################################
    
    Kmax =  floor(p/5)
    
    #########################################
    #####     Segmentation   ################
    #########################################
    
    
    DP=  .F_SegKern_Norm(G, Kmax,n)
    
    ########################################
    ######## Estimate gamma.0 #############
    ########################################
    
    Gamma.0 = max(c(0, median(diag(G[-1,])) ) )
    #############################################################################
    ######################## Model Selection Criterion ##########################
    #############################################################################
    #S = 0.7
    loglik=(-1/2)*(DP$J.est+n*p*log(2*pi))
    loglik = loglik[1:Kmax]
    Kseq=1:Kmax
    
    pen=5*Kseq+2*Kseq*log(p/Kseq)
    
    K.PenJ=.selection1(loglik,pen,S)
    
    #############################################################################
    tau.hat = c(0,DP$t.est[K.PenJ,1:K.PenJ],p)
    
    
    Results.rho = matrix(NA,K.PenJ,9)
    Results.rho = data.frame(Results.rho)
    colnames(Results.rho) =c('CHR','Start','End','Rho',"length","first gene","last gene","p.value","genes")
    for(i in 1:K.PenJ){
      
      l = length((tau.hat[i]+1):tau.hat[i+1])
      B = (sum(G[(tau.hat[i]+1):tau.hat[i+1],(tau.hat[i]+1):tau.hat[i+1]]) -  l)/(l^2-l)
      
      
      
      Results.rho[i,1] = CHR
      Results.rho[i,2] = tau.hat[i]+1
      Results.rho[i,3] = tau.hat[i+1]
      Results.rho[i,4] = B
      Results.rho[i,5] = l
      Results.rho[i,6] = genes[Results.rho[i,2]]
      Results.rho[i,7] = genes[Results.rho[i,3]]
      ####################################################### 
      ####### calculating pvalues ###########################
      #######################################################
      data.loc = EXP[,Results.rho[i,2]:Results.rho[i,3]]
      mean.loc = apply(data.loc,1,mean)
      mm = (t(data.loc)%*%mean.loc)/n
      mp= mean(mm)
      khi.loc.gamma0 = mean(mean.loc*mean.loc)*(n*dim(data.loc)[2])/(1+(dim(data.loc)[2]-1)*Gamma.0)
      Results.rho[i,8] = pchisq(khi.loc.gamma0, n-1, ncp=0, lower.tail = F, log.p = F) 
      Results.rho[i,9] = paste(genes[Results.rho[i,2]:Results.rho[i,3]],collapse=",")
    }
    
    
    list(Results=Results.rho,rho0 = Gamma.0, likelihood = loglik[K.PenJ], K = K.PenJ)   
    
  }
  
}


##############################################################################################################
#################################  penalization criterion Marc    ############################################
##############################################################################################################
### Summary: performing the model selection criterion by Lavielle

### Input: Lv:              log-likelihood for different no of segments (Kseq)
###        Kseq:            penalty function
###        S:               threshold. default S=.7


.selection1<-function(Lv, Kseq,S)
{
  Kmax=length(Kseq)
  J=-Lv
  Jtild=(J[Kmax]-J)/(J[Kmax]-J[1])*(Kseq[Kmax]-Kseq[1])+1
  D=diff(diff(Jtild))
  #plot(D)
  if (length(which(D>=S))>0){  
    rg <- max(which(D>=S))+1
  }else {
    rg <- 1
  }
  rg=Kseq[rg]
  Kh=which(rg==Kseq)
  return(Kh)
}


##############################################################################################################
#################################  DP algortihm                  ############################################
##############################################################################################################

### Summary: computes the change points given a cost matrix matD and a maximum number of segments Kmax.
###          J.est[K]    = minimum contrast for a model with K segments (-J.est is the log-likelihood)
###          t.test[K,:] = coordinates of the change points for a model with K segments

### Input: matD:            cost matrix
###        Kmax:            maximum number of segments to be considered

.DynProg<-function(matD,Kmax){
  
  N<-dim(matD)[1]
  if (Kmax>N){
    cat("Kmax ", Kmax, "is greater than N ",N,"\n")
    cat("Kmax is supposed to be equal to N :", N,"\n")
  }  
  
  I<-matrix(Inf,Kmax,N)
  t<-matrix(0,Kmax,N)   
  I[1,]=matD[1,]
  matD=t(matD)
  
  if (Kmax>2) {
    
    for (k in 2:(Kmax-1)){      
      for (L in k:N){                
        I[k,L]<-min(I[(k-1),1:(L-1)]+matD[L,2:L])
        
        if(I[k,L]!=Inf){
          t[(k-1),L]<-which.min(I[(k-1),1:(L-1)]+matD[L,2:L])
        } else {
          t[(k-1),L] = Inf
        } # end else
      } #end L      
    }#end k
  } #end K
  
  
  I[Kmax,N]<-min(I[(Kmax-1),1:(N-1)]+matD[N,2:N])
  
  if(I[Kmax,N]!=Inf){
    t[(Kmax-1),N]<-which.min(I[(Kmax-1),1:(N-1)]+matD[N,2:N])
  } #end if
  
  
  # *** computation of breakpoint instants ***
  
  t.est<-matrix(0,Kmax,Kmax)
  diag(t.est)<-N
  for (K in 2:Kmax){
    for (k in seq(K-1,1,by=-1)){
      if(t.est[K,k+1]!=0){            
        t.est[K,k]<-t[k,t.est[K,k+1]]
      } #end if
    } #end k
  } #end K
  t.est[which(is.na(t.est))]=Inf
  list(J.est = I[,N],t.est = t.est)
  #cat("Kmax", Kmax,"\n")
}
##############################################################################################################
#################################  calculate the cost matrix     ############################################
##############################################################################################################

### Summary: calculates the cost matrix for the DP.


### Input: H:               correlation matrix
###        Kmax:            maximum number of segments to be considered
###        p:               number of patients

.F_SegKern_Norm <- function(H, Kmax,p)
{
  options(warn=-1)
  cat('Block diagonal segmentation:')
  CumDiagH = cumsum(diag(H))
  n = length(CumDiagH)
  
  # Cumulated sum within rectangles
  cat(' rectangles,')
  # R[i, j] = sum_{1 <= u <= i} sum_{1 <= v <= j} H[i, j]
  R = apply(H, 2, cumsum)
  R = t(apply(R, 1, cumsum))
  
  # Cumulated sums within diagonal blocks
  cat(' blocks,')
  # B[i, j] = sum_{i <= u <= j} sum_{u <= v <= j} H[i, j]
  RR = matrix(0, (n+1), (n+1))
  RR[(2:(n+1)), (2:(n+1))] = R
  diagRR = diag(RR)
  RRjj = matrix(rep(diagRR, each=(n+1)), (n+1), (n+1))
  RRii = matrix(rep(diagRR, (n+1)), (n+1), (n+1))
  BB = RRjj + RRii - RR - t(RR)
  BB[1, ] = c(0, R[1, ])
  B = BB[1:n, 2:(n+1)]
  
  # Next line : if correction for symmetry is needed
  B[lower.tri(B)] = 0
  
  B = B + t(B) - diag(diag(B))
  
  ## introduced these lines because line 1 was not sum(H[i:j,i:j])
  B[1,] = diag(R)
  B[,1]= diag(R)
  
  # Segment length
  Lg = matrix(rep((1:n), each=n), n, n) - matrix(rep(0:(n-1), n), n, n)   
  
  # Cost matrix
  cat(' costs,')
  
  L = Lg + (Lg-1)*log(-(B-Lg*Lg)/((Lg-1)*Lg)) + log(B/Lg) ## on log scale
  L = p*L
  diag(L) = Inf #n
  L[lower.tri(L)] = Inf
  # restriction of having segments of length two
  L[seq(from = n+1,to = (n+1)*(n-1),by = n+1 )] = Inf
  
  # Segmentation
  cat(' segmentation',sep="\n")
  DP = .DynProg(L, Kmax)
  J.est = DP$J.est; t.est = DP$t.est
  
  return(list(J.est=J.est, t.est=t.est))
}

