TDDclust <- function(data,numClust,lambda,Th,niter,T0,simAnn,alpha,data1,verbose=TRUE){ 

  B <- 5 ; 
  Km <- 2 ; 
  norm <- 0 ; 
  acc <- c() 
  improv_kl <- list() ; 
  indiv_trimmed <- list() ; 
  NN1_aux <- list()
  if(verbose){
   cat("Step (1): \n Initial partition with PAM:")
  }
  
  pcc <- pam(data,numClust)
  if(verbose){
   print(table(pcc$clustering))
   cat("Medoids:")
   print(rownames(data[pcc$id.med,]))   
  }
  Y <- pcc$med   
  
  prs <- 0 
  move <- 0 
  
  if(verbose){
   cat("Step (2): \n")     
  }
  NN <- NNDDVQA1(data,pcc,lambda,norm) 
  NNa <- NN 
  NNold_old <- NN
    
  Cost <- NNa$Kmata   
  if(verbose){
   cat("Objective function value for this partition:")      
   print(mean(Cost))
  }
  Ls <- sum(NN$Kmata) 
  
  
  qa <- 0 
  
  if(niter > 1){#niter - number of iterations, default 5 
    if(verbose){
     cat("Iteration for niter")   
     print(niter)
    }
    kl <- 0 
    
    while(kl < (niter-1)){ 
      if(verbose){
       cat("Iteration for kl:") 
       print(kl)
      }
      kl <- kl + 1 
       
      if(niter<=5 & kl>=(niter-5)){ 
        if(verbose){
         cat("\n")  
         cat("Steps (3),(4),(5),(6) (for niter <= 5 | kl>=(niter-5)) NNDDVQE: \n") 
        }
        NN <- NNDDVQE(data,Y,0,lambda,NN,Km,Th,alpha,kl,NNold_old,verbose)
        
        if(sum(table(NN$NN[1,])) == dim(data)[1]){
         if(verbose){
          cat("The optimal partition is provided by PAM")
         }
         Y <- as.integer(rownames(Y))
         return(list(NN=NN$NN,cases=Y,Cost=NN$Cost,discarded="None",klBest=0))
        }
        
      }
      
      if(niter>5 & kl!=0){ 
        if(verbose){
         cat("\n")  
         cat("Steps (3),(4),(5),(6) (for niter > 5 & kl!=0) NNDDVQEstart: \n")
        }
        NN <- NNDDVQEstart(data,Y,T0,lambda,NN,Km,Th,alpha,kl,NNold_old,verbose) 
        
        if(sum(table(NN$NN[1,])) == dim(data)[1]){
          if(verbose){
           cat("The optimal partition is provided by PAM")
          }
          Y <- as.integer(rownames(Y))
          return(list(NN=NN$NN,cases=Y,Cost=NN$Cost,discarded="None",klBest=0))
        }
        
        T0 <-NN$T0 * simAnn 
        Th <- Th * simAnn 
      }
      
      aux1 <- list()    ; Kmata_aux <- list() ; pcc_aux1 <- list() 
      aux_mat <- list() ; DDi_aux <- list()   ; DD_aux <- list()
      
      NNold_old <- NN 
      if(verbose){
       cat("Optimal partition when finishing the kl iteration: \n") 
       print(table(NNold_old$Nuvec))
       cat("\n") 
      }
      NN1_aux[[kl]] <- NN
      NNN <- NN
      
      if(length(NN$Nuvec) < dim(data1)[1]){
        indiv_trimmed[[kl]] <- NN$trimmed
        improv_kl[[kl]] <- NN$improv
        
        for(trimm in NN$trimmed){
          for(nu in 1:numClust){
            aux1[[nu]] <- append(NN$Nuvec,nu,trimm-1)
            pcc_aux1[[nu]] <- pamsil(data1,aux1[[nu]],numClust)
            
            DDi_aux[[nu]] <- matrix(0,dim(data1)[1],numClust) 
            Nvec <- rep(1,dim(data1)[1]) 
            for(kz in (1:numClust)){ 
              for(ky in (1:dim(data1)[1])){ 
                Xmatr <- data1[aux1[[nu]] == kz,] 
                Nvecr <- Nvec[aux1[[nu]] == kz]
                DDi_aux[[nu]][ky,kz] <- DDfcnadj(Xmatr,Nvecr,data1[ky,]) 
              } 
            }
            
            num <- sample((1:numClust)[-nu],1)
            aux2 <- append(NN$NN[2,],num,trimm-1)
            aux3 <- append(NN$NN[3,],(1:numClust)[-c(nu,num)],trimm-1)
            aux_mat[[nu]] <- as.matrix(rbind(aux1[[nu]],aux2,aux3))
            
            DD_aux[[nu]] <- DDcalc2(DDi_aux[[nu]],aux_mat[[nu]],numClust,Km,norm)
            
            Kmata_aux[[nu]] <- pcc_aux1[[nu]]$sil * (1 - lambda) + lambda * DD_aux[[nu]]$DD 
            
            acc[nu] <- Kmata_aux[[nu]][trimm,]
            
            nu_max <- which.max(acc)
          }
          
          NN$Kmata <- Kmata_aux[[nu_max]]
          NN$pcc <- pcc_aux1[[nu_max]] 
          NN$NN <- aux_mat[[nu_max]]
          NN$DDi <- DDi_aux[[nu_max]]
          NN$DD <- DD_aux[[nu_max]]
          NN$Nuvec <- aux1[[nu_max]]
        }  
      }          
            
      move <- c(move,NN$ct) 
      prs <- c(prs,NN$prs) 
      Ls <- c(Ls,sum(NN$Kmata)) 
      
      NN$Cost <- NN1_aux[[kl]]$Cost 
      
      if(verbose){
       cat("Cost of the partition when finishing the kl iteration: \n") 
       print(NN$Cost) 
       cat("\n")
      }
      data=data1 
      
      Nvec <- rep(1,dim(data)[1]) 
      for (kk in (1:numClust)){ 
        Xmatr <- data[NN$Nuvec==kk,] 
        Nvecr <- Nvec[NN$Nuvec==kk] 
        yr <- initW2(Xmatr,Nvecr) # Initial y for Weiszfeld algorithm 
        W <- Weisziteradj(Xmatr,Nvecr,yr,B) 
        Y[kk,] <- W$y 
      }
      rm(W,Nvecr,Xmatr,yr)
      
      if(verbose){
       cat("Deepest women: (modified Weiszfeld algorithm)\n")
       print(rownames(Y))
       cat("\n")  
      }
    }
    
     ma <- max(unlist(improv_kl),na.rm=TRUE) 
     indiv_trimmed1 <- indiv_trimmed[[ma+1]] 
     NN1_aux1 <- NN1_aux[[ma+1]] 
     
     data <- data[-NN1_aux1$trimmed,]
     Nvec <- rep(1,dim(data)[1])
    
     ntt <- matrix(0,numClust,1) 
         
     if(NN1_aux1$ct==0 & T0==0){
      qa <- qa + 1 
      if(qa > 5){ 
        kl <- niter + 1 
      } 
     }
    
     if(NN1_aux1$ct==0 & T0>0 & kl>=5){ 
      T0 <- 0 
     } 
    
     if(NN1_aux1$stopnow==1 & T0<.001){
      kl <- niter + 1 
     }  
    
    for(kk in (1:numClust)){ 
      Xmatr <- data[NN1_aux1$Nuvec==kk,] 
      Nvecr <- Nvec[NN1_aux1$Nuvec==kk] 
      yr <- initW2(Xmatr,Nvecr) 
      W <- Weisziteradj(Xmatr,Nvecr,yr,B) 
      Y[kk,] <- W$y 
    } 

    NNa <- NN1_aux1
    DD <- NNa$DD$DD
    Cost <- NNa$Cost
    NN <- NNa$NN 
  }
  
  Y <- as.integer(rownames(Y))

  return(list(NN=NN,cases=Y,DD=DD,Cost=Cost,discarded=indiv_trimmed1,klBest=ma))
} 
