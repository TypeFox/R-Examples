graddetmfxt_ext <- function(model_switch,axt,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db){
  
  #------------------- Gradients for optimization module
  #  Looks at the gradient of det(FIM) with respect to time (xt).
  #  problems can arise when xt goes negative. So only do forward
  #  differencing. The gradient is calculated with the fact that some xt
  #  are grouped to have the same value according to the matrix G 
  #  (actually they are treated as the same variable). The
  #  derivative are calculated with:
  #  d det(A)/dx = det(A) * tr(A^-1 *dA/dX)
  
  
  iParallelN = (poped.db$settings$parallel$bParallelSG==1) + 1 #1 if no parallel, 2 if parallel
  
  if((iParallelN == 2)){
    designsin = cell(1,0)
  }
  
  n = get_fim_size(poped.db)
  m=size(ni,1)
  gdmf=matrix(1,m,size(xt,2))
  
  it = 1
  for(p in 1:iParallelN){
    
    if((p==2)){
      #Execute parallel designs
      #designout = execute_parallel(designsin,poped.db)
      stop("Parallel execution not yet implemented in PopED for R")
      designout = designsin
    }
    
    if((iParallelN==1)){
      returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
      mft <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
    } else {
      if((p==1)){
        designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
      } else {
        mft = designout[[it]]$FIM
        it = it+1
      }
    }
    
    if((iParallelN==1 || p==2)){
      #If we have a prior
      if(all(size(poped.db$settings$prior_fim)==size(mft))){
        mft = mft + poped.db$settings$prior_fim
      }
      imft=inv(mft)
      if((isinf(imft[1,1]))){
        imft = zeros(size(mft))
      }
    }
    
    for(k in 1:max(max(max(poped.db$design_space$G_xt)),0)){
      tmp = matrix(1,size(xt,1),size(xt,2))*k
      inters = (poped.db$design_space$G_xt==tmp)
      if((sum(sum(inters))!=0) ){#If we have a time-point defined here (accord. to G)
        xt_plus = xt+poped.db$settings$hgd*inters
        if((iParallelN==1)){
          returnArgs <-  mftot(model_switch,groupsize,ni,xt_plus,x,a,bpop,d,sigma,docc,poped.db) 
          mf_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
        } else {
          if((p==1)){
            designsin = update_designinlist(designsin,groupsize,ni,xt_plus,x,a,-1,0)
          } else {
            mf_plus = designout[[it]]$FIM
            it = it+1
          }
        }
        
        if((iParallelN==1 || p==2)){
          #If we have a prior
          if(all(size(poped.db$settings$prior_fim)==size(mft))){
            mf_plus = mf_plus + poped.db$settings$prior_fim
          }
          ir=(mf_plus-mft)/poped.db$settings$hgd
          
          s=0 #Calc the tr(A^-1 * dA/dX) for some X
          for(ct2 in 1:n){
            s=s+imft[ct2,,drop=F]%*%ir[,ct2,drop=F]
          }
          
          if((s==0) ){#The model doesn't depend on t e$g. PD with Placebo dose, fix the xt-gradient to a small value
            s = 1e-12
          }
          
          for(i in 1:size(xt,1)){
            for(j in 1:size(xt,2)){
              if((inters[i,j]==1 && axt[i,j]!=0)){
                gdmf[i,j]=s
              }
            }
          }
        }
      }
    }
  }
  gdmf=gdmf*det(mft)
  return(list( gdmf= gdmf,poped.db=poped.db)) 
}


