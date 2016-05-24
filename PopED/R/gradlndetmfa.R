gradlndetmfa <- function(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db){

n = get_fim_size(poped.db)
m=size(ni,1)
gdmf=matrix(1,m,size(a,2))

iParallelN = (poped.db$settings$parallel$bParallelSG==1) + 1 #1 if no parallel, 2 if parallel

if((iParallelN == 2)){
    designsin = cell(1,0)
    it=1
}
for(p in 1:iParallelN){
    if((p==2)){
        #Execute parallel designs
      stop("Parallel execution not yet implemented in PopED for R")#
      designout = designsin
      #designout = execute_parallel(designsin,poped.db)
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
        if(all(size(poped.db$settings$prior_fim) == size(mft))){
            mft = mft + poped.db$settings$prior_fim
        }
        imft=inv(mft)
        if((isinf(imft[1,1]))){
            imft = zeros(size(mft))
        }
    }
    
    

for(i in 1:m){
    if((groupsize[i]==0)){
        gdmf[i,1:ni[i]]=zeros(1,ni(i))
    } else {
       if((!isempty(x))){
        x_i = t(x[i,,drop=F])
       } else {
        x_i =  zeros(0,1)
       }
       if((iParallelN ==1)){
            returnArgs <- mf_all(t(model_switch[i,1:ni[i,drop=F],drop=F]),t(xt[i,1:ni[i,drop=F],drop=F]),x_i,t(a[i,,drop=F]),bpop,d,sigma,docc,poped.db) 
mf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
           mfb=groupsize[i]*mf_tmp
       } else {
           if((p==1)){
               designsin = update_designinlist(designsin,1,ni,xt,x,a,-1,i)
           } else {
               mf_tmp = designout[[it]]$FIM
               it = it+1
               mfb=groupsize[i]*mf_tmp
           }
       }
        for(ct1 in 1:size(poped.db$design$a,2)){
            if((aa[i,ct1]!=0)){
                a_plus=a
                a_plus[i,ct1]=a_plus[i,ct1]+poped.db$settings$hgd
                if((iParallelN ==1)){
                     returnArgs <- mf_all(t(model_switch[i,1:ni[i,drop=F],drop=F]),t(xt[i,1:ni[i,drop=F],drop=F]),x_i,t(a_plus[i,,drop=F]),bpop,d,sigma,docc,poped.db) 
mf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
                 } else {
                    if((p==1)){
                        designsin = update_designinlist(designsin,1,ni,xt,x,a_plus,-1,i)
                    } else {
                        mf_tmp = designout[[it]]$FIM
                        it = it+1
                    }
                }
                if((iParallelN ==1 || p==2)){
                    mf_plus = groupsize[i]*mf_tmp
                    ir=(mf_plus-mfb)/poped.db$settings$hgd
                    s=0
                    for(ct2 in 1:n){
                        s=s+imft[ct2,,drop=F]%*%ir[,ct2,drop=F]
                    }
                    if((s!=0)){
                        gdmf[i,ct1]=s
                    } else {  #The model doesn't depend on a, e$g. PD is only dependent on time and not dose, fix the a-gradient to a small value JN
                        gdmf[i,ct1]=1e-12
                    }
                }
            }
        }
    }
}
}
ret=gdmf
return(list( ret= ret,poped.db=poped.db)) 
}





