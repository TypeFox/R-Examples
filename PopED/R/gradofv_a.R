## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradofv_a <- function(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db){

#Input: the prior FIM or (empty) and all the other things to calculate the
#grad with for a
#Return a vector that is the gradient

if((!isempty(poped.db$settings$ed_penalty_pointer)) ){#If a penalty function is used
    #[ni, xt, model_switch, x, a, bpop, n, d, maxxt, minxt, maxa,mina]=downsizing_general_design
    na = size(a,2)
    m=size(ni,1)
    gdmf=zeros(m,na)
     returnArgs <- ed_mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$covd,sigma,docc,poped.db) 
fmf <- returnArgs[[1]]
mft <- returnArgs[[2]]
poped.db <- returnArgs[[3]]
    for(i in 1:m){
        if((groupsize[i]==0)){
            gdmf[i,1:ni[i]]=zeros(1,ni(i))
        } else {
            for(ct1 in 1:na){
                if((aa[i,ct1]!=0)){
                    a_plus=a
                    a_plus[i,ct1]=a_plus[i,ct1]+poped.db$settings$hgd
                     returnArgs <-  ed_mftot(model_switch,groupsize,ni,xt,x,a_plus,bpop,d,poped.db$parameters$covd,sigma,docc,poped.db) 
fmf_tmp <- returnArgs[[1]]
mft_plus <- returnArgs[[2]]
poped.db <- returnArgs[[3]]
                    tmp=(mft_plus-mft)/poped.db$settings$hgd
                    if((tmp==0)){
                        gdmf[i,ct1]=1e-12
                    } else {
                        gdmf[i,ct1]=tmp
                    }
                }
            }
        }
    }
    ofv_grad = gdmf
    return
}

if((poped.db$settings$ofv_calc_type==1) ){#determinant
    if((!poped.db$design_space$bUseGrouped_a)){
         returnArgs <- graddetmfa(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
ofv_grad <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
    } else {
         returnArgs <- graddetmfa_ext(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
ofv_grad <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
    }
    return
}

if((poped.db$settings$ofv_calc_type==2) ){#A-Optimal Design
    if((!poped.db$design_space$bUseGrouped_a)){
         returnArgs <-  gradtrmfa(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
ofv_grad <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
        return
    } else {
        fprintf('Warning: grad mf for grouped a on A-optimal design is not implemented analytically\nNumerical difference is used instead\n')
    }
}

if((poped.db$settings$ofv_calc_type==3) ){#S-Optimal Design
    stop(sprintf('S-optimal design for a is not implemented yet!'))
}

if((poped.db$settings$ofv_calc_type==4) ){#Log determinant
    if((!poped.db$design_space$bUseGrouped_a)){
         returnArgs <- gradlndetmfa(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
ofv_grad <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
    } else {
         returnArgs <- gradlndetmfa_ext(model_switch,aa,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
ofv_grad <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
    }
    return
}

#All other types of criterions, i$e. Ds-optimal design, CV-optimal etc.
iParallelN = (poped.db$settings$parallel$bParallelSG==1) + 1 #1 if no parallel, 2 if parallel

if((iParallelN == 2)){
    designsin = cell(1,0)
    it=1
}

if((!poped.db$design_space$bUseGrouped_a) ){#If not a grouped criterion
    na = size(a,2)
    m=size(ni,1)
    gdmf=zeros(m,na)

    for(p in 1:iParallelN){
        if((p==2)){
            #Execute parallel designs
          stop("Parallel execution not yet implemented in PopED for R")
          designout = designsin
          #designout = execute_parallel(designsin,poped.db)
        }
        if((iParallelN==1)){
             returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
mf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
        } else {
            if((p==1)){
                designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
            } else {
                mf_tmp = designout[[it]]$FIM
                it = it+1
            }
        }

        if((iParallelN==1 || p==2)){
            mft = ofv_fim(mf_tmp,poped.db)
        }

        for(i in 1:m){
            if((groupsize[i]==0)){
                gdmf[i,1:ni[i]]=zeros(1,ni(i))
            } else {
                for(ct1 in 1:na){
                    if((aa[i,ct1]!=0)){
                        a_plus=a
                        a_plus[i,ct1]=a_plus[i,ct1]+poped.db$settings$hgd

                        if((iParallelN ==1)){
                             returnArgs <-  mftot(model_switch,groupsize,ni,xt,x,a_plus,bpop,d,sigma,docc,poped.db) 
fmf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
                        } else {
                            if((p==1)){
                                designsin = update_designinlist(designsin,groupsize,ni,xt,x,a_plus,-1,0)
                            } else {
                                fmf_tmp = designout[[it]]$FIM
                                it = it+1
                            }
                        }

                        if((iParallelN ==1 || p==2)){
                            mft_plus = ofv_fim(fmf_tmp,poped.db)
                            tmp=(mft_plus-mft)/poped.db$settings$hgd
                            if((tmp==0)){
                                tmp=1e-12
                            }
                            gdmf[i,ct1]=tmp
                        }
                    }
                }
            }
        }
    }
    ofv_grad = gdmf
    return
}

if((poped.db$design_space$bUseGrouped_a) ){#If grouped a
    m=size(ni,1)
    gdmf=matrix(1,m,size(a,2))

    for(p in 1:iParallelN){
        if((p==2)){
            #Execute parallel designs
          stop("Parallel execution not yet implemented in PopED for R")
          #designout = execute_parallel(designsin,poped.db)
        }
        if((iParallelN==1)){
             returnArgs <- mftot(model_switch,groupsize,ni,xt,x,a,bpop,d,sigma,docc,poped.db) 
mf_tmp <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
        } else {
            if((p==1)){
                designsin = update_designinlist(designsin,groupsize,ni,xt,x,a,-1,0)
            } else {
                mf_tmp = designout[[it]]$FIM
                it = it+1
            }
        }

        if((iParallelN==1 || p==2)){
            mft = ofv_fim(mf_tmp,poped.db)
        }

        for(k in 1:max(max(max(poped.db$design_space$G_a)),0)){
            tmp = matrix(1,size(a,1),size(a,2))*k
            inters = (poped.db$design_space$G_a==tmp)
            if((sum(sum(inters))!=0) ){#If we have a time-point defined here (accord. to G)
                a_plus = a+poped.db$settings$hgd*inters
                if((iParallelN ==1)){
                     returnArgs <-  mftot(model_switch,groupsize,ni,xt,x,a_plus,bpop,d,sigma,docc,poped.db) 
mf_plus <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
                } else {
                    if((p==1)){
                        designsin = update_designinlist(designsin,groupsize,ni,xt,x,a_plus,-1,0)
                    } else {
                        mf_plus = designout[[it]]$FIM
                        it = it+1
                    }
                }

                if((iParallelN ==1 || p==2)){
                    mft_plus = ofv_fim(mf_plus,poped.db)
                    s=(mft_plus-mft)/poped.db$settings$hgd

                    if((s==0) ){#The model doesn't depend on a e$g.
                        s = 1e-12
                    }

                    for(i in 1:size(a,1)){
                        for(j in 1:size(a,2)){
                            if((inters[i,j]==1 && aa[i,j]!=0)){
                                gdmf[i,j]=s
                            }
                        }
                    }
                }
            }
        }
        ofv_grad = gdmf
    }
}
return(list( ofv_grad= ofv_grad,poped.db =poped.db )) 
}
