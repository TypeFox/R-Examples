## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m2 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db){
# M2 derivative of the vectorized variance w$r.t. bpops
# 
# the output is a matrix with dimensions (ind_samps^2 X nbpop)
# create a (n^2 x nbpop) matrix
ns=size(xt_ind,1)^2
dv_dbeta=zeros(ns,sum(poped.db$parameters$notfixed_bpop))

k=1
if((poped.db$settings$m2_switch[1]==0 || poped.db$settings$m2_switch[1] == 1)){
#     if (poped.db$settings$m2_switch[1] == 0 && (poped.db$settings$gradff_switch[1]==0 || poped.db$settings$gradfg_switch[1] == 0 || poped.db$settings$hle_switch==0))
#         stop(sprintf('Complex differentiation to derive M2 cant be used with complex differentiaton of gradff, gradfg or hle\n'))
#     }
    
    for(i in 1:poped.db$parameters$nbpop){
        if((poped.db$parameters$notfixed_bpop[i]==1)){
            bpop_plus=bpop
            
            #if (poped.db$settings$m2_switch[1] == 1) #Central approximation or complex always uses central
                bpop_plus[i]=bpop_plus[i]+poped.db$settings$hm2
                bpop_minus=bpop
                bpop_minus[i]=bpop_minus[i]-poped.db$settings$hm2
                
                if((poped.db$settings$bCalculateEBE)){
                    #zeros(size(b_ind)[1],size(b_ind)[2])
                    start_bind = t(b_ind)
                    b_ind_plus = ind_estimates(poped.db$mean_data,bpop_plus,d,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
                    b_ind_minus = ind_estimates(poped.db$mean_data,bpop_minus,d,sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
                } else {
                    b_ind_plus = b_ind
                    b_ind_minus = b_ind
                }
                
                 returnArgs <-  v(model_switch,xt_ind,x,a,bpop_plus,b_ind_plus,bocc_ind,d,sigma,docc,poped.db) 
v_plus <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
                 returnArgs <-  v(model_switch,xt_ind,x,a,bpop_minus,b_ind_minus,bocc_ind,d,sigma,docc,poped.db) 
v_minus <- returnArgs[[1]]
poped.db <- returnArgs[[2]]
                dv=v_plus-v_minus
                if((!isempty(dv))){
                    #Central
                    ir=dv/(2*poped.db$settings$hm2)
                    ir=reshape_matlab(ir,ns,1)
                    dv_dbeta[,k]=ir
                }
            #else #Complex differentiation
            #    stop(sprintf('Complex derivative cant be used to derive M2\n'))
#                 bpop_plus[i]=bpop_plus[i]+j*poped.db$settings$hm2
#                 [v_plus,poped.db] = v(model_switch,xt_ind,x,a,bpop_plus,b_ind,bocc_ind,d,sigma,docc,poped.db)
#                 dv=v_plus
#                 if (!isempty(dv))
#                     #Complex
#                     ir=Im(dv)/(poped.db$settings$hm2)
#                     ir=reshape_matlab(ir,ns,1)
#                     dv_dbeta[,k]=ir
#                 }
           # }
            k=k+1
        }
    }
} else {
    if((poped.db$settings$m2_switch[1]==20) ){#Analytic derivative
      stop("Analytic derivatives not implemented in PopED for R")
      #         dv_dbeta=analytic_dvary_dbpop[model_switch,xt_ind,x,a,bpop,d]
#         for(i in poped.db$parameters$nbpop:-1:1){
#             if((poped.db$parameters$notfixed_bpop[i]==0)){
#                 dv_dbeta[,i]=matrix(0,0,0)
#             }
#         }
    } else {
        if((poped.db$settings$m2_switch[1]==30) ){#Automatic differentiation (INTLab), only works with "normal" variance
          stop("Automatic differentiation not implemented in PopED for R")
#              returnArgs <-  m2_ad_wrapper(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,sigma,docc,poped.db) 
# dv_dbeta <- returnArgs[[1]]
# poped.db <- returnArgs[[2]]
        } else {
         stop(sprintf('Unknown derivative option for m2'))
        }
    }
}
return(list( dv_dbeta= dv_dbeta,poped.db=poped.db)) 
}

  

