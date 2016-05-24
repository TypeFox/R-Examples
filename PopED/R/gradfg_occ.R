gradfg_occ <- function(x,a,bpop,b_ind,bocc_ind,currentOcc,poped.db){
#
#
# size: (number of g's x NumOccVariables)
#
# deriv of g's w$r.t. bocc's eval at b_ind, bocc_ind at occ currentOcc
#

dfg_db0=zeros(poped.db$parameters$ng,size(bocc_ind,1))

#Central approximation
if((poped.db$settings$gradfg_switch[1] == 1)){
    for(k in 1:size(bocc_ind,1)){
        bocc_plus = bocc_ind
        bocc_minus = bocc_ind
        bocc_plus[k,currentOcc] = bocc_plus[k,currentOcc] + poped.db$settings$hlg
        bocc_minus[k,currentOcc] =bocc_minus[k,currentOcc] - poped.db$settings$hlg
        dfg_db0[,k]=(feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_plus)-feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_minus))/(2.0*poped.db$settings$hlg)
    }
} else {
    #Complex approximation
    if((poped.db$settings$gradfg_switch[1] == 0)){
        for(k in 1:size(bocc_ind,1)){
            bocc_plus = bocc_ind
            bocc_plus[k,currentOcc] = complex(real=bocc_plus[k,currentOcc], imaginary=poped.db$settings$hlg)
            dfg_db0[,k]=Im(feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_plus))/poped.db$settings$hlg
        }
    } else {
        #Calculate the analytic solution at b=b_ind
        if((poped.db$settings$gradfg_switch[1] == 20)){
            stop("Automatic calculation of analytic derivatives not currently implemented in PopED for R")
            #dfg_db0 = analytic_dfg_db1(x,a,bpop,b_ind)
        } else {
            if((poped.db$settings$gradfg_switch[1] == 30) ){#Automatic differentiation (INTLab)
              stop("Automatic differentiation not currently implemented in PopED for R")
#                 bocc_init = gradientinit(bocc_ind)
#                 val = poped.db$model$fg_pointer(x,a,bpop,b_ind,bocc_init)
#                 dfg_db0 = val$dx(,(currentOcc-1)*poped.db$parameters$NumDocc+1:currentOcc*poped.db$parameters$NumDocc)
            } else {
                stop(sprintf('Unknown derivative option for gradfg_occ'))
            }
        }
    }
}
return( dfg_db0) 
}

