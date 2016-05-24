## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

gradff <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db){
  #----------Model linearization with respect to random var.
  #
  # size of return is (samples per individual x number of g's)
  #
  # derivative of model w$r.t. g eval at b=b_ind
  #
  #
  dff_dg0=zeros(size(xt_ind,1),poped.db$parameters$ng)
  fg0=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind)
  
  epsi0 = zeros(1,length(poped.db$parameters$notfixed_sigma))
  
  #Central approximation
  if((poped.db$settings$gradff_switch[1] == 1)){
    for(i in 1:poped.db$parameters$ng){
      g_plus=fg0
      g_minus=fg0
      g_plus[i]=g_plus[i]+poped.db$settings$hlf
      g_minus[i]=g_minus[i]-poped.db$settings$hlf
      returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,poped.db) 
      ferror_plus <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_minus,epsi0,poped.db) 
      ferror_minus <- returnArgs[[1]]
      poped.db <- returnArgs[[2]]
      dff_dg0[,i]=(ferror_plus-ferror_minus)/(2.0*poped.db$settings$hlf)
    }
  } else {
    #Complex approximation
    if((poped.db$settings$gradff_switch[1]==0)){
      for(i in 1:poped.db$parameters$ng){
        g_plus=fg0
        g_plus[i] = complex(real=g_plus[i],imaginary=poped.db$settings$hlf)
        returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,poped.db) 
        ferror_pointer_plus <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        dff_dg0[,i]=Im(ferror_pointer_plus)/poped.db$settings$hlf
      }
    } else {
      if((poped.db$settings$gradff_switch[1] == 20) ){#Analytic derivative
        for(k in 1:size(xt_ind,1)){
          dff_dg0[k,] = eval(sprintf('analytic_dff_dg%d(model_switch,xt_ind[k],fg0)',model_switch[k]))
        }
      } else {
        if((poped.db$settings$gradff_switch[1] == 30) ){#Calculate using automatic differentiation (INTLab)
          stop("Automatic differentiation not currently implemented in PopED for R")
          #                 fg_init = gradientinit(fg0)
          #                  returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg_init,epsi0,poped.db) 
          # val <- returnArgs[[1]]
          # poped.db <- returnArgs[[2]]
          #                 dff_dg0 = val$dx
        } else {
          stop(sprintf('Unknown derivative option for gradff'))
        }
      }
    }
  }
  return(list( dff_dg0= dff_dg0,poped.db=poped.db)) 
}
