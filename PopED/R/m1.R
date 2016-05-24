## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

m1 <- function(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,d,poped.db){
  #
  # function computes the derivative of the
  # linerarized model function w$r.t. bpop
  # for an individual
  #
  # the output is a matrix with dimensions (ind_samps X nbpop)
  
  
  df_dbeta = zeros(size(xt_ind,1),sum(poped.db$parameters$notfixed_bpop))
  
  epsi0 = zeros(1,length(poped.db$parameters$notfixed_sigma))
  
  h = poped.db$settings$hm1
  
  # create linearized model
  if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==3) ){#FO, FOI
    b_ind=zeros(poped.db$parameters$NumRanEff,1)
  }
  
  if((poped.db$settings$m1_switch[1] == 1)){
    #Central approximation
    k=1
    for(i in 1:poped.db$parameters$nbpop){
      if((poped.db$parameters$notfixed_bpop[i]==1)){
        bpop_plus=bpop
        bpop_minus=bpop
        bpop_plus[i]=bpop_plus[i]+h
        bpop_minus[i]=bpop_minus[i]-h
        if((poped.db$settings$bCalculateEBE)){
          start_bind = t(b_ind)
          b_ind_plus = ind_estimates(poped.db$mean_data,bpop_plus,d,poped.db$parameters$sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
          b_ind_minus = ind_estimates(poped.db$mean_data,bpop_minus,d,poped.db$parameters$sigma,start_bind,(poped.db$settings$iApproximationMethod==2),FALSE,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db)
          
        } else {
          b_ind_plus = b_ind
          b_ind_minus = b_ind
        }
        
        g_plus=feval(poped.db$model$fg_pointer,x,a,bpop_plus,b_ind_plus,bocc_ind)
        g_minus=feval(poped.db$model$fg_pointer,x,a,bpop_minus,b_ind_minus,bocc_ind)
        
        
        if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==3 || (isempty(b_ind) && isempty(bocc_ind))) ){#FO, FOI
          returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,poped.db) 
          ferror_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_minus,epsi0,poped.db) 
          ferror_minus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          if((poped.db$settings$bUseSecondOrder)){
            hess_eta_plus = zeros(length(xt_ind),1)
            hess_eta_minus = zeros(length(xt_ind),1)
            for(o in 1:length(xt_ind)){
              hessian_eta_plus = hessian_eta_complex(model_switch[o],xt_ind[o],x,a,bpop_plus,b_ind,bocc_ind,poped.db)
              hessian_eta_minus = hessian_eta_complex(model_switch[o],xt_ind[o],x,a,bpop_minus,b_ind,bocc_ind,poped.db)
              hess_eta_plus[o] = 1/2*trace_matrix(hessian_eta_plus*d)
              hess_eta_minus[o] = 1/2*trace_matrix(hessian_eta_minus*d)
            }
            ferror_plus = ferror_plus+hess_eta_plus
            ferror_minus = ferror_minus+hess_eta_minus
          }
          df_dbeta[,k]=(ferror_plus-ferror_minus)/(2.0*h)
        } else { #FOCE, FOCEI
          returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,poped.db) 
          ferror_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop_plus,b_ind_plus,bocc_ind,poped.db) 
          l_plus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_minus,epsi0,poped.db) 
          ferror_minus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop_minus,b_ind_minus,bocc_ind,poped.db) 
          l_minus <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          
          occ_add_plus = zeros(size(xt_ind,1), 1)
          occ_add_minus = zeros(size(xt_ind,1), 1)
          
          if((isempty(b_ind)) ){#No IIV present
            l_plus = zeros(size(xt_ind,1), 1)
            l_minus = zeros(size(xt_ind,1),1)
          } else {
            l_plus = l_plus%*%b_ind_plus
            l_minus = l_minus%*%b_ind_minus
          }
          if(poped.db$parameters$NumOcc!=0){
            for(m in 1:poped.db$parameters$NumOcc){
              returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop_plus,b_ind,bocc_ind,m,poped.db) 
              l_plus_occ <- returnArgs[[1]]
              poped.db <- returnArgs[[2]]
              returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop_minus,b_ind,bocc_ind,m,poped.db) 
              l_minus_occ <- returnArgs[[1]]
              poped.db <- returnArgs[[2]]
              occ_add_plus=occ_add_plus+l_plus_occ*(bocc_ind[,m])
              occ_add_minus=occ_add_minus+l_minus_occ*(bocc_ind[,m])
            }
          }
          df_dbeta[,k]=((ferror_plus-(l_plus+occ_add_plus))-(ferror_minus-(l_minus+occ_add_minus)))/(2*h)
        }
        k=k+1
      }
    }
  } else {
    #Complex derivative
    if((poped.db$settings$m1_switch[1] == 0)){
      k=1
      for(i in 1:poped.db$parameters$nbpop){
        if((poped.db$parameters$notfixed_bpop[i]==1)){
          bpop_plus=bpop
          bpop_plus[i] = complex(real=bpop_plus[i],imaginary=h)
          g_plus=feval(poped.db$model$fg_pointer,x,a,bpop_plus,b_ind,bocc_ind)
          
          if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==3) ){#FO, FOI
            returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,poped.db) 
            ferror_tmp <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            df_dbeta[,k] = Im(ferror_tmp)/h
          } else { #FOCE, FOCEI
            returnArgs <- feval(poped.db$model$ferror_pointer,model_switch,xt_ind,g_plus,epsi0,poped.db) 
            ferror_tmp <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            #dLinMatrixL/dbpop, dLinMatrixL_occ must be central difference to assure
            #that complex step can be used within Linmatrix
            bpop_plus_c = bpop
            bpop_minus_c = bpop
            bpop_plus_c[i]=bpop_plus_c[i]+h
            bpop_minus_c[i]=bpop_minus_c[i]-h
            returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop_plus_c,b_ind,bocc_ind,poped.db) 
            l_plus <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            returnArgs <- LinMatrixL(model_switch,xt_ind,x,a,bpop_minus_c,b_ind,bocc_ind,poped.db) 
            l_minus <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            dL_dbpop = ((l_plus-l_minus))/(2*h)
            occ_add_plus = zeros(size(xt_ind,1), 1)
            occ_add_minus = zeros(size(xt_ind,1), 1)
            if((isempty(b_ind)) ){#No IIV present
              dL_dbpop = zeros(size(xt_ind,1), 1)
            } else {
              dL_dbpop = dL_dbpop*b_ind
            }
            for(m in 1:poped.db$parameters$NumOcc){
              returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop_plus_c,b_ind,bocc_ind,m,poped.db) 
              l_plus_occ <- returnArgs[[1]]
              poped.db <- returnArgs[[2]]
              returnArgs <- LinMatrixL_occ(model_switch,xt_ind,x,a,bpop_minus_c,b_ind,bocc_ind,m,poped.db) 
              l_minus_occ <- returnArgs[[1]]
              poped.db <- returnArgs[[2]]
              occ_add_plus=occ_add_plus+l_plus_occ*(bocc_ind[,m])
              occ_add_minus=occ_add_minus+l_minus_occ*(bocc_ind[,m])
            }
            df_dbeta[,k] = Im(ferror_tmp)/h-(dL_dbpop+(occ_add_plus-occ_add_minus)/(2*h))
          }
          k=k+1
        }
      }
    } else {
      if((poped.db$settings$m1_switch[1] == 20) ){#Analytic derivative
        df_dbeta_tmp = zeros(size(xt_ind,1),length(poped.db$parameters$notfixed_bpop))
        for(k in 1:size(xt_ind,1)){
          df_dbeta_tmp[k,] = eval(sprintf('analytic_dff_dbpop%d(model_switch,xt_ind[k],x,a,bpop,b_ind)',model_switch[k]))
        }
        m=1
        for(i in 1:poped.db$parameters$nbpop){
          if((poped.db$parameters$notfixed_bpop[i]==1)){
            df_dbeta[,m] = df_dbeta_tmp(,i)
            m=m+1
          }
        }
      } else {
        if((poped.db$settings$m1_switch[1] == 30) ){#Automatic differentiation using INTLab
          if((poped.db$settings$Engine$Type==2) ){#FreeMat
            stop(sprintf('Automatic differentiation is not available in PopED with FreeMat'))
          }
          if((poped.db$settings$iApproximationMethod==0 || poped.db$settings$iApproximationMethod==3 || (isempty(b_ind) && isempty(bocc_ind))) ){#FO, FOI
            stop("Automatic differentiation not currently implemented in PopED for R")
            #                   bpop_init = gradientinit(bpop)
            #                     fg_init=feval(poped.db$model$fg_pointer,x,a,bpop_init,b_ind,bocc_ind)
            #                      returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg_init,epsi0,poped.db) 
            # val <- returnArgs[[1]]
            # poped.db <- returnArgs[[2]]
            #                     df_dbeta = val$dx
            #                     for(i in poped.db$parameters$nbpop:-1:1){
            #                         if((poped.db$parameters$notfixed_bpop[i]==0)){
            #                             df_dbeta[,i]=matrix(0,0,0)
            #                         }
            #                     }
          } else { #FOCE, FOCEI
            stop("Automatic differentiation not currently implemented in PopED for R")
            #bpop_init = gradientinit(bpop)
            #             fg_init=feval(poped.db$model$fg_pointer,x,a,bpop_init,b_ind,bocc_ind)
            #             returnArgs <-  feval(poped.db$model$ferror_pointer,model_switch,xt_ind,fg_init,epsi0,poped.db) 
            #             val <- returnArgs[[1]]
            #             poped.db <- returnArgs[[2]]
            #             returnArgs <-  dLinMatrixL_dbpop[model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db] 
            #             cellDeriv <- returnArgs[[1]]
            #             L <- returnArgs[[2]]
            #             poped.db <- returnArgs[[3]]
            #             returnArgs <-  dLinMatrixL_occ_dbpop[model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db] 
            #             cellDerivOcc <- returnArgs[[1]]
            #             L_occ <- returnArgs[[2]]
            #             poped.db <- returnArgs[[3]]
            #             o = 1
            #             for(k in 1:poped.db$parameters$nbpop){
            #               if((poped.db$parameters$notfixed_bpop[k]==1)){
            #                 if((isempty(cellDeriv)) ){#Add linmatrix
            #                   l_tmp = zeros(size(xt_ind,1),1)
            #                 } else {
            #                   l_tmp = cellDeriv[[k]]*b_ind
            #                 }
            #                 occ_add = zeros(size(xt_ind,1),1)
            #                 for(m in 1:poped.db$parameters$NumOcc ){#Add occcasion
            #                   occ_add=occ_add+cellDerivOcc[[m,k]]*(bocc_ind(,m))
            #                 }
            #                 df_dbeta[,o] = val$dx(,k) - (l_tmp+occ_add)
            #                 o=o+1
            #               }
            #             }
          }
          
        } else {
          stop(sprintf('Unknown derivative option for m1'))
        }
      }
    }
  }
  return(list( df_dbeta= df_dbeta,poped.db=poped.db)) 
}

