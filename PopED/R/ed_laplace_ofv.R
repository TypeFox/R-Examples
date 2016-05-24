#' Evaluate the expectation of determinant the Fisher Information Matrix (FIM)
#' using the Laplace approximation.
#' 
#' Compute the expectation of the \code{det(FIM)} using the Laplace
#' approximation to the expectation. Computations are made based on the model,
#' parameters, distributions of parameter uncertainty, design and methods
#' defined in the PopED database or as arguments to the funciton.
#' 
#' This computation follows the method outlined in Dodds et al, 
#' "Robust Population Pharmacokinetic Experiment Design" JPP, 2005, equation 16.
#' 
#' Typically this function will not be run by the user.  Instead use \code{\link{evaluate.e.ofv.fim}}.
#' 
#' @param x The design parameters to compute the gradient on.
#' @inheritParams evaluate.fim
#' @inheritParams create.poped.database
#' @inheritParams Doptim
#' @param xtopto the sampling times
#' @param xopto the discrete design variables
#' @param optxt If sampling times are optimized
#' @param opta If continuous design variables are optimized
#' @param aopto the continuous design variables
#' @param method If 0 then use an optimization routine translated from poped code written in MATLAB to
#'        optimize the parameters in the Laplace approximation.  If 1 then use \code{\link{optim}} to compute both
#'        k and the hessian of k (see Dodds et al, JPP, 2005 for more information). If 2 then use \code{\link{fdHess}}
#'        to compute the hessian.
#' @param return_gradient Should the gradient be returned.
#'   
#' @return The FIM and the hessian of the FIM.
#'   
#' @family FIM
#' @family E-family
#' @example tests/testthat/examples_fcn_doc/examples_ed_laplace_ofv.R
#' @export
#' @keywords internal
# @importFrom nlme fdHess
 
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker
## right now function only works for normal and log-normal priors

ed_laplace_ofv <- function(model_switch,groupsize,ni,xtopto,xopto,aopto,
                           bpopdescr,ddescr,covd,sigma,docc,poped.db,
                           method=1,
                           return_gradient=FALSE,
                           optxt=poped.db$settings$optsw[2], 
                           opta=poped.db$settings$optsw[4],
                           x=c())
{
  
  if(any(ddescr[,1,drop=F]!=0&ddescr[,1,drop=F]!=4)){
    stop(sprintf('Only lognormal prior is supported for random effects!')) 
  }
  if(any(bpopdescr[,1,drop=F]!=0 & 
           bpopdescr[,1,drop=F]!=4 & 
           #bpopdescr[,1,drop=F]!=2 & 
           bpopdescr[,1,drop=F]!=1)){
    #stop(sprintf('Only uniform, normal and lognormal priors are supported for fixed effects!')) 
    stop(sprintf('Only normal and lognormal priors are supported for fixed effects!')) 
    
  }
  Engine = list(Type=1,Version=version$version.string)
  
  x2=x
  if(!isempty(x)){
    if(optxt){
      notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
      if(poped.db$design_space$bUseGrouped_xt){
        xtopto[notfixed]=x[poped.db$design_space$G_xt[notfixed]]
        x[1:numel(unique(poped.db$design_space$G_xt[notfixed]))]=matrix(0,0,0)
      } else {
        xtopto[notfixed]=x[1:numel(xtopto[notfixed])]
        x=x[-c(1:numel(xtopto[notfixed]))]
      }
    }
    if(opta){
      notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
      if(poped.db$design_space$bUseGrouped_a){
        aopto[notfixed]=x[poped.db$design_space$G_a[notfixed]]
      } else {
        aopto[notfixed]=x
      }
    }
    x=x2    
  }
  
  # alpha parameter vector
  alpha_k=matrix(c(bpopdescr[bpopdescr[,1]!=0,2], ddescr[ddescr[,1]!=0,2]),ncol=1,byrow=T)
  
  ## do log transformation of ln bpop and d parameters for unconstrained optimization 
  if(length(alpha_k)>sum(bpopdescr[,1]!=0)){
    d_index <- (sum(bpopdescr[,1]!=0)+1):length(alpha_k)
  } else {
    d_index <- NULL
  }  
  bpop_index=1:sum(bpopdescr[,1]!=0)
  unfixed_bpop <- bpopdescr[bpopdescr[,1]!=0,,drop=F]
  exp_index=c(unfixed_bpop[,1]==4,d_index==d_index)
  alpha_k_log=alpha_k
  if(any(exp_index)) alpha_k_log[exp_index]=log(alpha_k[exp_index])
  #alpha_k_log[d_index]=log(alpha_k[d_index])
  trans  <- function(x){ 
    x[exp_index] <- exp(x[exp_index])
    return(x)
  }
  #trans  <- function(x) matrix(c(x[bpop_index],exp(x[d_index])),ncol=1,byrow=T)
  
  #calc initial k value and gradient
  returnArgs <- calc_k(alpha_k,model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                       return_gradient=T) 
  f_k <- returnArgs[[1]]
  gf_k <- returnArgs[[2]]
  
  #   ## test f_k value
  #   fim <- det(evaluate.fim(poped.db))
  #   # assuming all normal distributions and only first 3 fixed effects
  #   p <- prod(dnorm(bpopdescr[1:3,2],mean=bpopdescr[1:3,2],sd=sqrt(bpopdescr[1:3,3]))) 
  #   k_test <- -log(fim*p)
  #   k_test==f_k
  
  #   ## test gf_k
  #   alpha_k_plus <- alpha_k+rbind(0.00001,0,0)
  #   returnArgs.1 <- calc_k(alpha_k_plus,model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
  #                        return_gradient=F) 
  #   f_k_plus <- returnArgs.1[[1]]
  # 
  #   gf_k_test.1  <- (f_k_plus - f_k)/0.00001
  
  #transform gradient for ds (log(alpha))
  #   if(!isempty(d_index)){
  #     gf_k[d_index]=gf_k[d_index]*exp(alpha_k_log[d_index])
  #   }
  if(!isempty(exp_index)){
    gf_k[exp_index]=gf_k[exp_index]*exp(alpha_k_log[exp_index])
  }
  if(isnan(f_k)){
    f=0
    gf=zeros(size(x))
    if(return_gradient) return(list( f= f,gf=gf))
    return(list(f=f))
  }
  
  ###################
  ## minimization of k(alpha) 
  ###################
  
  if(method==0){ ## sebastian method
    #initialize optimization variables
    dim=length(alpha_k)
    H_k=diag_matlab(dim)
    B_k=H_k
    niter=0
    while(norm(gf_k,type="2")>0.001){ 	# while inner conv. krit. not met
      #determine search direction for line search
      p_k=-H_k%*%gf_k
      f_name  <- "calc_k"
      #f_options <- list(trans(alpha),model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine)
      f_options <- list("replace",model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine)
      returnArgs <- line_search_uc(alpha_k_log,f_k,gf_k,p_k,f_name,f_options,exp_index)
      alpha_k1_log <- returnArgs[[1]]
      f_k1 <- returnArgs[[2]]
      s_k=alpha_k1_log-alpha_k_log
      if(max(abs(t(s_k))/max(matrix(c(t(alpha_k1_log), matrix(1,1,length(t(alpha_k1_log)))),nrow=2,byrow=T)))<1e-3){ 
        # check this that it is the same as in matlab
        break
      }
      returnArgs <- calc_k(trans(alpha_k1_log),model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                           return_gradient=T) 
      f_k1 <- returnArgs[[1]]
      gf_k1 <- returnArgs[[2]]
      #transform gradient for ds (log(alpha))
      #     if(!isempty(d_index)){
      #       gf_k1[d_index]=gf_k1[d_index]*exp(alpha_k1_log(d_index))
      #     }
      if(!isempty(exp_index)){
        gf_k1[exp_index]=gf_k1[exp_index]*exp(alpha_k1_log[exp_index])
      }
      y_k=gf_k1-gf_k
      rho_k=1/(t(y_k)%*%s_k)
      rho_k <- rho_k[,,drop=T]
      if((t(y_k)%*%s_k)/(-t(gf_k)%*%s_k) > .Machine$double.eps){
        H_k=(diag_matlab(dim)-rho_k*s_k%*%t(y_k))%*%H_k%*%(diag_matlab(dim)-rho_k*y_k%*%t(s_k))+rho_k*s_k%*%t(s_k)
      }
      alpha_k_log=alpha_k1_log
      gf_k=gf_k1
      f_k=f_k1
      niter=niter+1
    }
    alpha_k=trans(alpha_k_log)
    
    #if the number of iterations is smaller than the dimension of the problem
    #we have to calculate the hessian explicitly
    if((niter<length(B_k)||poped.db$settings$iEDCalculationType==1)){
      hess=hesskalpha2(alpha_k, model_switch,groupsize,ni,xtopto,xopto,aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,1e-6,Engine)
      detHessPi=det(hess)*(2*pi)^(-length(hess))
    } else {
      temp=matrix(1,size(gf_k))
      #temp[d_index]=1/alpha_k[d_index]
      temp[exp_index]=1/alpha_k[exp_index]
      
      iH_k=inv(H_k)
      hess=iH_k*(temp%*%t(temp))
      detHessPi=det(hess)*(2*pi)^(-length(hess))
    }
  } else { # end sebastian method
    
    #     priordescr <- rbind(bpopdescr,ddescr)
    #     priordescr <- priordescr[priordescr[,1]==2,]
    #     lb <- priordescr[,2]-priordescr[,3]/2
    #     ub <- priordescr[,2]+priordescr[,3]/2
    
    
    ## minimize K(alpha_k)
    
    output <- optim(alpha_k, 
                    function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
                                       aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                                       return_gradient=F)[["k"]],
                    #gr=function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
                    #                     aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                    #                     return_gradient=T)[["grad_k"]],
                    #method="L-BFGS-B",
                    #method="BFGS",
                    #method="Brent",
                    #lower=-100000000000,
                    #upper=0,
                    #lower=0,
                    #lower=lb,
                    #upper=ub,
                    hessian=T)
    hess <- output$hessian
    f_k <- output$value
    
    ## Different hessian and gradient calculation
    if(method==2){ 
      k_vals <- nlme::fdHess(output$par,
                             function(x) calc_k(x,model_switch,groupsize,ni,xtopto,xopto,
                                                aopto,bpopdescr,ddescr,covd,sigma,docc,poped.db,Engine,
                                                return_gradient=F)[["k"]]) 
      hess <- k_vals$Hessian
    }  
  }
  #f=Re(-exp(-f_k)/sqrt(detHessPi))
  det_hess_pi <- det(hess*(2*pi)^(-1))
  if(det_hess_pi < 0) {
    warning("The laplace OFV is ", NaN, " because det(hessian)<0.")
    f <- NaN
  } else {
    f <- sqrt(det_hess_pi)^(-1)*exp(-f_k)
  }
  
  if(return_gradient){
    bpop=bpopdescr[,2,drop=F]
    bpop[bpopdescr[,1,drop=F]!=0]=alpha_k[1:sum(bpopdescr[,1,drop=F]!=0),drop=F]
    d=ddescr[,2,drop=F]
    d[ddescr[,1]!=0]=alpha_k[sum(bpopdescr[,1,drop=F]!=0)+1:end,drop=F]
    d=getfulld(d,covd)
    
    gradxt=matrix(0,0,0)
    grada=matrix(0,0,0)
    if((optxt==TRUE)){
      notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
      gradxt=-gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
      gradxt=gradxt(notfixed)
      if(poped.db$design_space$bUseGrouped_xt){
        index=unique(poped.db$design_space$G_xt)
        gradxt=gradxt(index)
      }
    }
    if((opta==TRUE)){
      notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
      grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)                
      grada=grada(notfixed)
      if(poped.db$design_space$bUseGrouped_a){
        index=unique(poped.db$design_space$G_a)
        grada=grada(index)
      }
    }
    dkdxt=matrix(c(gradxt,grada),nrow=1,byrow=T)
    
    h_alpha=1e-4
    tensor=array(0,dim=c(length(alpha_k), length(alpha_k), length(dkdxt)))
    for(i in 1:length(alpha_k)){
      for(j in 1:i){
        alpha_plus_plus=alpha_k
        alpha_plus_plus[i]=alpha_plus_plus[i]+h_alpha
        alpha_plus_plus[j]=alpha_plus_plus[j]+h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_plus_plus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$design_space$bUseGrouped_xt){
            index=unique(poped.db$design_space$G_xt)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$design_space$bUseGrouped_a){
            index=unique(poped.db$design_space$G_a)
            grada=grada(index)
          }
        }
        dkdxt_plus_plus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        alpha_minus_plus=alpha_k
        alpha_minus_plus[i]=alpha_minus_plus[i]-h_alpha
        alpha_minus_plus[j]=alpha_minus_plus[j]+h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_minus_plus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_minus_plus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$design_space$bUseGrouped_xt){
            index=unique(poped.db$design_space$G_xt)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$design_space$bUseGrouped_a){
            index=unique(poped.db$design_space$G_a)
            grada=grada(index)
          }
        }
        dkdxt_minus_plus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        alpha_plus_minus=alpha_k
        alpha_plus_minus[i]=alpha_plus_minus[i]+h_alpha
        alpha_plus_minus[j]=alpha_plus_minus[j]-h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_plus_minus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_plus_minus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$design_space$bUseGrouped_xt){
            index=unique(poped.db$design_space$G_xt)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$design_space$bUseGrouped_a){
            index=unique(poped.db$design_space$G_a)
            grada=grada(index)
          }
        }
        dkdxt_plus_minus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        alpha_minus_minus=alpha_k
        alpha_minus_minus[i]=alpha_minus_minus[i]-h_alpha
        alpha_minus_minus[j]=alpha_minus_minus[j]-h_alpha
        bpop=bpopdescr[,2,drop=F]
        bpop[bpopdescr[,1,drop=F]!=0]=alpha_minus_minus[1:sum(bpopdescr[,1,drop=F]!=0)]
        d=ddescr[,2,drop=F]
        d[ddescr[,1]!=0]=alpha_minus_minus[sum(bpopdescr[,1,drop=F]!=0)+1:end]
        d=getfulld(d,covd)
        if((optxt==TRUE)){
          notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
          gradxt=t(gradlndetmfxt(model_switch,xtopto,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db))
          gradxt=gradxt(notfixed)
          if(poped.db$design_space$bUseGrouped_xt){
            index=unique(poped.db$design_space$G_xt)
            gradxt=gradxt(index)
          }
        }
        if((opta==TRUE)){
          notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
          grada=-gradlndetmfa(model_switch,matrix(1,size(aopto)),groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc,poped.db)
          grada=grada(notfixed)
          if(poped.db$design_space$bUseGrouped_a){
            index=unique(poped.db$design_space$G_a)
            grada=grada(index)
          }
        }
        dkdxt_minus_minus=matrix(c(gradxt,grada),nrow=1,byrow=T)
        
        tensor[i,j,]=((dkdxt_plus_plus-dkdxt_plus_minus-dkdxt_minus_plus+dkdxt_minus_minus))/(4*h_alpha^2)
        tensor[j,i,]=tensor[i,j,]
      }
    }
    ddetHessdxt=zeros(length(dkdxt),1)
    for(i in 1:length(dkdxt)){
      ddetHessdxt[i]=detHessPi*trace_matrix(inv(hess)*(-tensor[,,i]))
    }
    
    gf=Re(exp(-f_k)*(2*detHessPi*dkdxt+ddetHessdxt)/(2*detHessPi^(3/2)))
    
  }
  if(return_gradient) return(list( f= f,gf=gf))
  return(list(f=f))
}

calc_k <- function(alpha, model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,
                   ddescr,covd,sigma,docc,poped.db,Engine,return_gradient=F){
  bpop=bpopdescr[,2,drop=F]
  bpop[bpopdescr[,1,drop=F]!=0]=alpha[1:sum(bpopdescr[,1,drop=F]!=0),drop=F]
  d=ddescr[,2,drop=F]
  d[ddescr[,1]==4]=alpha[(sum(bpopdescr[,1,drop=F]!=0)+1):length(alpha),drop=F]
  d=getfulld(d,covd)
  retargs=mftot(model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,sigma,docc,poped.db)
  fim <- retargs$ret
  if((!return_gradient)){
    #tryCatch(log(det(fim)), warning = function(w) browser())
    det_fim <- det(fim)
    if(det_fim<0) det_fim <- 0
    k=-log_prior_pdf(alpha, bpopdescr, ddescr)-log(det_fim)
    grad_k=matrix(0,0,0)
  } else {
    returnArgs <- log_prior_pdf(alpha, bpopdescr, ddescr,return_gradient=T) 
    logp <- returnArgs[[1]]
    grad_p <- returnArgs[[2]]
    returnArgs <- dfimdalpha(alpha,model_switch,groupsize,ni,xtoptn,xoptn,aoptn,bpopdescr,ddescr,covd,sigma,docc,poped.db,1e-6) 
    d_fim <- returnArgs[[1]]
    fim <- returnArgs[[2]]
    ifim <- inv(fim)
    dim(ifim) <- c(length(ifim),1)
    gradlogdfim=t(reshape_matlab(d_fim,length(fim),length(grad_p)))%*%ifim
    grad_k=-(gradlogdfim+grad_p)
    ## if not positive definite set grad_k=zeros(length(alpha),1)
    
    #tryCatch(log(det(fim)), warning = function(w) browser())
    det_fim <- det(fim)
    if(det_fim<0) det_fim <- 0
    k=-logp-log(det_fim)
  }
  return(list( k= k, grad_k= grad_k)) 
}
