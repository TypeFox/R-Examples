#Get the emperical bayes estimates for one individual
#Written for PopED by JN
ind_estimates <- function(data,bpop,d,sigma,start_bind,bInter,bUDDLike,model_switch,xt_ind,x,a,b_ind,bocc_ind,poped.db){
  
  b_i = t(start_bind)
  c1 = length(t(start_bind))/2*log(2*pi)
  c2 = 1/2*log(det(d))
  c3 = inv(d)
  
  if((bInter==FALSE && bUDDLike==FALSE) ){#Calculate only one variance for all samples
    #eps = zeros(size(tdata,1),size(sigma,2))
    h = LinMatrixH(model_switch,xt_ind,x,a,bpop,b_ind,bocc_ind,poped.db)
    res_var = diag_matlab(diag_matlab(h%*%sigma%*%t(h)))
    lC=solve(chol(res_var),diag_matlab(length(res_var))) #Calculate cholesky factorization
    det_res_var = det(res_var)
  } else {
    lC = 0
    det_res_var = 0
  }
  
  #This can be changed to any optimization method in MATLAB or user written
  #Could scale b_i by b_i/sqrt(OM_i)*0.01 because fminsearch uses
  #max(5%*b_i,1E-05)
  #b_z = fminsearch(@(b_z) min_function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped.db,c1,c2,c3,lC,det_res_var,b_z),b_i,optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',1000))
  opt_fun <- function(b_z){
    min_function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped.db,c1,c2,c3,lC,det_res_var,b_z)
  }
  b_z = optim(b_i,opt_fun)
  eta = b_z
  
  return( eta ) 
}

#This is the function that should be minimized, w$r.t eta
min_function <- function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped.db,c1,c2,c3,lC,det_res_var,b_ind,return_deriv=F){
  bAnalytic = FALSE
  li = ind_likelihood(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped.db,lC,det_res_var,b_ind) #Individual log likelihood
  ret =c1+c2+1/2*t(b_ind)*c3*b_ind-li
  if((return_deriv)){
    if((!bAnalytic)){
      dret=zeros(length(d),1)
      heta=1e-8
      for(i in 1:length(d)){
        bind_plus=b_ind
        bind_plus[i]=b_ind(i)+heta
        li_plus = ind_likelihood(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped.db,lC,det_res_var,bind_plus) #Individual log likelihood
        dret[i]=(li_plus-li)/heta
      }
      #}# else {
      #  dret = sum(dmodel_deta(tdata,cdata,theta,eta))
    }
    #dret=c3*eta-dret
  }
  return(list( ret= ret, dret= dret)) 
}
