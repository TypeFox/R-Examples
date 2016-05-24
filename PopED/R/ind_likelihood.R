#Calculates the individual ln likelihood
# function li=posthoc_dist(data,xt,x,a,bpop,d,bhat,...
#                          EPS,model_switch,...
#                          hlf,hlg,hle,...
#                          INTER)
# % from davidian and giltinan p 173 
# dmat= diag(d);
# 
# fg_bhat=fg(x,a,bpop,bhat);
# ff_bhat = ff(model_switch,xt,fg_bhat);
# 
# res = data - ff_bhat;
# 
# if INTER==true
# h=LinMatrixH_foce(model_switch,xt,x,a,bpop,EPS,hle,bhat);
# else
#   h=LinMatrixH(model_switch,xt,x,a,bpop,EPS,hle);
# end
# RR = diag(diag(h*diag(EPS)*h'));
# 
#   %li = log(det(dmat))+ (bhat'/dmat)*bhat+ ...
#           li = (bhat'/dmat)*bhat+ ...
#        log(det(RR))+ (res'/RR)*res;        

ind_likelihood <- function(data,bpop,d,sigma,bInter,bUDDLike,model_switch,xt_ind,x,a,bocc_ind,poped.db,lC,det_res_var,b_ind){
  # % from davidian and giltinan p 173 
  if((!bUDDLike)){
    #browser()
    fg_bhat=feval(poped.db$model$fg_pointer,x,a,bpop,b_ind,bocc_ind)
    ipred = feval(poped.db$model$ff_pointer,model_switch,xt_ind,fg_bhat,poped.db)[["y"]]
    res = data-ipred#Individual residuals
    
    if((bInter)){
      #For Cases WITH interaction, linearize around eta = eta^
      #eps = zeros(size(tdata),1),size(sigma,2))
      #eps = zeros(size(t(data),1),size(sigma,2))
      h = LinMatrixH(t(model_switch),
                     t(xt_ind),
                     t(x),
                     t(a),
                     bpop,
                     b_ind,
                     bocc_ind,
                     poped.db)["y"] #The covariance for this individual
      res_var = diag_matlab(diag_matlab(h%*%sigma%*%t(h)))
      lC=solve(chol(res_var),diag_matlab(length(res_var)))
      det_res_var = det(res_var)
    } else {
      #Cases WITHOUT interaction, linearize around eta = 0
      #  h = LinMatrixH(tdata,cdata,theta,zeros(size(eta)),eps) #The covariance for this individual
    }
    
    R=(t(res)%*%lC)
    li = -1/2*log(det_res_var)-1/2*R%*%t(R) # + const
    
  } else {
    #%UDD likelihood
    #li=sum(model(tdata,cdata,theta,eta))
    stop("User defined likelihood not implemented for PopED in R")  
  }
  return( li ) 
}
