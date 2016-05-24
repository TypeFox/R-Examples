drawdensity <-
function(fit, d = "best", xl, xu, ql = NA, qu = NA, lp = F, ex = NA, sf = 3, ind = T, lpw = 1){

  if(is.na(ql)==F & (ql <0 | ql>1 )){stop("Lower feedback quantile must be between 0 and 1")}
  if(is.na(qu)==F & (qu <0 | qu>1 )){stop("Upper feedback quantile must be between 0 and 1")}
  
  if(nrow(fit$vals)>1 & is.na(ex)==T & lp==F){
    plotgroup(fit, xl, xu, d, bw = T)
  }
  
  if(nrow(fit$vals)>1 & lp==T){
    plotlinearpool(fit, xl, xu, ql, qu , d, ind, lpw)
  }
  
  if(nrow(fit$vals)>1 & is.na(ex)==F){
    plotsingle(fit, d, xl, xu, ql, qu, sf, ex)
  }
  
  if(nrow(fit$vals)==1){
    plotsingle(fit, d, xl, xu, ql, qu, sf, ex=1)
  }	
}
