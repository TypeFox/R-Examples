R_mstep2 <- function(invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,
tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_){
.Call( "R_mstep_cpp",invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,
tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_)
}

