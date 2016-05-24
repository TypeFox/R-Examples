
ret_hess <- function(final_par,A,S,F,A_fixed,A_est,S_fixed,S_est){
  mult = RAMmult(par=final_par,A,S,F,A_fixed,A_est,S_fixed,S_est)
  retH = hessian(par=final_par,ImpCov=mult$ImpCov,A,A_fixed,A_est,
                          S,S_fixed,S_est,F)
  retH
}

