print.jzs_med <-
function(x,...){
  print(x[!(names(x)%in%c("alpha_samples","beta_samples","tau_prime_samples","rho",
                          "prob_alpha","prob_beta","prob_tau_prime_samples",
                          "ab_samples","jagssamples","jagssamplesA","jagssamplesTB"))])
}
