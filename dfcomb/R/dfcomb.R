CombIncrease_sim = function(ndose_a1, ndose_a2, p_tox, target, target_min, target_max, prior_tox_a1, prior_tox_a2, n_cohort, 
                        cohort, tite=FALSE, time_full=0, poisson_rate=0, nsim, c_e=0.85, c_d=0.45, c_stop=0.95, n_min=6, 
                        seed = 14061991){
  n_min=n_min-1                      
  c_d = 1-c_d
  dim_ptox = dim(p_tox)
  
  if(dim_ptox[1] != ndose_a1 || dim_ptox[2] != ndose_a2){
    stop("Wrong dimension of the matrix for true toxicity probabilities.")
  }
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }
  
  
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  target = as.double(target)[1]
  target_min = as.double(target_min)[1]
  target_max = as.double(target_max)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  n_cohort = as.integer(n_cohort)[1]
  cohort = as.integer(cohort)[1]
  tite = as.logical(tite)[1]
  time_full = as.double(time_full)[1]
  poisson_rate = as.double(poisson_rate)[1]
  nsim = as.integer(nsim)[1]
  c_e = as.double(c_e)[1]
  c_d = as.double(c_d)[1]
  c_stop = as.double(c_stop)[1]
  n_min = as.integer(n_min)[1]
  seed = as.integer(seed)[1]
  
  for(a1 in 1:ndose_a1){
    if(prior_tox_a1[a1] < 0 || prior_tox_a1[a1] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(a2 in 1:ndose_a2){
    if(prior_tox_a2[a2] < 0 || prior_tox_a2[a2] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 2 is not comprised between 0 and 1.")
    }
  }
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox[a1,a2] < 0 || p_tox[a1,a2] > 1){
        stop("At least one of the initial guessed toxicity probability is not comprised between 0 and 1.")
      }
    }
  }
  p_tox_na = matrix(NA, nrow=ndose_a1+1, ncol=ndose_a2+1)
  p_tox_na[1:ndose_a1, 1:ndose_a2] = p_tox
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox[a1,a2] >
         min(1,p_tox_na[a1+1,a2],p_tox_na[a1,a2+1],p_tox_na[a1+1,a2+1],na.rm=TRUE)){
        stop("The partial ordering between true toxicity probabilities is not satisfied.")
      }
    }
  }

  p_tox = as.double(p_tox)  
  inconc = as.double(numeric(1))
  n_pat_dose = as.double(numeric(ndose_a1*ndose_a2))
  rec_dose = as.double(numeric(ndose_a1*ndose_a2))
  n_tox_dose = as.double(numeric(ndose_a1*ndose_a2))

  
  # Appeler fonction C
  logistic = .C(C_logistic_sim, tite=tite, ndose_a1=ndose_a1, ndose_a2=ndose_a2, time_full=time_full, poisson_rate=poisson_rate,
    p_tox=p_tox, target=target, target_max=target_max, target_min=target_min, prior_tox_a1=prior_tox_a1, prior_tox_a2=prior_tox_a2,
    n_cohort=n_cohort, cohort=cohort, nsim=nsim, c_e=c_e, c_d=c_d, c_stop=c_stop, n_min=n_min, seed=seed,
    rec_dose=rec_dose, n_pat_dose=n_pat_dose, n_tox_dose=n_tox_dose, inconc=inconc)

  nsim = logistic$nsim

  inconc=logistic$inconc*100
  rec_dose=logistic$rec_dose*100
  n_pat_dose=logistic$n_pat_dose
  n_tox_dose=logistic$n_tox_dose
  
  # Reformat outputs 
  p_tox= matrix(p_tox,nrow=ndose_a1)
  rec_dose=matrix(rec_dose,nrow=ndose_a1)
  n_pat_dose=matrix(n_pat_dose,nrow=ndose_a1)
  n_tox_dose=matrix(n_tox_dose,nrow=ndose_a1)
  p_tox_p = t(p_tox)[ndose_a2:1,]
  rec_dose_p = t(rec_dose)[ndose_a2:1,]
  n_pat_dose_p = t(n_pat_dose)[ndose_a2:1,]
  n_tox_dose_p = t(n_tox_dose)[ndose_a2:1,]
  dimnames(p_tox_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(rec_dose_p) = list("Agent 2 " = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_pat_dose_p) = list("Agent 2"=ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_tox_dose_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  pat_tot = round(sum(n_pat_dose),1)

  res = list(call = match.call(),
             tite=tite, 
             ndose_a1=ndose_a1, 
             ndose_a2=ndose_a2, 
             time_full=time_full, 
             poisson_rate=poisson_rate,
             p_tox=p_tox,
             p_tox_p=p_tox_p, 
             target=target, 
             target_min=target_min,
             target_max=target_max, 
             prior_tox_a1=prior_tox_a1, 
             prior_tox_a2=prior_tox_a2,
             n_cohort=n_cohort, 
             cohort=cohort, 
             pat_tot=pat_tot,
             nsim=nsim, 
             c_e=c_e, 
             c_d=c_d, 
             c_stop=c_stop, 
             n_min=n_min, 
             seed=seed,
             rec_dose=rec_dose, 
             n_pat_dose=n_pat_dose, 
             n_tox_dose=n_tox_dose, 
             rec_dose_p=rec_dose_p, 
             n_pat_dose_p=n_pat_dose_p, 
             n_tox_dose_p=n_tox_dose_p,
             inconc=inconc)
             
  class(res) = "CombIncrease_sim"

  return(res)
}



print.CombIncrease_sim = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 

  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("True toxicities:", x$p_tox_p)
  print_rnd("Percentage of Selection:", x$rec_dose_p)
  print_rnd("Number of patients:" , x$n_pat_dose_p)
  print_rnd("Number of toxicities:", x$n_tox_dose_p)

  cat(paste("Percentage of inconclusive trials:\t",x$inconc,"\n",sep=""), sep="")
  nmin_c = which(x$cohort*(0:x$n_cohort) >= x$n_min)[1]
  cat("The minimum number of cohorts to stop the trial is:\t", nmin_c-1, "\n")
  cat("\n", "\n")              
  cat("Number of simulations:\t", x$nsim, "\n")
  cat("Cohort size:\t", x$cohort, "\n")
  cat("Number of cohort planned:\t", x$n_cohort, "\n")
  cat("Total patients accrued:\t", x$pat_tot, "\n")
  cat("Toxicity target:\t", x$target, "\n")
  cat("Targeted toxicity interval:\t [", x$target_min, ",", x$target_max, "]\n")
  cat("Prior toxicity probabilities for agent 1:\n")
  print(round(x$prior_tox_a1, digits = dgt))
  cat("Prior toxicity probabilities for agent 2:\n")
  print(round(x$prior_tox_a2, digits = dgt))
  cat("Escalation threshold:\t", x$c_e, "\n")
  cat("Deescalation threshold:\t", 1-x$c_d, "\n")
  cat("Stopping threshold:\t", x$c_stop, "\n")
  if (x$tite) {
    cat("Toxicity is a time-to-event \n")
    cat("Full follow-up time:\t", x$time_full, "\n")
    cat("Patient arrival is modeled as a Poisson process with rate:", x$poisson_rate, "\n") 
  }
  else{
    cat("Toxicity is not a time-to-event but binary \n")
  }
}



CombIncrease_next = function(ndose_a1, ndose_a2, target, target_min, target_max, prior_tox_a1, prior_tox_a2, in_startup = TRUE, final, 
                         pat_incl, dose_adm1, dose_adm2, tite=FALSE, toxicity, time_full=0, time_tox=0, time_follow=0,
                         c_e=0.85, c_d=0.45, c_stop=0.95, n_min){  
                         
  if(tite == TRUE) {
    toxicity = as.numeric(time_tox < time_follow)
  }
  if(pat_incl > 0) {
    cdose1 = dose_adm1[pat_incl]
    cdose2 = dose_adm2[pat_incl]
  }
  else {
    cdose1 = 0
    cdose2 = 0
  }
 
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }
  
  n_toxicity = length(toxicity)
  n_time_follow = length(time_follow)
  n_time_tox = length(time_tox)
  n_dose_adm1 = length(dose_adm1)
  n_dose_adm2 = length(dose_adm2)
  if(tite==FALSE && n_toxicity != pat_incl){
    stop("The entered vector of observed toxicities is of wrong length.")
  }
  if(tite==TRUE && n_time_follow != pat_incl){
    stop("The entered vector for patients' follow-up time is of wrong length.")
  }
  if(tite==TRUE && n_time_tox != pat_incl){
    stop("The entered vector for patients' time-to-toxicity is of wrong length.")
  }
  if(n_dose_adm1 != pat_incl){
    stop("The entered vector for patients' dose of agent 1 is of wrong length.")
  }
  if(n_dose_adm2 != pat_incl){
    stop("The entered vector for patients' dose of agent 2 is of wrong length.")
  }
  
  tite = as.logical(tite)
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  time_full = as.double(time_full)[1]
  target = as.double(target)[1] 
  target_max = as.double(target_max)[1]
  target_min = as.double(target_min)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  final = as.logical(final)
  c_e = as.double(c_e)[1]
  c_d = as.double(c_d)[1]
  c_stop = as.double(c_stop)[1]
  n_min = as.integer(n_min)[1]
  pat_incl = as.integer(pat_incl)[1]
  cdose1 = as.integer(cdose1-1)
  cdose2 = as.integer(cdose2-1)
  in_startup = as.logical(in_startup)[1]
  dose_adm1 = as.integer(dose_adm1-1)
  dose_adm2 = as.integer(dose_adm2-1)
  time_tox = as.double(time_tox)
  time_follow = as.double(time_follow)
  toxicity = as.logical(toxicity)
  
  for(i in 1:ndose_a1){
    if(prior_tox_a1[i] < 0 || prior_tox_a1[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(i in 1:ndose_a2){
    if(prior_tox_a2[i] < 0 || prior_tox_a2[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 2 is not comprised between 0 and 1.")
    }
  }
  if(target < 0 || target > 1){stop("The toxicity target is not comprised between 0 and 1.")}
  if(target_max < 0 || target_max > 1){stop("The maximum of the targeted toxicity interval is not comprised between 0 and 1.")}
  if(target_min < 0 || target_min > 1){stop("The minimum of the targeted toxicity interval is not comprised between 0 and 1.")}
  
                
  inconc = as.logical(numeric(1))
  pi = as.double(numeric(ndose_a1*ndose_a2))
  ptox_inf = as.double(numeric(ndose_a1*ndose_a2))
  ptox_inf_targ = as.double(numeric(ndose_a1*ndose_a2))
  ptox_targ = as.double(numeric(ndose_a1*ndose_a2))
  ptox_sup_targ = as.double(numeric(ndose_a1*ndose_a2))
  
  logistic = .C(C_logistic_next, tite=tite, ndose_a1=ndose_a1, ndose_a2=ndose_a2, time_full=time_full, target=target, 
                target_max=target_max, target_min=target_min, prior_tox_a1=prior_tox_a1, prior_tox_a2=prior_tox_a2,
                final=final, c_e=c_e, c_d=c_d, c_stop=c_stop, n_min=n_min, pat_incl=pat_incl, cdose1=cdose1,
                cdose2=cdose2, dose_adm1=dose_adm1, dose_adm2=dose_adm2, time_tox=time_tox,
                time_follow=time_follow, toxicity=toxicity, in_startup=in_startup, inconc=inconc, pi=pi,
                ptox_inf=ptox_inf, ptox_inf_targ=ptox_inf_targ, ptox_targ=ptox_targ, ptox_sup_targ=ptox_sup_targ,
                NAOK=TRUE)

  # Reformat outputs       
  in_startup = logistic$in_startup
  cdose1=logistic$cdose1+1
  cdose2=logistic$cdose2+1
  dose_adm1=dose_adm1+1 
  dose_adm2=dose_adm2+1 
  pi=matrix(logistic$pi, nrow=ndose_a1)
  ptox_inf=matrix(logistic$ptox_inf, nrow=ndose_a1)
  ptox_inf_targ=matrix(logistic$ptox_inf_targ, nrow=ndose_a1)
  ptox_targ=matrix(logistic$ptox_targ, nrow=ndose_a1)
  ptox_sup_targ=matrix(logistic$ptox_sup_targ, nrow=ndose_a1)
        
  n_pat_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  n_tox_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  for(i in 1:pat_incl){
    n_pat_comb[dose_adm1[i],dose_adm2[i]] = n_pat_comb[dose_adm1[i],dose_adm2[i]]+1
    n_tox_comb[dose_adm1[i],dose_adm2[i]] = n_tox_comb[dose_adm1[i],dose_adm2[i]]+toxicity[i]
  }
  n_pat_comb_p = t(n_pat_comb)[ndose_a2:1,]
  n_tox_comb_p = t(n_tox_comb)[ndose_a2:1,]
  pi_p = t(pi)[ndose_a2:1,]
  ptox_inf_p = t(ptox_inf)[ndose_a2:1,]
  ptox_inf_targ_p = t(ptox_inf_targ)[ndose_a2:1,]
  ptox_targ_p = t(ptox_targ)[ndose_a2:1,]
  ptox_sup_targ_p = t(ptox_sup_targ)[ndose_a2:1,]
  dimnames(n_pat_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(n_tox_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(pi_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_inf_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_inf_targ_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_targ_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_sup_targ_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  startup_in = ifelse(in_startup, "NO", "YES")

  res = list(call = match.call(),
             tite=tite, 
             ndose_a1=ndose_a1, 
             ndose_a2=ndose_a2, 
             time_full=time_full, 
             target=target, 
             target_max=target_max, 
             target_min=target_min, 
             prior_tox_a1=prior_tox_a1, 
             prior_tox_a2=prior_tox_a2,
             final=final, 
             c_e=c_e, 
             c_d=c_d, 
             c_stop=c_stop, 
             n_min=n_min, 
             pat_incl=pat_incl,
             cdose1=cdose1,
             cdose2=cdose2, 
             startup_in=startup_in, 
             dose_adm1=dose_adm1, 
             dose_adm2=dose_adm2, 
             time_tox=time_tox,
             time_follow=time_follow, 
             toxicity=toxicity, 
             inconc=logistic$inconc, 
             n_pat_comb=n_pat_comb,
             n_tox_comb=n_tox_comb,
             pi=pi,
             ptox_inf=ptox_inf, 
             ptox_inf_targ=ptox_inf_targ, 
             ptox_targ=ptox_targ,
             ptox_sup_targ=ptox_sup_targ,
             n_pat_comb_p=n_pat_comb_p,
             n_tox_comb_p=n_tox_comb_p,
             pi_p=pi_p,
             ptox_inf_p=ptox_inf_p, 
             ptox_inf_targ_p=ptox_inf_targ_p, 
             ptox_targ_p=ptox_targ_p,
             ptox_sup_targ_p=ptox_sup_targ_p)
             
  class(res) = "CombIncrease_next"

  return(res)
}




print.CombIncrease_next = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 
     
  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("Number of patients:" , x$n_pat_comb_p)
  print_rnd("Number of toxicities:", x$n_tox_comb_p)
  print_rnd("Toxicity prob:", x$pi_p)
  print_rnd("P(toxicity prob < target):", x$ptox_inf_p)
  print_rnd("Prob underdosing:", x$ptox_inf_targ_p)
  print_rnd("Prob targeted interval:", x$ptox_targ_p)
  print_rnd("Prob overdosing:", x$ptox_sup_targ_p)
  
  cat("Start-up phase ended:\t", x$startup_in,"\n")
  if(!x$inconc){
    if(x$final){
      cat(paste("RECOMMENDED COMBINATION at the end of the trial:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
    else{
      cat(paste("NEXT RECOMMENDED COMBINATION:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
  }
  else{
    cat(paste("THE DOSE-FINDING PROCESS SHOULD BE STOPPED WITHOUT COMBINATION RECOMMENDATION\n",sep=""), sep="")
  }
  cat("\n", "\n")
  
  cat("Number of patients included:\t", x$pat_incl, "\n")
  cat("Toxicity target:\t", x$target, "\n")
  cat("Targeted toxicity interval:\t [", x$target_min, ",", x$target_max, "]\n")
  cat("Prior toxicity probabilities for agent 1:\n")
  print(round(x$prior_tox_a1, digits = dgt))
  cat("Prior toxicity probabilities for agent 2:\n")
  print(round(x$prior_tox_a2, digits = dgt))
  cat("The minimum number of patients to stop the trial is:\t", x$n_min, "\n")
  cat("Escalation threshold:\t", x$c_e, "\n")
  cat("Deescalation threshold:\t", 1-x$c_d, "\n")
  cat("Stopping threshold:\t", x$c_stop, "\n")
  if (x$tite) {
    cat("Toxicity is a time-to-event \n")
    cat("Full follow-up time:\t", x$time_full, "\n")
  }
  else{
    cat("Toxicity is not a time-to-event but binary \n")
  }
}


#############################################

CombPlateau_sim = function(ndose_a1, ndose_a2, p_tox, p_eff, tox_max, eff_min, prior_tox_a1, prior_tox_a2, prior_eff_a1, prior_eff_a2, n, 
                        cohort_start=3, cohort=3, time_full, poisson_rate, cycle=0, nsim, c_tox=0.85, c_eff=0.10, seed = 2174892, threads=0){
  dim_ptox = dim(p_tox)
  dim_peff = dim(p_eff)
  if(dim_ptox[1] != ndose_a1 || dim_ptox[2] != ndose_a2){
    stop("Wrong dimension of the matrix for true toxicity probabilities.")
  }
  if(dim_peff[1] != ndose_a1 || dim_peff[2] != ndose_a2){
    stop("Wrong dimension of the matrix for true efficacy probabilities.")
  }
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }
  n_prior_eff_a1 = length(prior_eff_a1)
  if(n_prior_eff_a1 != ndose_a1){
    stop("The entered vector of initial guessed efficacy probabities for agent 1 is of wrong length.")
  }
  n_prior_eff_a2 = length(prior_eff_a2)
  if(n_prior_eff_a2 != ndose_a2){
    stop("The entered vector of initial guessed efficacy probabities for agent 2 is of wrong length.")
  }
  
  
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]  
  tox_max = as.double(tox_max)[1]
  eff_min = as.double(eff_min)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  prior_eff_a1 = as.double(prior_eff_a1)
  prior_eff_a2 = as.double(prior_eff_a2) 
  n = as.integer(n)[1]
  cohort_start = as.integer(cohort_start)[1]
  cohort = as.integer(cohort)[1]
  time_full = as.double(time_full)[1]
  poisson_rate = as.double(poisson_rate)[1]
  cycle = as.double(cycle)[1]
  nsim = as.integer(nsim)[1]
  c_tox = as.double(c_tox)[1]
  c_eff = as.double(c_eff)[1]
  threads = as.integer(threads)[1]
  seed = as.integer(seed)[1]
  
  for(a1 in 1:ndose_a1){
    if(prior_tox_a1[a1] < 0 || prior_tox_a1[a1] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 1 is not comprised between 0 and 1.")
    }
    if(prior_eff_a1[a1] < 0 || prior_eff_a1[a1] > 1){
      stop("At least one of the initial guessed efficacy probability for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(a2 in 1:ndose_a2){
    if(prior_tox_a2[a2] < 0 || prior_tox_a2[a2] > 1){
      stop("At least one of the initial guessed toxicity probability for agent 2 is not comprised between 0 and 1.")
    }
    if(prior_eff_a2[a2] < 0 || prior_eff_a2[a2] > 1){
      stop("At least one of the initial guessed efficacy probability for agent 2 is not comprised between 0 and 1.")
    }
  }
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox[a1,a2] < 0 || p_tox[a1,a2] > 1){
        stop("At least one of the initial guessed toxicity probability is not comprised between 0 and 1.")
      }
      if(p_eff[a1,a2] < 0 || p_eff[a1,a2] > 1){
        stop("At least one of the initial guessed efficacy probability is not comprised between 0 and 1.")
      }
    }
  }
  p_tox_na = matrix(NA, nrow=ndose_a1+1, ncol=ndose_a2+1)
  p_tox_na[1:ndose_a1, 1:ndose_a2] = p_tox
  p_eff_na = matrix(NA, nrow=ndose_a1+1, ncol=ndose_a2+1)
  p_eff_na[1:ndose_a1, 1:ndose_a2] = p_eff
  for(a1 in 1:ndose_a1){
    for(a2 in 1:ndose_a2){
      if(p_tox[a1,a2] >
         min(1,p_tox_na[a1+1,a2],p_tox_na[a1,a2+1],p_tox_na[a1+1,a2+1],na.rm=TRUE)){
        stop("Toxicity probabilities are not increasing with the dose of both agents.")
      }
      if(p_eff[a1,a2] >
         min(1,p_eff_na[a1+1,a2],p_eff_na[a1,a2+1],p_eff_na[a1+1,a2+1],na.rm=TRUE)){
        stop("Efficacy probabilities are not increasing or plateaus.")
      }
    }
  }

  p_tox = as.double(p_tox)
  p_eff = as.double(p_eff)  
  inconc = as.integer(numeric(1))
  n_pat_dose = as.integer(numeric(ndose_a1*ndose_a2))
  rec_dose = as.integer(numeric(ndose_a1*ndose_a2))
  n_tox_dose = as.integer(numeric(ndose_a1*ndose_a2))
  n_eff_dose = as.integer(numeric(ndose_a1*ndose_a2))
  duration = as.double(numeric(1))
   
  # Appeler fonction C
  plateau = .C(C_plateau_sim, ndose_a1=ndose_a1, ndose_a2=ndose_a2, p_tox=p_tox, p_eff=p_eff, tox_max=tox_max,
                 eff_min=eff_min, prior_tox_a1=prior_tox_a1,prior_tox_a2=prior_tox_a2, prior_eff_a1=prior_eff_a1,
                 prior_eff_a2=prior_eff_a2, poisson_rate=poisson_rate,n=n, time_full=time_full, cycle=cycle,
                 cohort_start=cohort_start, cohort=cohort, nsim=nsim, c_tox=c_tox, c_eff=c_eff, seed=seed, threads=threads,
                 inconc=inconc, n_pat_dose=n_pat_dose, rec_dose=rec_dose, n_tox_dose=n_tox_dose, n_eff_dose=n_eff_dose, 
                 duration=duration)

  nsim = plateau$nsim
  inconc=plateau$inconc/nsim*100
  rec_dose=plateau$rec_dose/nsim*100
  n_pat_dose=plateau$n_pat_dose/nsim
  n_tox_dose=plateau$n_tox_dose/nsim
  n_eff_dose=plateau$n_eff_dose/nsim
  duration=plateau$duration/nsim
  
  # Reformat outputs 
  p_tox= matrix(p_tox,nrow=ndose_a1)
  p_eff= matrix(p_eff,nrow=ndose_a1)
  rec_dose=matrix(rec_dose,nrow=ndose_a1) 
  n_pat_dose=matrix(n_pat_dose,nrow=ndose_a1) 
  n_tox_dose=matrix(n_tox_dose,nrow=ndose_a1) 
  n_eff_dose=matrix(n_eff_dose,nrow=ndose_a1)
  p_tox_p = t(p_tox)[ndose_a2:1,]
  p_eff_p = t(p_eff)[ndose_a2:1,]
  rec_dose_p = t(rec_dose)[ndose_a2:1,]
  n_pat_dose_p = t(n_pat_dose)[ndose_a2:1,]
  n_tox_dose_p = t(n_tox_dose)[ndose_a2:1,]
  n_eff_dose_p = t(n_eff_dose)[ndose_a2:1,]
  dimnames(p_tox_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(p_eff_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(rec_dose_p) = list("Agent 2 " = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_pat_dose_p) = list("Agent 2"=ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_tox_dose_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  dimnames(n_eff_dose_p) = list("Agent 2" = ndose_a2:1, "Agent 1" = 1:ndose_a1)
  pat_tot = round(sum(n_pat_dose),1)

  res = list(call = match.call(),
             ndose_a1=ndose_a1, 
             ndose_a2=ndose_a2, 
             time_full=time_full, 
             poisson_rate=poisson_rate,
             p_tox=p_tox,
             p_eff=p_eff,
             p_tox_p=p_tox_p, 
             p_eff_p=p_eff_p,
             tox_max=tox_max, 
             eff_min=eff_min, 
             prior_tox_a1=prior_tox_a1, 
             prior_tox_a2=prior_tox_a2,
             prior_eff_a1=prior_eff_a1, 
             prior_eff_a2=prior_eff_a2,
             n=n, 
             cohort_start=cohort_start,
             cohort=cohort, 
             pat_tot=pat_tot,
             cycle=cycle,
             nsim=nsim, 
             c_tox=c_tox, 
             c_eff=c_eff, 
             seed=seed,
             rec_dose=rec_dose, 
             n_pat_dose=n_pat_dose, 
             n_tox_dose=n_tox_dose,
             n_eff_dose=n_eff_dose, 
             rec_dose_p=rec_dose_p, 
             n_pat_dose_p=n_pat_dose_p, 
             n_tox_dose_p=n_tox_dose_p,
             n_eff_dose_p=n_eff_dose_p,
             inconc=inconc,
             duration=duration)
             
  class(res) = "CombPlateau_sim"

  return(res)                     
}



print.CombPlateau_sim = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 

  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("True toxicities:", x$p_tox_p)
  print_rnd("True efficacies:", x$p_eff_p)
  print_rnd("Percentage of Selection:", x$rec_dose_p)
  print_rnd("Number of patients:" , x$n_pat_dose_p)
  print_rnd("Number of toxicities:", x$n_tox_dose_p)
  print_rnd("Number of efficacies:", x$n_eff_dose_p)

  cat(paste("Percentage of inconclusive trials:\t",x$inconc,"\n",sep=""), sep="")
  cat("\n", "\n")              
  cat("Number of simulations:\t", x$nsim, "\n")
  cat("Cohort size for the start-up phase:\t", x$cohort_start, "\n")
  cat("Cohort size for the model-based phase:\t", x$cohort, "\n")
  cat("Total sample size:\t", x$n, "\n")
  cat("Total patients accrued:\t", x$pat_tot, "\n")
  cat("Toxicity upper bound:\t", x$tox_max, "\n")
  cat("Efficacy lower bound:\t", x$eff_min, "\n")
  cat("Prior toxicity probabilities for agent 1:\n")
  print(round(x$prior_tox_a1, digits = dgt))
  cat("Prior toxicity probabilities for agent 2:\n")
  print(round(x$prior_tox_a2, digits = dgt))
  cat("Prior efficacy probabilities for agent 1:\n")
  print(round(x$prior_eff_a1, digits = dgt))
  cat("Prior efficacy probabilities for agent 2:\n")
  print(round(x$prior_eff_a2, digits = dgt))
  cat("Toxicity threshold:\t", x$c_tox, "\n")
  cat("Efficacy threshold:\t", x$c_eff, "\n")
  cat("Toxicity is binary \n")
  cat("Efficacy is a time-to-event \n")
  cat("Full follow-up time:\t", x$time_full, "\n")
  cat("Minimum waiting time between two dose cohorts is of one toxicity cycle of:\t", x$cycle, "\n")
  cat("Patient arrival is modeled as a Poisson process with rate:", x$poisson_rate, "\n")
  cat("Trial mean duration:", x$duration, "\n") 
}


                     
CombPlateau_next = function(ndose_a1, ndose_a2, tox_max, eff_min, prior_tox_a1, prior_tox_a2, prior_eff_a1, prior_eff_a2, stage, in_startup,
cohort_start=3, cohort, pat_incl, dose_adm1, dose_adm2, toxicity, time_full, time_prog, time_follow, cycle=0, c_tox=0.85, c_eff=0.10){ 
   
  ndose_a1 = as.integer(ndose_a1)[1]
  ndose_a2 = as.integer(ndose_a2)[1]
  if(pat_incl > 0) {
    cdose1 = dose_adm1[pat_incl]
    cdose2 = dose_adm2[pat_incl]
  }
  else {
    cdose1 = 0
    cdose2 = 0
  }
  time_cur = 0
  time_incl = -time_follow
  time_full = as.double(time_full)[1]
  tox_max = as.double(tox_max)[1]
  eff_min = as.double(eff_min)[1]
  prior_tox_a1 = as.double(prior_tox_a1)
  prior_tox_a2 = as.double(prior_tox_a2)
  prior_eff_a1 = as.double(prior_eff_a1)
  prior_eff_a2 = as.double(prior_eff_a2)
  stage = as.integer(stage)[1]
  c_tox = as.double(c_tox)[1]
  c_eff = as.double(c_eff)[1]
  cycle = as.double(cycle)[1]
  pat_incl = as.integer(pat_incl)[1]
  cohort_start = as.integer(cohort_start)[1]
  cohort = as.integer(cohort)[1]

  n_toxicity = length(toxicity)
  n_time_follow = length(time_follow)
  n_time_prog = length(time_prog)
  n_dose_adm1 = length(dose_adm1)
  n_dose_adm2 = length(dose_adm2)
  if(n_toxicity != pat_incl){
    stop("The entered vector of observed toxicities is of wrong length.")
  }
  if(n_time_follow != pat_incl){
    stop("The entered vector for patients' follow-up time is of wrong length.")
  }
  if(n_time_prog != pat_incl){
    stop("The entered vector for patients' time-to-progression is of wrong length.")
  }
  if(n_dose_adm1 != pat_incl){
    stop("The entered vector for patients' dose of agent 1 is of wrong length.")
  }
  if(n_dose_adm2 != pat_incl){
    stop("The entered vector for patients' dose of agent 2 is of wrong length.")
  }  
    
  n_pat_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  n_tox_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  n_prog_comb = matrix(0, nrow=ndose_a1, ncol=ndose_a2)
  for(i in 1:pat_incl){
    n_pat_comb[dose_adm1[i],dose_adm2[i]] = n_pat_comb[dose_adm1[i],dose_adm2[i]]+1
    n_tox_comb[dose_adm1[i],dose_adm2[i]] = n_tox_comb[dose_adm1[i],dose_adm2[i]]+toxicity[i]
    n_prog_comb[dose_adm1[i],dose_adm2[i]] = n_prog_comb[dose_adm1[i],dose_adm2[i]]+(time_prog[i]<time_follow[i])
  }
 
  n_prior_tox_a1 = length(prior_tox_a1)
  if(n_prior_tox_a1 != ndose_a1){
    stop("The entered vector of initial guessed toxicity probabities for agent 1 is of wrong length.")
  }
  n_prior_tox_a2 = length(prior_tox_a2)
  if(n_prior_tox_a2 != ndose_a2){
    stop("The entered vector of initial guessed toxicity probabities for agent 2 is of wrong length.")
  }
  n_prior_eff_a1 = length(prior_eff_a1)
  if(n_prior_eff_a1 != ndose_a1){
    stop("The entered vector of initial guessed efficacy probabities for agent 1 is of wrong length.")
  }
  n_prior_eff_a2 = length(prior_eff_a2)
  if(n_prior_eff_a2 != ndose_a2){
    stop("The entered vector of initial guessed efficacy probabities for agent 2 is of wrong length.")
  }
  
  for(i in 1:ndose_a1){
    if(prior_tox_a1[i] < 0 || prior_tox_a1[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 1 is not comprised between 0 and 1.")
    }
    if(prior_eff_a1[i] < 0 || prior_eff_a1[i] > 1){
      stop("At least one of the initial guessed efficacy for agent 1 is not comprised between 0 and 1.")
    }
  }
  for(i in 1:ndose_a2){
    if(prior_tox_a2[i] < 0 || prior_tox_a2[i] > 1){
      stop("At least one of the initial guessed toxicity for agent 2 is not comprised between 0 and 1.")
    }
    if(prior_eff_a2[i] < 0 || prior_eff_a2[i] > 1){
      stop("At least one of the initial guessed efficacy for agent 2 is not comprised between 0 and 1.")
    }
  }
  if(tox_max < 0 || tox_max > 1){stop("The maximum toxicity is not comprised between 0 and 1.")}
  if(eff_min < 0 || eff_min > 1){stop("The minimum efficacy is not comprised between 0 and 1.")}

  if(stage != 0 && stage != 1 && stage != 2){stop("The stage must be 0, 1 or 2.")}

  cdose1 = as.integer(cdose1-1)
  cdose2 = as.integer(cdose2-1)
  time_cur = as.double(time_cur)[1]
  dose_adm1 = as.integer(dose_adm1-1)
  dose_adm2 = as.integer(dose_adm2-1)
  time_prog = as.double(time_prog)
  time_incl = as.double(time_incl)
  toxicity = as.logical(toxicity)
                
  pi = as.double(numeric(ndose_a1*ndose_a2))
  ptox_sup = as.double(numeric(ndose_a1*ndose_a2))
  resp = as.double(numeric(ndose_a1*ndose_a2))
  qeff_min = as.double(numeric(ndose_a1*ndose_a2))
  proba_tau = as.double(numeric(ndose_a2))
  in_startup = as.logical(in_startup)[1]

  plateau = .C(C_plateau_next, ndose_a1=ndose_a1, ndose_a2=ndose_a2, tox_max=tox_max, eff_min=eff_min, 
               prior_tox_a1=prior_tox_a1, prior_tox_a2=prior_tox_a2, prior_eff_a1=prior_eff_a1, prior_eff_a2=prior_eff_a2,
               time_full=time_full, cycle=cycle, cohort_start=cohort_start, cohort=cohort,
               c_tox=c_tox, c_eff=c_eff, cdose1=cdose1, cdose2=cdose2, in_startup=in_startup, time_cur=time_cur,
               pat_incl=pat_incl, dose_adm1=dose_adm1, dose_adm2=dose_adm2, 
               time_prog=time_prog, time_incl=time_incl, toxicity=toxicity,
               stage=stage,

               pi=pi, ptox_sup=ptox_sup, resp=resp,
               qeff_min=qeff_min, proba_tau=proba_tau,
               NAOK=TRUE)
                
  # Reformat outputs
  in_startup = plateau$in_startup
  cdose1=plateau$cdose1+1
  cdose2=plateau$cdose2+1
  dose_adm1=plateau$dose_adm1+1 
  dose_adm2=plateau$dose_adm2+1 
  pi=matrix(plateau$pi, nrow=ndose_a1)
  ptox_sup=matrix(plateau$ptox_sup, nrow=ndose_a1)
  resp=matrix(plateau$resp, nrow=ndose_a1)
  qeff_min=matrix(plateau$qeff_min, nrow=ndose_a1)
  proba_tau=plateau$proba_tau
         
  n_pat_comb_p = t(n_pat_comb)[ndose_a2:1,]
  n_tox_comb_p = t(n_tox_comb)[ndose_a2:1,]
  n_prog_comb_p = t(n_prog_comb)[ndose_a2:1,]
  pi_p = t(pi)[ndose_a2:1,]
  ptox_sup_p = t(ptox_sup)[ndose_a2:1,]
  resp_p = t(resp)[ndose_a2:1,]
  qeff_min_p = t(qeff_min)[ndose_a2:1,]
  proba_tau_p = t(proba_tau)[ndose_a2:1]
  dimnames(n_pat_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(n_tox_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(n_prog_comb_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(pi_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(ptox_sup_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(resp_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  dimnames(qeff_min_p) = list("Agent 2"=ndose_a2:1, "Agent 1"=1:ndose_a1)
  names(proba_tau_p) = c(ndose_a2:1)
  startup_in = ifelse(in_startup, "NO", "YES")

  res = list(call = match.call(),
             ndose_a1=ndose_a1, 
             ndose_a2=ndose_a2, 
             time_full=time_full, 
             cycle=cycle, 
             tox_max=tox_max, 
             eff_min=eff_min, 
             prior_tox_a1=prior_tox_a1, 
             prior_tox_a2=prior_tox_a2,
             prior_eff_a1=prior_eff_a1, 
             prior_eff_a2=prior_eff_a2,
             stage=stage, 
             c_tox=c_tox, 
             c_eff=c_eff,
             cohort_start=cohort_start,
             cohort=cohort, 
             pat_incl=pat_incl, 
             cdose1=cdose1,
             cdose2=cdose2, 
             startup_in=startup_in, 
             dose_adm1=dose_adm1, 
             dose_adm2=dose_adm2, 
             time_prog=time_prog,
             time_follow=time_follow, 
             toxicity=toxicity, 
             n_pat_comb=n_pat_comb,
             n_tox_comb=n_tox_comb,
             n_prog_comb=n_prog_comb,
             pi=pi,
             ptox_sup=ptox_sup, 
             resp=resp, 
             qeff_min=qeff_min,
             proba_tau=proba_tau,
             n_pat_comb_p=n_pat_comb_p,
             n_tox_comb_p=n_tox_comb_p,
             n_prog_comb_p=n_prog_comb_p,
             pi_p=pi_p,
             ptox_sup_p=ptox_sup_p, 
             resp_p=resp_p, 
             qeff_min_p=qeff_min_p,
             proba_tau_p=proba_tau_p)
             
  class(res) = "CombPlateau_next"

  return(res)
}




print.CombPlateau_next = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 
     
  print_rnd= function (hd, x) {cat(hd, "\n"); print(round(x, digits = dgt)); cat("\n")}
  print_rnd("Number of patients:" , x$n_pat_comb_p)
  print_rnd("Number of toxicities:", x$n_tox_comb_p)
  print_rnd("Number of progressions:", x$n_prog_comb_p)
  print_rnd("Toxicity prob:", x$pi_p)
  print_rnd("P(toxicity prob > tox_max):", x$ptox_sup_p)
  print_rnd("Efficacy prob:", x$resp_p)
  print_rnd("P(efficacy prob < eff_min):", x$qeff_min_p)
  print_rnd("Prob plateau: Agent 2", x$proba_tau_p)
  
  cat("Start-up phase ended:\t", x$startup_in,"\n")
  if(x$cdose1>0){
    if(x$stage==2){
      cat(paste("RECOMMENDED COMBINATION at the end of the trial:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
    else{
      cat(paste("NEXT RECOMMENDED COMBINATION:\t (",x$cdose1, ",", x$cdose2, ")\n",sep=""), sep="")
    }
  }
  else{
    cat(paste("THE DOSE-FINDING PROCESS SHOULD BE STOPPED WITHOUT COMBINATION RECOMMENDATION\n",sep=""), sep="")
  }
  cat("\n", "\n")
  
  cat("Number of patients included:\t", x$pat_incl, "\n")
  cat("Maximum toxicity:\t", x$tox_max, "\n")
  cat("Minimum efficacy:\t", x$eff_min, "\n")
  cat("Prior toxicity probabilities for agent 1:\n")
  print(round(x$prior_tox_a1, digits = dgt))
  cat("Prior toxicity probabilities for agent 2:\n")
  print(round(x$prior_tox_a2, digits = dgt))
  cat("Prior efficacy probabilities for agent 1:\n")
  print(round(x$prior_eff_a1, digits = dgt))
  cat("Prior efficacy probabilities for agent 2:\n")
  print(round(x$prior_eff_a2, digits = dgt))
  
  cat("Toxicity threshold:\t", x$c_tox, "\n")
  cat("Efficacy threshold:\t", x$c_eff, "\n")
  cat("Toxicity is binary \n")
  cat("Efficacy is a time-to-event \n")
  cat("Full follow-up time:\t", x$time_full, "\n")
  cat("Minimum waiting time between two dose cohorts is of one toxicity cycle of:\t", x$cycle, "\n") 

}



               
                    
                    
