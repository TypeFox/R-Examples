mtaBin_sim = function(ngroups=1, ndose, p_tox, p_eff, tox_max, eff_min, prior_tox, prior_eff, poisson_rate=1, n,
                   cohort_start=3, cohort=3, tite=TRUE, time_full=0, method="MTA-RA", s_1=function(n_cur){0.2}, 
                   s_2=0.07, cycle=0, nsim, c_tox=0.90, c_eff=0.40, seed=8, threads=0){
  val_s_1 = sapply(0:(n-1),s_1)
  n_ptox = length(p_tox)
  n_prior_tox = length(prior_tox)
  n_rate = length(poisson_rate)

  if(n_ptox != ndose){
    stop("The entered vector for true toxicity probabities is of wrong length.")
  }
  if(n_prior_tox != ndose){
    stop("The entered vector of initial guessed toxicity probabities is of wrong length.")
  }
  if(n_rate!= ngroups){
    stop("The entered vector of poisson parameters for groups inclusion is of wrong length.")
  }
  
  n_peff = length(p_eff)
  n_prior_eff = length(prior_eff)
  
  if( (is.matrix(p_eff) == FALSE && (n_peff != ndose || ngroups != 1)) ||
      (is.matrix(p_eff) == TRUE && (dim(p_eff)[1] != ngroups || dim(p_eff)[2] != ndose)) ){
    stop("The entered vector or matrix for true efficacy probabities is of wrong length.")
  }
  if( (is.matrix(prior_eff) == FALSE && (n_prior_eff != ndose || ngroups != 1)) ||
      (is.matrix(prior_eff) == TRUE && (dim(prior_eff)[1] != ngroups || dim(prior_eff)[2] != ndose)) ){
    stop("The entered vector or matrix of initial guessed efficacy probabities is of wrong length.")
  }

  ngroups = as.integer(ngroups)[1]
  ndose = as.integer(ndose)[1]
  p_tox = as.double(p_tox)
  p_eff = as.double(p_eff)
  tox_max = as.double(tox_max)[1]
  eff_min = as.double(eff_min)[1]
  prior_tox = as.double(prior_tox)
  prior_eff = as.double(prior_eff)
  poisson_rate = as.double(poisson_rate)
  n = as.integer(n)[1]
  tite = as.logical(tite)
  time_full = as.double(time_full)[1]
  cycle = as.double(cycle)[1]
  cohort_start = as.integer(cohort_start)[1]
  cohort = as.integer(cohort)[1]
  nsim = as.integer(nsim)[1]
  c_tox = as.double(c_tox)[1]
  c_eff = as.double(c_eff)[1]
  val_s_1 = as.double(val_s_1)
  s_2 = as.double(s_2)[1]
  seed = as.integer(seed)[1]
  threads = as.integer(threads)[1]

  for(i in 1:n_ptox){
    if(p_tox[i] < 0 || p_tox[i] > 1 || p_eff[i] < 0 || p_eff[i] > 1){
      stop("At least one of the true toxicity or efficacy probability is not comprised between 0 and 1.")
    }
  }
  for(i in 1:ndose){
    if(prior_tox[i] < 0 || prior_tox[i] > 1){
      stop("At least one of the initial guessed toxicity probability is not comprised between 0 and 1.")
    }
  }
  for(i in 1:(ndose*ngroups)){
    if(prior_eff[i] < 0 || prior_eff[i] > 1){
      stop("At least one of the initial guessed efficacy probability is not comprised between 0 and 1.")
    }
  }
  if(tox_max < 0 || tox_max > 1){stop("The toxicity upper bound is not comprised between 0 and 1.")}
  if(eff_min < 0 || eff_min > 1){stop("The efficacy lower bound is not comprised between 0 and 1.")}
  
  inconc = as.integer(numeric(ngroups))
  n_pat_dose = as.integer(numeric(ngroups*ndose))
  rec_dose = as.integer(numeric(ngroups*ndose))
  n_pat_tot = as.integer(0)
  n_tox = as.integer(numeric(ngroups*ndose))
  n_eff = as.integer(numeric(ngroups*ndose))
  tox_tot = as.integer(numeric(ngroups))
  eff_tot = as.integer(numeric(ngroups))
  n_pat_mtd = as.integer(numeric(ngroups))
  duration = as.double(0)
  
  # Appeler fonction C
  if(method=="MTA-RA"){
    mta = .C(C_dfmta_simu, tite=tite, TRUE, ngoups=ngroups, ndose=ndose, p_tox=p_tox, p_eff=p_eff, tox_max=tox_max, eff_min=eff_min,
                     prior_tox=prior_tox, prior_eff=prior_eff, poisson_rate=poisson_rate, n=n,
                     time_full=time_full, cycle=cycle, cohort_start=cohort_start, cohort=cohort, nsim=nsim, c_tox=c_tox, 
                     c_eff=c_eff, val_s_1=val_s_1, s_2=s_2, seed=seed, threads=threads,

                     inconc=inconc, n_pat_dose=n_pat_dose, rec_dose=rec_dose, n_pat_tot=n_pat_tot, n_tox=n_tox, n_eff=n_eff,
                     tox_tot=tox_tot, eff_tot=eff_tot, n_pat_mtd=n_pat_mtd, duration=duration)
  }
  else if(tite == TRUE && method=="MTA-PM"){
    mta = .C(C_dfmta_simu, tite=tite, FALSE, ngoups=ngroups, ndose=ndose, p_tox=p_tox, p_eff=p_eff, tox_max=tox_max, eff_min=eff_min,
                     prior_tox=prior_tox, prior_eff=prior_eff, poisson_rate=poisson_rate, n=n,
                     time_full=time_full, cycle=cycle, cohort_start=cohort_start, cohort=cohort, nsim=nsim, c_tox=c_tox, 
                     c_eff=c_eff, val_s_1=val_s_1, s_2=s_2, seed=seed, threads=threads,

                     inconc=inconc, n_pat_dose=n_pat_dose, rec_dose=rec_dose, n_pat_tot=n_pat_tot, n_tox=n_tox, n_eff=n_eff,
                     tox_tot=tox_tot, eff_tot=eff_tot, n_pat_mtd=n_pat_mtd, duration=duration)
  }
  else{stop("Method not found.")}

  nsim = mta$nsim
  
  inconc=mta$inconc/nsim*100
  n_pat_dose=mta$n_pat_dose/nsim
  rec_dose=mta$rec_dose/nsim*100
  n_pat_tot=mta$n_pat_tot/nsim
  n_tox=mta$n_tox/nsim
  n_eff=mta$n_eff/nsim
  tox_tot=mta$tox_tot/nsim
  eff_tot=mta$eff_tot/nsim
  n_pat_mtd=mta$n_pat_mtd/nsim
  duration=mta$duration/nsim

  res = list(call = match.call(),
             method=method,
             tite=tite,
             ngroups=ngroups,
             ndose=ndose,
             p_tox=p_tox,
             p_eff=matrix(p_eff, nrow=ngroups),
             tox_max=tox_max, 
             eff_min=eff_min,
             prior_tox=prior_tox, 
             prior_eff=matrix(prior_eff, nrow=ngroups),
             poisson_rate=poisson_rate,
             n=n,
             time_full=time_full, 
             cycle=cycle,
             cohort_start=cohort_start, 
             cohort=cohort, 
             nsim=nsim, 
             c_tox=c_tox, 
             c_eff=c_eff, 
             seed=seed,
             inconc=inconc, 
             n_pat_dose=matrix(n_pat_dose, nrow=ngroups), 
             rec_dose=matrix(rec_dose, nrow=ngroups), 
             n_pat_tot=n_pat_tot, 
             n_tox=matrix(n_tox, nrow=ngroups), 
             n_eff=matrix(n_eff, nrow=ngroups),
             duration=duration)
             
  class(res) = "mtaBin_sim"

  return(res)
}



print.mtaBin_sim = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 

  tab_res = rbind(x$p_tox, x$p_eff, x$prior_tox, x$prior_eff, x$rec_dose, x$n_pat_dose, x$n_tox, x$n_eff)
  
  dimnames(tab_res) = list(
      c("True toxicities", paste("True efficacies for group ",1:x$ngroups,sep=""), "Prior toxicities", paste("Prior efficacies for group ",1:x$ngroups,sep=""), 
                    paste("Percentage of Selection for group ",1:x$ngroups,sep=""), paste("Number of patients for group ",1:x$ngroups,sep=""), 
                    paste("Number of toxicities for group ",1:x$ngroups,sep=""), paste("Number of efficacies for group ",1:x$ngroups,sep="")),
      doses=1:x$ndose)
  print(round(tab_res, digits = dgt))   
  cat("\n") 
  cat(paste("Percentage of inconclusive trials for group ",1:x$ngroups,":\t",x$inconc,"\n",sep=""), sep="")
  cat("Allocation method:\t", x$method)
  cat("\n", "\n")              
  cat("Number of simulations:\t", x$nsim, "\n")
  cat("Total patients accrued:\t", x$n_pat_tot, "\n")
  cat("Toxicity upper bound:\t", x$tox_max, "\n")
  cat("Efficacy lower bound:\t", x$eff_min, "\n")
  cat(paste("Patient arrival for group ", 1:x$ngroups,  " is modeled as a Poisson process with rate: ", x$poisson_rate, " that is in mean ", 1/x$poisson_rate, 
      " patients during a full follow-up time\n", sep=""), sep="")
  cat("Toxicity threshold:\t", x$c_tox, "\n")
  cat("Efficacy threshold:\t", x$c_eff, "\n")
  cat("Cohort size start-up phase:\t", x$cohort_start, "\n")
  cat("Cohort size model phase:\t", x$cohort, "\n")
  if (x$tite) {
    cat("Efficicy is a time-to-event \n")
    cat("Trial mean duration:\t", x$duration, "\n")
    cat("Full follow-up time:\t", x$time_full, "\n")
    cat("Minimum waiting time between two dose cohorts is of one cycle of ", x$cycle, "\n") 
  }
  else{
    cat("Efficicy is not a time-to-event but binay \n")
  }
}


mtaBin_next = function(ngroups=1, group_cur=1, ndose, prior_tox, prior_eff, tox_max, eff_min, cohort_start, cohort, 
                    final=FALSE, method="MTA-RA", s_1=function(n_cur){0.2}, s_2=0.07, group_pat=rep(1,length(id_dose)), id_dose, toxicity, tite=TRUE, 
                    efficacy=0, time_follow=0, time_eff=0, time_full=0, cycle=0, c_tox=0.90, c_eff=0.40, seed = 8){
    
  # time_eff = time-to-efficacy with +infini if no efficacy 
  # Fill the efficacy for the group for which the estimation are not done or put NA  
  pat_incl_group = c()
  for(ng in 1:ngroups){
    pat_incl_group = c(pat_incl_group, length(which(group_pat==ng)))
  }
  pat_tot = sum(pat_incl_group)
  pat_group = which(group_pat==group_cur)
  group = c(group_pat, group_cur)
  if(pat_incl_group[group_cur] > 0) {
    cdose = id_dose[pat_group[pat_incl_group[group_cur]]]
  }
  else {
    cdose = 0
  }
  time_cur = 0
  time_incl = -time_follow
  val_s_1 = s_1(pat_incl_group[group_cur])
  if(tite == TRUE) {
    efficacy = as.numeric(time_eff < time_follow)
  }

  if(method=="MTA-RA"){
    ra_pm = TRUE
  }
  else if (method=="MTA-PM"){    
    ra_pm = FALSE
  }              
  else{
    stop("Method not found.")
  }  

  n_prior_tox = length(prior_tox)
  if(n_prior_tox != ndose){
    stop("The entered vector of initial guessed toxicity probabities is of wrong length.")
  }
  
  n_prior_eff = length(prior_eff)
  if(n_prior_eff != ndose){
    stop("The entered vector of initial guessed efficacy probabities is of wrong length.")
  }
  
  n_toxicity = length(toxicity)
  n_efficacy = length(efficacy)
  n_time_follow = length(time_follow)
  n_time_eff = length(time_eff)
  n_group_pat = length(group_pat)
  n_id_dose = length(id_dose)
  if(n_toxicity != pat_tot){
    stop("The entered vector of observed toxicities is of wrong length.")
  }
  if(n_id_dose != pat_tot){
    stop("The entered vector of administered dose levels is of wrong length.")
  }
  if(tite==FALSE && n_efficacy != pat_tot){
    stop("The entered vector of observed efficacies is of wrong length.")
  }
  if(tite==TRUE && n_time_follow != pat_tot){
    stop("The entered vector for patients' follow-up time is of wrong length.")
  }
  if(tite==TRUE && n_time_eff != pat_tot){
    stop("The entered vector for patients' time-to-efficacy is of wrong length.")
  }
  if(n_group_pat != pat_tot){
    stop("The entered vector for patients' group is of wrong length.")
  }
  
  ngroups = as.integer(ngroups)[1]
  ndose = as.integer(ndose)[1]
  cdose = as.integer(cdose-1)[1]
  prior_tox = as.double(prior_tox)
  prior_eff = as.double(prior_eff)
  tox_max = as.double(tox_max)[1]
  eff_min = as.double(eff_min)[1]
  cohort_start = as.integer(cohort_start)[1]
  cohort = as.integer(cohort)[1]
  final = as.logical(final)
  id_dose = as.integer(id_dose-1)
  toxicity = as.logical(toxicity)
  efficacy = as.logical(efficacy)
  tite = as.logical(tite)
  time_follow = as.double(time_follow)
  time_eff = as.double(time_eff)
  group = as.integer(group-1)
  time_full = as.double(time_full)[1]
  cycle = as.double(cycle)[1]
  pat_incl_group = as.integer(pat_incl_group)
  c_tox = as.double(c_tox)[1]
  c_eff = as.double(c_eff)[1]
  time_cur = as.double(time_cur)[1]
  val_s_1 = as.double(val_s_1)[1]
  s_2 = as.double(s_2)[1]
  seed = as.integer(seed)[1]
  in_startup = as.logical(numeric(1))
  
  if(group_cur < 1 || group_cur > ngroups){
    stop("The group number for which the next optimal dose should be estimated is not comprised between 1 and ngroups.")
  }
  for(i in 1:ndose){
    if(prior_tox[i] < 0 || prior_tox[i] > 1 || prior_eff[i] < 0 || prior_eff[i] > 1){
      stop("At least one of the initial guessed toxicity or efficacy probability is not comprised between 0 and 1.")
    }
  }
  if(tox_max < 0 || tox_max > 1){stop("The toxicity upper bound is not comprised between 0 and 1.")}
  if(eff_min < 0 || eff_min > 1){stop("The efficacy lower bound is not comprised between 0 and 1.")}
  

  pi = numeric(ndose)
  ptox_inf = numeric(ndose)
  resp = numeric(ndose)
  qeff_inf = numeric(ndose)
  proba_tau = numeric(ndose)
  
  mta = .C(C_dfmta_next, tite=tite, ra_pm=ra_pm , ngroups=ngroups, ndose=ndose, tox_max=tox_max, eff_min=eff_min,
           prior_tox=prior_tox, prior_eff=prior_eff, time_full=time_full, cycle=cycle, cohort_start=cohort_start, 
           cohort=cohort, c_tox=c_tox, c_eff=c_eff, val_s_1=val_s_1, s_2=s_2, seed=seed, cdose=cdose, 
           time_cur=time_cur, pat_incl_group=pat_incl_group, id_dose=id_dose, group=group, time_eff=time_eff, 
           time_incl=time_incl, efficacy=efficacy, toxicity=toxicity, final=final,
           in_startup=in_startup, pi=pi, ptox_inf=ptox_inf, resp=resp, qeff_inf=qeff_inf, proba_tau=proba_tau,
           NAOK=TRUE)
           
  pi=mta$pi
  ptox_inf=mta$ptox_inf
  resp=mta$resp
  qeff_inf=mta$qeff_inf 
  proba_tau=mta$proba_tau
  group=group+1
  id_dose=id_dose+1
  cdose=mta$cdose+1
  in_startup=mta$in_startup
           
  pat_tot = sum(pat_incl_group) 
  n_pat_dose_tot = numeric(ndose)
  n_tox_tot = numeric(ndose)
  n_eff = numeric(ndose)
  if(ngroups > 1){
    n_pat_dose = matrix(0, nrow=2, ncol=ndose)
    n_tox = matrix(0, nrow=2, ncol=ndose)
    for(i in 1:pat_tot){
      n_pat_dose_tot[id_dose[i]] = n_pat_dose_tot[id_dose[i]]+1
      n_pat_dose[ifelse(group_pat[i]==group_cur,1,2),id_dose[i]] = n_pat_dose[ifelse(group_pat[i]==group_cur,1,2),id_dose[i]]+1
      n_tox[ifelse(group_pat[i]==group_cur,1,2),id_dose[i]] = n_tox[ifelse(group_pat[i]==group_cur,1,2),id_dose[i]]+toxicity[i]
      n_tox_tot[id_dose[i]] = n_tox_tot[id_dose[i]]+toxicity[i]
      if(group_pat[i]==group_cur){
        n_eff[id_dose[i]] = n_eff[id_dose[i]]+efficacy[i]
      }
    }
  
    if(!in_startup){
      tab_tox = rbind(prior_tox, n_pat_dose_tot, n_tox[1,], n_tox[2,], pi, ptox_inf)
      tab_eff = rbind(prior_eff, n_pat_dose[1,], n_eff, resp, 1-qeff_inf, proba_tau)
      tab_res = rbind(tab_tox,tab_eff)
  
      dimnames(tab_res) = list(
          c("Prior toxicities", 
          "Total number of patients included",
          paste("Number of observed toxicities in group ", group_cur, sep=""), 
          "Number of observed toxicities in other groups", 
          "Estimated toxicity probabilities", 
          paste("P(toxicity probability < ", tox_max, ")", sep=""),    
          "Prior efficacies", 
          paste("Number of patients included in group ", group_cur, sep=""),
          paste("Number of observed efficacies in group ", group_cur, sep=""),
          "Estimated efficacy probabilities",
          paste("P(efficacy probability > to ", eff_min, ")", sep=""),
          paste("Plateau probabilities in group ", group_cur, sep="")),
          doses=1:ndose)
    }
    else{
      tab_res = rbind(prior_tox, prior_eff, n_pat_dose[1,], n_tox[1,], n_eff)
      dimnames(tab_res) = list(c("Prior toxicities", "Prior efficacies", 
      paste("Number of patients included in group ", group_cur, sep=""),
      paste("Number of observed toxicities in group ", group_cur, sep=""),
      paste("Number of observed efficacies in group ", group_cur, sep="")),
      doses=1:ndose)
    }
  }
  else{
    for(i in 1:pat_tot){
      n_pat_dose_tot[id_dose[i]] = n_pat_dose_tot[id_dose[i]]+1
      n_tox_tot[id_dose[i]] = n_tox_tot[id_dose[i]]+toxicity[i]
      n_eff[id_dose[i]] = n_eff[id_dose[i]]+efficacy[i]
    }
  
    if(!in_startup){
      tab_res = rbind(prior_tox, prior_eff, n_pat_dose_tot, n_tox_tot, pi, ptox_inf, n_eff, resp, 1-qeff_inf, proba_tau)
  
      dimnames(tab_res) = list(
          c("Prior toxicities", 
          "Prior efficacies",
          "Number of patients included",
          "Number of observed toxicities", 
          "Estimated toxicity probabilities", 
          paste("P(toxicity probability < ", tox_max, ")", sep=""),   
          "Number of observed efficacies",
          "Estimated efficacy probabilities",
          paste("P(efficacy probability > to ", eff_min, ")", sep=""),
          "Plateau probabilities"),
          doses=1:ndose)
    }
    else{
      tab_res = rbind(prior_tox, prior_eff, n_pat_dose_tot, n_tox_tot, n_eff)
      dimnames(tab_res) = list(c("Prior toxicities", "Prior efficacies", 
      "Number of patients included",
      "Number of observed toxicities",
      "Number of observed efficacies"),
      doses=1:ndose)
    }
  }

  res = list(call = match.call(),
             tite=tite, 
             ra_pm=ra_pm , 
             ngroups=ngroups, 
             ndose=ndose, 
             tox_max=tox_max, 
             eff_min=eff_min,
             prior_tox=prior_tox, 
             prior_eff=prior_eff, 
             time_full=time_full, 
             cycle=cycle, 
             cohort_start=cohort_start, 
             cohort=cohort, 
             final = final,
             c_tox=c_tox,
             c_eff=c_eff,
             seed=seed, 
             cdose=cdose, 
             in_startup=in_startup, 
             pat_incl_group=pat_incl_group, 
             id_dose=id_dose, 
             group_cur=group_cur,
             group_pat=group_pat,
             group=group, 
             time_eff=time_eff, 
             time_follow=time_follow,
             efficacy=efficacy, 
             toxicity=toxicity,
             pi=pi, 
             ptox_inf=ptox_inf, 
             resp=resp, 
             qeff_inf=qeff_inf, 
             proba_tau=proba_tau,
             pat_tot=pat_tot, 
             n_pat_dose_tot=n_pat_dose_tot,
             n_tox_tot=n_tox_tot,
             n_eff=n_eff,
             tab_res=tab_res)
             
  class(res) = "mtaBin_next"

  return(res)
}




print.mtaBin_next = function(x, dgt = 2, ...) {
  cat("Call:\n")
  print(x$call) 
   
  print(round(x$tab_res, digits = dgt))   
  cat("\n") 
  cat("Current Group for dose determination:\t", x$group_cur,"\n")
  cat("Start-up phase ended:\t", ifelse(x$in_startup, "NO", "YES"),"\n")
  if(x$cdose>0){
    if(x$final){
      cat(paste("RECOMMENDED DOSE at the end of the trial for group ",x$group_cur,":\t",x$cdose,"\n",sep=""), sep="")
    }
    else{
      cat(paste("NEXT RECOMMENDED DOSE for group ",x$group_cur,":\t",x$cdose,"\n",sep=""), sep="")
    }
  }
  else{
    cat(paste("THE DOSE-FINDING PROCESS SHOULD BE STOPPED for group ",x$group_cur," WITHOUT DOSE RECOMMENDATION\n",sep=""), sep="")
  }
  cat("\n", "\n")
  cat("Number of groups:\t", x$ngroups,"\n")
  cat("Number of patients included in group ", x$group_cur, ":\t", x$pat_incl_group[x$group_cur], "\n")
  cat("Total number of patients included:\t", x$pat_tot, "\n")
  cat("Maximum sample size reached:\t", x$final,"\n")
  cat("Allocation method:\t", ifelse(x$ra_pm==TRUE, "MTA-RA", "MTA-PM"),"\n")
               
  if(!x$in_startup){
    cat("Toxicity upper bound:\t", x$tox_max, "\n")
    cat("Efficacy lower bound:\t", x$eff_min, "\n")
    cat("Toxicity threshold:\t", x$c_tox, "\n")
    cat("Efficacy threshold:\t", x$c_eff, "\n")
  } 
  if (x$tite) {
    cat("Efficacy is a time-to-event \n")
    cat("Full follow-up time:\t", x$time_full, "\n")
    cat("Minimum waiting time between two dose cohorts is of one cycle of ", x$cycle, "\n") 
  }
  else{
    cat("Efficacy is binary \n")
  }
}
