powerLmer <-
function(varcomp,nval,alpha=0.05,nsim=100,ml=F) {
  if (missing(varcomp)) stop("Need the variance component vector")  #dam,sire,dxs,res
  if (missing(nval)) stop("Need the sample size vector") #dam,sire,off
  print(time1<- Sys.time()) #start time
dam<- NULL; sire<- NULL; famil<- NULL
damK<- varcomp[1];sireK<- varcomp[2];dxsK<- varcomp[3];resK<- varcomp[4]
damN<- nval[1];sireN<- nval[2];offN<- nval[3]
N<- damN*sireN*offN
sim_dat<- data.frame(dam_var=numeric(nsim), dam_pval=numeric(nsim), sire_var=numeric(nsim), sire_pval=numeric(nsim),
  dam.sire_var=numeric(nsim),dam.sire_pval=numeric(nsim),residual=numeric(nsim))
for (i in 1:nsim) {
  print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK))
  sireR<- rnorm(sireN,0,sd=sqrt(sireK))
  dxsR<- rnorm(damN*sireN,0,sd=sqrt(dxsK))
  resR<- rnorm(N,0,sd=sqrt(resK)) #noise
observ <- expand.grid(off=1:offN, dam=1:damN, sire=1:sireN)
observ$famil<- rep(1:(damN*sireN),each=offN)
observ<- within(observ, { resp<- damR[dam] + sireR[sire] + dxsR[famil] + resR } )
if (ml == F) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire), observ) }
if (ml == T) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire), observ, REML=F) }
comp<- data.frame(effect= as.data.frame(VarCorr(m))$grp[1:3],variance= colSums(diag(VarCorr(m)))[1:3])
sim_dat$dam_var[i]<- comp$variance[which(comp$effect=="dam")]
sim_dat$sire_var[i]<- comp$variance[which(comp$effect=="sire")]
sim_dat$dam.sire_var[i]<- comp$variance[which(comp$effect=="dam:sire")]
sim_dat$residual[i]<- attr(VarCorr(m),"sc")^2
p_rand<- randLmer(model=m,observ=observ)
sim_dat$dam_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam)")]
sim_dat$sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | sire)")]
sim_dat$dam.sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam:sire)")]
} #end loop
pwr_res<- data.frame(term=c("dam","sire","dam.sire","residual"), n= c(nval[c(1,2)],damN*sireN,NA),
  var_in= varcomp,
  var_out=c(mean(sim_dat$dam_var),mean(sim_dat$sire_var),mean(sim_dat$dam.sire_var),mean(sim_dat$residual)),
  power=c(sum(sim_dat$dam_pval < alpha)/nsim,sum(sim_dat$sire_pval < alpha)/nsim,
   sum(sim_dat$dam.sire_pval < alpha)/nsim,NA))
  print(Sys.time()- time1) #end time
  invisible(pwr_res)  #after time
}
