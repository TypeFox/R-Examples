powerGlmer <-
function(varcomp,nval,fam_link,alpha=0.05,nsim=100,poisLog=NULL) {
  if (missing(varcomp)) stop("Need the variance component vector")  #dam,sire,dxs
  if (missing(nval)) stop("Need the sample size vector") #dam,sire,off
  if (missing(fam_link)) stop("Need the family(link) for the glmer")
  if(paste(fam_link)[2]== "log" && is.null(poisLog)) stop("Need the poisLog variance component")
  print(time1<- Sys.time()) #start time
dam<- NULL; sire<- NULL; famil<- NULL
damK<- varcomp[1];sireK<- varcomp[2];dxsK<- varcomp[3]
  if(paste(fam_link)[2]== "logit") { resK<- (pi^2)/3 }
  if(paste(fam_link)[2]== "probit") { resK<- 1 }
  if(paste(fam_link)[2]== "sqrt") { resK<- 0.25 }
  if(paste(fam_link)[2]== "log") { resK<- poisLog }
damN<- nval[1];sireN<- nval[2];offN<- nval[3]
N<- damN*sireN*offN
sim_dat<- data.frame(dam_var=numeric(nsim), dam_pval=numeric(nsim), sire_var=numeric(nsim), sire_pval=numeric(nsim),
  dam.sire_var=numeric(nsim),dam.sire_pval=numeric(nsim),residual=numeric(nsim))
for (i in 1:nsim) {
  print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK))
  sireR<- rnorm(sireN,0,sd=sqrt(sireK))
  dxsR<- rnorm(damN*sireN,0,sd=sqrt(dxsK))
observ <- expand.grid(off=1:offN, dam=1:damN, sire=1:sireN)
observ$famil<- rep(1:(damN*sireN),each=offN)
observ<- within(observ, { resp<- damR[dam] + sireR[sire] + dxsR[famil] } ) #do not add resR
  if(paste(fam_link)[2]== "logit") { observ<- within(observ, {resp2<- rbinom(nrow(observ),prob=plogis(resp),size=1)}) }
  if(paste(fam_link)[2]== "probit") { observ<- within(observ, {resp2<- rbinom(nrow(observ),prob=pnorm(resp),size=1)}) } 
  if(paste(fam_link)[2]== "sqrt") { observ<- within(observ, {resp2<- rpois(nrow(observ),resp^2)}) }
  if(paste(fam_link)[2]== "log") { observ<- within(observ, {resp2<- rpois(nrow(observ),exp(resp))}) }  
m<- glmer(resp2~ (1|dam) + (1|sire) + (1|dam:sire), family=fam_link, data=observ)
comp<- data.frame(effect= as.data.frame(VarCorr(m))$grp,variance= colSums(diag(VarCorr(m))))
sim_dat$dam_var[i]<- comp$variance[which(comp$effect=="dam")]
sim_dat$sire_var[i]<- comp$variance[which(comp$effect=="sire")]
sim_dat$dam.sire_var[i]<- comp$variance[which(comp$effect=="dam:sire")]
if(paste(fam_link)[2]!= "log") { sim_dat$residual[i]<- resK }
if(paste(fam_link)[2]== "log") { 
    rand.formula <- reformulate(sapply(findbars(formula(m)),function(x) paste0("(", deparse(x), ")")),response=".")
    null.m <- update(m, rand.formula)
    sim_dat$residual[i]<- log(1/exp(fixef(null.m)[1]) + 1) }
p_rand<- randGlmer(model=m,observ=observ,fam_link=fam_link)
sim_dat$dam_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam)")]
sim_dat$sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | sire)")]
sim_dat$dam.sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam:sire)")]
} #end loop
 pwr_res<- data.frame(term=c("dam","sire","dam.sire","residual"), n= c(nval[c(1,2)],damN*sireN,NA),
  var_in= c(varcomp,resK),
  var_out=c(mean(sim_dat$dam_var),mean(sim_dat$sire_var),mean(sim_dat$dam.sire_var),mean(sim_dat$residual)),
  power=c(sum(sim_dat$dam_pval < alpha)/nsim,sum(sim_dat$sire_pval < alpha)/nsim,
  sum(sim_dat$dam.sire_pval < alpha)/nsim,NA))
  print(Sys.time()- time1) #end time
  invisible(pwr_res)  #after time
}
