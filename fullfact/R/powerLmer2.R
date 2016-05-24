powerLmer2 <-
function(varcomp,nval,alpha=0.05,nsim=100,position=NULL,block=NULL,ml=F) {
  if (missing(varcomp)) stop("Need the variance component vector")  #dam,sire,dxs,res,(position),(block)
  if (missing(nval)) stop("Need the sample size vector") #dam,sire,off
  print(time1<- Sys.time()) #start time
dam<- NULL; sire<- NULL; famil<- NULL
damK<- varcomp[1];sireK<- varcomp[2];dxsK<- varcomp[3];resK<- varcomp[4]
damN<- nval[1];sireN<- nval[2];offN<- nval[3]
if (is.null(position) && is.null(block)) {
sim_dat<- data.frame(dam_var=numeric(nsim), dam_pval=numeric(nsim), sire_var=numeric(nsim), sire_pval=numeric(nsim),
  dam.sire_var=numeric(nsim),dam.sire_pval=numeric(nsim),residual=numeric(nsim))
for (i in 1:nsim) {
print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK));sireR<- rnorm(sireN,0,sd=sqrt(sireK));dxsR<- rnorm(damN*sireN,0,sd=sqrt(dxsK))
observ <- expand.grid(off=1:offN, dam=1:damN, sire=1:sireN)
observ$famil<- rep(1:(damN*sireN),each=offN)
N<- nrow(observ); resR<- rnorm(N,0,sd=sqrt(resK)) #noise
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
} #end nothing loop
pwr_res<- data.frame(term=c("dam","sire","dam.sire","residual"), n= c(nval[c(1,2)],damN*sireN,NA),
  var_in= varcomp,
  var_out=c(mean(sim_dat$dam_var),mean(sim_dat$sire_var),mean(sim_dat$dam.sire_var),mean(sim_dat$residual)),
  power=c(sum(sim_dat$dam_pval < alpha)/nsim,sum(sim_dat$sire_pval < alpha)/nsim,
   sum(sim_dat$dam.sire_pval < alpha)/nsim,NA))
} #end nothing
if (!is.null(position) && is.null(block)) { posK<- varcomp[5]; posN<- nval[4]
if (offN/position != posN) stop("Sample size of position (in nval) does not match offspring / position (in position)")
sim_dat<- data.frame(dam_var=numeric(nsim), dam_pval=numeric(nsim), sire_var=numeric(nsim), sire_pval=numeric(nsim),
  dam.sire_var=numeric(nsim),dam.sire_pval=numeric(nsim),pos_var=numeric(nsim),pos_pval=numeric(nsim),
  residual=numeric(nsim))
for (i in 1:nsim) {
print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK));sireR<- rnorm(sireN,0,sd=sqrt(sireK));dxsR<- rnorm(damN*sireN,0,sd=sqrt(dxsK))
  posR<- rnorm(posN,0,sd=sqrt(posK))
observ <- expand.grid(off=1:offN, dam=1:damN, sire=1:sireN)
observ$famil<- rep(1:(damN*sireN),each=offN)
observ$position<- rep(1:posN,each=offN/position)
N<- nrow(observ);resR<- rnorm(N,0,sd=sqrt(resK)) #noise
observ<- within(observ, { resp<- damR[dam] + sireR[sire] + dxsR[famil] + posR[position] + resR } )
  if (ml == F) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire) + (1|position), observ) }
  if (ml == T) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire) + (1|position), observ, REML=F) }
comp<- data.frame(effect= as.data.frame(VarCorr(m))$grp[1:4],variance= colSums(diag(VarCorr(m)))[1:4])
sim_dat$dam_var[i]<- comp$variance[which(comp$effect=="dam")]
sim_dat$sire_var[i]<- comp$variance[which(comp$effect=="sire")]
sim_dat$dam.sire_var[i]<- comp$variance[which(comp$effect=="dam:sire")]
sim_dat$pos_var[i]<- comp$variance[which(comp$effect=="position")]
sim_dat$residual[i]<- attr(VarCorr(m),"sc")^2
p_rand<- randLmer(model=m,observ=observ)
sim_dat$dam_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam)")]
sim_dat$sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | sire)")]
sim_dat$dam.sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam:sire)")]
sim_dat$pos_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | position)")]
} #end position only loop
pwr_res<- data.frame(term=c("dam","sire","dam.sire","position","residual"), n= c(nval[c(1,2)],damN*sireN,nval[4],NA),
  var_in= varcomp[c(1:3,5,4)],
  var_out=c(mean(sim_dat$dam_var),mean(sim_dat$sire_var),mean(sim_dat$dam.sire_var),mean(sim_dat$pos_var),
  mean(sim_dat$residual)), power=c(sum(sim_dat$dam_pval < alpha)/nsim,sum(sim_dat$sire_pval < alpha)/nsim,
  sum(sim_dat$dam.sire_pval < alpha)/nsim,sum(sim_dat$pos_pval < alpha)/nsim,NA))
} #end position only
if (is.null(position) && !is.null(block)) { blocK<- varcomp[5]; blocN<- nval[4]
if (damN != block[1]*blocN) stop("Sample size of dams does not match block design")
if (sireN != block[2]*blocN) stop("Sample size of sires does not match block design")
sim_dat<- data.frame(dam_var=numeric(nsim), dam_pval=numeric(nsim), sire_var=numeric(nsim), sire_pval=numeric(nsim),
  dam.sire_var=numeric(nsim),dam.sire_pval=numeric(nsim),bloc_var=numeric(nsim),bloc_pval=numeric(nsim),
  residual=numeric(nsim))
for (i in 1:nsim) {
print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK));sireR<- rnorm(sireN,0,sd=sqrt(sireK));dxsR<- rnorm(block[1]*block[2]*blocN,0,sd=sqrt(dxsK))
  blocR<- rnorm(blocN,0,sd=sqrt(blocK))
dam0<- stack(as.data.frame(matrix(1:(block[1]*blocN),ncol=blocN,nrow=block[1])))
sire0<- stack(as.data.frame(matrix(1:(block[2]*blocN),ncol=blocN,nrow=block[2])))
observ0<- merge(dam0,sire0, by="ind")
levels(observ0[,1])<- 1:blocN; colnames(observ0)<- c("block","dam","sire")
observ0$famil<- 1:nrow(observ0)  #add family
observ<- do.call("rbind", replicate(offN,observ0,simplify=F)) #expand
N<- nrow(observ);resR<- rnorm(N,0,sd=sqrt(resK)) #noise
observ<- within(observ, { resp<- damR[dam] + sireR[sire] + dxsR[famil] + blocR[block] + resR } )
  if (ml == F) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire) + (1|block), observ) }
  if (ml == T) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire) + (1|block), observ, REML=F) }
comp<- data.frame(effect= as.data.frame(VarCorr(m))$grp[1:4],variance= colSums(diag(VarCorr(m)))[1:4])
sim_dat$dam_var[i]<- comp$variance[which(comp$effect=="dam")]
sim_dat$sire_var[i]<- comp$variance[which(comp$effect=="sire")]
sim_dat$dam.sire_var[i]<- comp$variance[which(comp$effect=="dam:sire")]
sim_dat$bloc_var[i]<- comp$variance[which(comp$effect=="block")]
sim_dat$residual[i]<- attr(VarCorr(m),"sc")^2
p_rand<- randLmer(model=m,observ=observ)
sim_dat$dam_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam)")]
sim_dat$sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | sire)")]
sim_dat$dam.sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam:sire)")]
sim_dat$bloc_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | block)")]
} #end block only loop
pwr_res<- data.frame(term=c("dam","sire","dam.sire","block","residual"),
  n= c(nval[c(1,2)],block[1]*block[2]*blocN,nval[4],NA), var_in= varcomp[c(1:3,5,4)],
  var_out=c(mean(sim_dat$dam_var),mean(sim_dat$sire_var),mean(sim_dat$dam.sire_var),mean(sim_dat$bloc_var),
  mean(sim_dat$residual)), power=c(sum(sim_dat$dam_pval < alpha)/nsim,sum(sim_dat$sire_pval < alpha)/nsim,
  sum(sim_dat$dam.sire_pval < alpha)/nsim,sum(sim_dat$bloc_pval < alpha)/nsim,NA))
} #end block only
if (!is.null(position) && !is.null(block)) { posK<- varcomp[5]; blocK<- varcomp[6]; posN<- nval[4]; blocN<- nval[5]
if (offN/position != posN) stop("Sample size of position (in nval) does not match offspring / position (in position)")
if (damN != block[1]*blocN) stop("Sample size of dams does not match block design")
if (sireN != block[2]*blocN) stop("Sample size of sires does not match block design")
sim_dat<- data.frame(dam_var=numeric(nsim), dam_pval=numeric(nsim), sire_var=numeric(nsim), sire_pval=numeric(nsim),
  dam.sire_var=numeric(nsim),dam.sire_pval=numeric(nsim),pos_var=numeric(nsim),pos_pval=numeric(nsim),
  bloc_var=numeric(nsim),bloc_pval=numeric(nsim),residual=numeric(nsim))
for (i in 1:nsim) {
print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK));sireR<- rnorm(sireN,0,sd=sqrt(sireK))
  dxsR<- rnorm(block[1]*block[2]*blocN,0,sd=sqrt(dxsK));posR<- rnorm(posN,0,sd=sqrt(posK))
  blocR<- rnorm(blocN,0,sd=sqrt(blocK))
dam0<- stack(as.data.frame(matrix(1:(block[1]*blocN),ncol=blocN,nrow=block[1])))
sire0<- stack(as.data.frame(matrix(1:(block[2]*blocN),ncol=blocN,nrow=block[2])))
observ0<- merge(dam0,sire0, by="ind")
levels(observ0[,1])<- 1:blocN; colnames(observ0)<- c("block","dam","sire")
observ0$famil<- 1:nrow(observ0)  #add family
observ1<- do.call("rbind", replicate(position,observ0,simplify=F));rm(observ0) #expand for position
observ1$position<- sample(rep(1:posN,each=position)) #random assignment
observ<- do.call("rbind", replicate(offN,observ1,simplify=F)); rm(observ1) #expand for offspring
observ$off<- rep(1:offN,each=block[1]*block[2]*blocN)
N<- nrow(observ);resR<- rnorm(N,0,sd=sqrt(resK)) #noise
observ<- within(observ, { resp<- damR[dam] + sireR[sire] + dxsR[famil] + posR[position] + blocR[block] + resR } )
  if (ml == F) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire) + (1|position) + (1|block), observ) }
  if (ml == T) { m<- lmer(resp~ (1|dam) + (1|sire) + (1|dam:sire) + (1|position) + (1|block), observ, REML=F) }
comp<- data.frame(effect= as.data.frame(VarCorr(m))$grp[1:5],variance= colSums(diag(VarCorr(m)))[1:5])
sim_dat$dam_var[i]<- comp$variance[which(comp$effect=="dam")]
sim_dat$sire_var[i]<- comp$variance[which(comp$effect=="sire")]
sim_dat$dam.sire_var[i]<- comp$variance[which(comp$effect=="dam:sire")]
sim_dat$pos_var[i]<- comp$variance[which(comp$effect=="position")]
sim_dat$bloc_var[i]<- comp$variance[which(comp$effect=="block")]
sim_dat$residual[i]<- attr(VarCorr(m),"sc")^2
p_rand<- randLmer(model=m,observ=observ)
sim_dat$dam_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam)")]
sim_dat$sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | sire)")]
sim_dat$dam.sire_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | dam:sire)")]
sim_dat$pos_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | position)")]
sim_dat$bloc_pval[i]<- p_rand$p.value[which(p_rand$term=="(1 | block)")]
} #end position and block loop
pwr_res<- data.frame(term=c("dam","sire","dam.sire","position","block","residual"),
  n= c(nval[c(1,2)],block[1]*block[2]*blocN,nval[c(4,5)],NA), var_in= varcomp[c(1:3,5,6,4)],
  var_out=c(mean(sim_dat$dam_var),mean(sim_dat$sire_var),mean(sim_dat$dam.sire_var),mean(sim_dat$pos_var),
  mean(sim_dat$ bloc_var),mean(sim_dat$residual)), power=c(sum(sim_dat$dam_pval < alpha)/nsim,
  sum(sim_dat$sire_pval < alpha)/nsim,sum(sim_dat$dam.sire_pval < alpha)/nsim,sum(sim_dat$pos_pval < alpha)/nsim,
  sum(sim_dat$bloc_pval < alpha)/nsim,NA))
} #end position and block
  print(Sys.time()- time1) #end time
  invisible(pwr_res)  #after time
}
