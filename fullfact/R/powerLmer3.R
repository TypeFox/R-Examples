powerLmer3 <-
function(var_rand,n_rand,design,remain,var_fix=NULL,n_fix=NULL,alpha=0.05,nsim=100,ml=F,ftest="LR",iter=NULL) {
  if (missing(var_rand)) stop("Need the random variance component vector")  #dam,sire,dxs,res,(others)
  if (missing(n_rand)) stop("Need the random sample size vector") #dam,sire,family,(random others)
  if (missing(design)) stop("Need the design data frame") #dam,sire,family,(random then fixed others)
  if (missing(remain)) stop("Need the remaining formula")
  if (!is.null(var_fix) && is.null(n_fix)) stop("Need the fixed sample size vector")
  if (is.null(var_fix) && !is.null(n_fix)) stop("Need the fixed variance component vector")
  if (ftest=="PB" && is.null(iter)) stop("Need the iter for ftest PB")
  print(time1<- Sys.time()) #start time
damK<- var_rand[1];sireK<- var_rand[2];dxsK<- var_rand[3];resK<- var_rand[4]
randK<- rep(0,length(var_rand)-4); for (i in 1:length(randK)) { randK[i]<- var_rand[4+i] }
damN<- n_rand[1];sireN<- n_rand[2];famN<- n_rand[3]
randN<- rep(0,length(n_rand)-3); for (i in 1:length(randN)) { randN[i]<- n_rand[3+i] }
rNames0<- colnames(design)[-c(1:3)]; rNames<- rNames0[1:length(randN)]
randWith0<- list(); for (w in 1:length(randN)) { randWith0[[w]]<- paste0("randR[[",w,"]][",rNames[w],"]") }
randWith<- noquote(paste(unlist(randWith0),"+ ",collapse=""))
sim_varR<- matrix(0,ncol=length(var_rand),nrow=nsim) #+ residual
sim_pvalR<- matrix(0,ncol=length(n_rand),nrow=nsim) #NA residual
if (!is.null(var_fix)) {
  fixK<- rep(0,length(var_fix)); for (i in 1:length(fixK)) { fixK[i]<- var_fix[i] }
  fixN<- rep(0,length(n_fix)); for (i in 1:length(fixN)) { fixN[i]<- n_fix[i] }
  fNames<- colnames(design)[-c(1:(3+length(randN)))]
fixWith0<- list(); for (w in 1:length(n_fix)) { fixWith0[[w]]<- paste0("fixR[[",w,"]][",fNames[w],"]") }
fixWith<- noquote(paste(unlist(fixWith0),"+ ",collapse=""))
sim_varF<- numeric(nsim) #total fix var in a column
sim_pvalF<- matrix(0,ncol=length(n_fix),nrow=nsim) #p-value each effect
} #end fixed effects
for (i in 1:nsim) {
  print(paste0("Starting simulation: ", i))
  damR<- rnorm(damN,0,sd=sqrt(damK));sireR<- rnorm(sireN,0,sd=sqrt(sireK));dxsR<- rnorm(famN,0,sd=sqrt(dxsK))
  randR<- list(); for (r in 1:length(randN)) { randR[[r]]<- rnorm(randN[r],0,sd=sqrt(randK[r])) }
  resR<- rnorm(nrow(design),0,sd=sqrt(resK)) #noise
observ<- design
if (is.null(var_fix)) { observ<- within(observ, {
 eval(parse(text=paste0("resp<- damR[dam] + sireR[sire] + dxsR[family] + ",randWith,"resR"))) } ) }
if (!is.null(var_fix)) {
 fixR<- list(); for (r in 1:length(fixN)) {
  if (fixN[r] == 1) { fixR[[r]]<- sort(rnorm(nrow(observ),0,sd=sqrt(fixK[r]))) }  #continous, sort ascending
  if (fixN[r] > 1) { fixR[[r]]<- rnorm(fixN[r],0,sd=sqrt(fixK[r])) }  #factor
} #end fixR
 observ<- within(observ, {
 eval(parse(text=paste0("resp<- damR[dam] + sireR[sire] + dxsR[family] + ",randWith, fixWith,"resR"))) } ) }
if (ml == F) { m<- lmer(formula= noquote(paste0("resp~ (1|dam) + (1|sire) + (1|dam:sire) +", remain)), observ) }
if (ml == T) { m<- lmer(formula= noquote(paste0("resp~ (1|dam) + (1|sire) + (1|dam:sire) +", remain)), observ, REML=F) }
n<- length(as.data.frame(VarCorr(m))$vcov)- 1 #minus residual
comp<- data.frame(effect= as.data.frame(VarCorr(m))$grp[1:n], effect2=as.data.frame(VarCorr(m))$var1[1:n],
  variance= as.data.frame(VarCorr(m))$vcov[1:n])
sim_varR[,1][i]<- comp$variance[which(comp$effect=="dam")]
sim_varR[,2][i]<- comp$variance[which(comp$effect=="sire")]
sim_varR[,3][i]<- comp$variance[which(comp$effect=="dam:sire")]
for (v in 1:length(randN)) { sim_varR[,(3+v)][i]<- comp$variance[which(comp$effect==rNames[v])] }
sim_varR[,(length(comp$effect)+1)][i]<- attr(VarCorr(m),"sc")^2
p_rand<- randLmer(model=m,observ=observ)
sim_pvalR[,1][i]<- p_rand$p.value[which(p_rand$term=="(1 | dam)")]
sim_pvalR[,2][i]<- p_rand$p.value[which(p_rand$term=="(1 | sire)")]
sim_pvalR[,3][i]<- p_rand$p.value[which(p_rand$term=="(1 | dam:sire)")]
for (p in 1:length(randN)) { sim_pvalR[,(3+p)][i]<- p_rand$p.value[which(p_rand$term==paste0("(1 | ",rNames[p],")"))] }
if (!is.null(var_fix)) {
  sim_varF[i]<- var(as.vector(fixef(m) %*% t(m@pp$X)))
if (ftest=="LR") { f_rand<- fixedLmer(model=m,observ=observ)$p.value } #p-value
if (ftest=="PB") { f_rand<- mixed(m,data=observ,method = "PB",args.test = list(nsim = iter))$anova_table[4] }
  sim_pvalF[i,]<- as.data.frame(f_rand)[,1]  }  #end fixed effects
} #end loop
pwr_rand<- matrix(0,ncol=5,nrow=length(var_rand))
pwr_rand[,1]<- c("dam","sire","dam.sire",rNames,"residual")
pwr_rand[,2]<- c(n_rand,NA)
pwr_rand[,3]<- c(var_rand[1:3],var_rand[-c(1:4)],var_rand[4])
for(x in 1:length(var_rand)) { pwr_rand[,4][x]<- mean(sim_varR[,x]) }
for(y in 1:length(n_rand)) { pwr_rand[,5][y]<- sum(sim_pvalR[,y] < alpha)/nsim }
pwr_rand[,5][length(n_rand)+1]<- NA            #residual
pwr_rand<- as.data.frame(pwr_rand); colnames(pwr_rand)<- c("term","n","var_in","var_out","power")
pwr_res<- pwr_rand
if (!is.null(var_fix)) {
pwr_fix1<- matrix(0,ncol=3,nrow=1)
  pwr_fix1[,1]<- "fix_eff"; pwr_fix1[,2]<- sum(var_fix); pwr_fix1[,3]<- mean(sim_varF)
  pwr_fix1<- as.data.frame(pwr_fix1); colnames(pwr_fix1)<- c("group","var_in","var_out")
pwr_fix2<- matrix(0,ncol=3,nrow=length(n_fix))
  pwr_fix2[,1]<- fNames; pwr_fix2[,2]<- n_fix
  for(y in 1:length(n_fix)) { pwr_fix2[,3][y]<- sum(sim_pvalF[,y] < alpha)/nsim }
  pwr_fix2<- as.data.frame(pwr_fix2); colnames(pwr_fix2)<- c("term","n","power")
  pwr_res<- list(group=pwr_fix1, fixed=pwr_fix2, random=pwr_rand)  } #end fixed
  print(Sys.time()- time1) #end time
  invisible(pwr_res)  #after time
}
