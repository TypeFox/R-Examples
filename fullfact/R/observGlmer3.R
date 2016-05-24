observGlmer3 <-
function(observ,dam,sire,response,fam_link,remain,quasi=F,iter=1000) {  #change to 1000
  if (missing(observ)) stop("Need the observed data frame")
  if (missing(dam)) stop("Need the dam column name")
  if (missing(sire)) stop("Need the sire column name")
  if (missing(response)) stop("Need the response column name")
  if (missing(fam_link)) stop("Need the family(link) for the glmer")
  if (missing(remain)) stop("Need the remaining formula")
  print(time1<- Sys.time()) #start time
  if(paste(fam_link)[2]== "logit") { m_err<- (pi^2)/3 }
  if(paste(fam_link)[2]== "probit") { m_err<- 1 }
  if(paste(fam_link)[2]== "sqrt") { m_err<- 0.25 }
 if(quasi == F) {
  m<- glmer(formula= noquote(paste(response,"~ (1| ",dam,") + (1| ",sire,") + (1|",dam,":",sire,") +",remain,sep="")),
    family=fam_link,data=observ) }
 if(quasi == T) {
  observ$dispersion<- as.factor(1:length(observ[,1]))
  m<- glmer(formula=
    noquote(paste(response,"~ (1| ",dam,") + (1| ",sire,") + (1|",dam,":",sire,") + (1|dispersion) + ",remain,sep="")),
    family=fam_link,data=observ) }
  if(paste(fam_link)[2]== "log") {
    #Get random effects names to generate null model
    rand.formula <- reformulate(sapply(findbars(formula(m)),function(x) paste0("(", deparse(x), ")")),response=".")
    #Generate null model (intercept and random effects only, no fixed effects)
    null.m <- update(m, rand.formula)
    m_err<- log(1/exp(fixef(null.m)[1]) + 1) }
  var_fixed <- var(as.vector(fixef(m) %*% t(m@pp$X)))  #total variation of fixed effects
  var_rand<- sum(as.data.frame(VarCorr(m))$vcov)
  tot<- sum(as.data.frame(VarCorr(m))$vcov,m_err,var_fixed)
  comp<- data.frame(effect= c(as.data.frame(VarCorr(m))$grp),effect2=c(as.data.frame(VarCorr(m))$var1),
    variance= c(as.data.frame(VarCorr(m))$vcov))
  other<- data.frame(component= c("Residual","Total"), variance=c(as.numeric(m_err),tot))
  comp$percent<- 100*comp$variance/tot; other$percent<- 100*other$variance/tot
 if (length(fixef(m)) > 1) {
   fix_eff<- data.frame(effect=c(paste(attributes(terms(m))$term.labels),"Fix_Tot"),
    variance=c(rep(NA,length(attributes(terms(m))$term.labels)),var_fixed))
   fix_eff$percent<- 100*fix_eff$variance/tot
   fix_eff$Chi.sq<- NA;fix_eff$p.value<- NA
 #afex package
   p_fixed<- mixed(m,data=observ,method = "PB",args.test = list(nsim = iter))$anova_table[c(1,4)]
 #fixed term matching
   f_term<- paste(attributes(terms(m))$term.labels)
 for (i in  1:length(f_term)) {
   fix_eff[,c(4,5)][which(fix_eff$effect==f_term[i]),]<- p_fixed[which(rownames(p_fixed)==f_term[i]),]   } }
   p_rand<- randGlmer(model=m,observ=observ,fam_link=fam_link)  #internal function
   p_fix<- fixedGlmer(model=m,observ=observ,fam_link=fam_link)
  levels(comp$effect)[which(levels(comp$effect)==dam)]<- "dam"
  levels(comp$effect)[which(levels(comp$effect)==sire)]<- "sire"
  levels(comp$effect)[which(levels(comp$effect)==noquote(paste(dam,":",sire,sep="")))]<- "dam:sire"
  levels(comp$effect)[which(levels(comp$effect)==noquote(paste(dam,".",sire,sep="")))] <- "dam:sire"  #weird period sometimes
  comp2<- data.frame(component=c("additive","nonadd","maternal"),variance=0,percent=0)
  comp2$variance<- c(4*comp$variance[which(comp$effect=="sire")],4*comp$variance[which(comp$effect=="dam:sire")],
    comp$variance[which(comp$effect=="dam")]- comp$variance[which(comp$effect=="sire")])
  comp2$percent<- 100*comp2$variance/tot
  if (length(fixef(m)) == 1) { var_obj<- list(random=comp,LRT=p_rand,other=other,calculation=comp2) }
  if (length(fixef(m)) > 1) { var_obj<- list(fixed=fix_eff,LRT.fixed=p_fix,random=comp,LRT.random=p_rand,other=other,calculation=comp2) }
   print(Sys.time()- time1) #end time
   invisible(var_obj)  #after time
}
