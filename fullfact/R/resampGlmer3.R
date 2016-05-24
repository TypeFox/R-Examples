resampGlmer3 <-
function(resamp,dam,sire,response,fam_link,start,end,remain,quasi=F) {
  if (missing(resamp)) stop("Need the resampled data frame")
  if (missing(dam)) stop("Need the dam column name")
  if (missing(sire)) stop("Need the sire column name")
  if (missing(response)) stop("Need the response column name")
  if (missing(fam_link)) stop("Need the family(link) for the glmer")
  if (missing(start)) stop("Need the starting model number")
  if (missing(end)) stop("Need the ending model number")
  if (missing(remain)) stop("Need the remain model formula with #")
  print(time1<- Sys.time()) #start time
  if(paste(fam_link)[2]== "logit") { m_err<- (pi^2)/3 }
  if(paste(fam_link)[2]== "probit") { m_err<- 1 }
  if(paste(fam_link)[2]== "sqrt") { m_err<- 0.25 }
  response2<- colnames(resamp)[grep(paste(response), colnames(resamp))]
  dam2<- colnames(resamp)[grep(paste(dam), colnames(resamp))]
  sire2<- colnames(resamp)[grep(paste(sire), colnames(resamp))]
  mdl<- matrix(0,ncol=1,nrow=length(response2))
    for (i in 1:length(response2)) { mdl[i,]<- gsub("#",i,remain)  }
  var_fixed<- matrix(0,ncol=1,nrow=length(response2))   #variance of fixed effects
 if (quasi == F) {
    rand_l<- nchar(remain)- nchar(gsub(")","",remain))  #number of random effects
    var_rand<- matrix(0,ncol=rand_l+3,nrow=length(response2)) }  #variance of random effects + 3 constant
 if (quasi == T) {
    rand_l<- nchar(remain)- nchar(gsub(")","",remain)) + 1  #number of random effects  + quasi
    var_rand<- matrix(0,ncol=rand_l+3,nrow=length(response2))
    resamp$dispersion<- as.factor(1:length(resamp[,1])) }
  res<- matrix(0,ncol=1,nrow=length(response2))   #for poisson(log)
  for (j in start:end) {
    print(paste("Working on model: ", j, sep=""))
 if (quasi == F) {
  m<- glmer(formula=
    noquote(paste(response2[j],"~ (1| ",dam2[j],") + (1| ",sire2[j],") + (1|",dam2[j],":",sire2[j],") +",mdl[j,],sep="")),
    family=fam_link,data=resamp) }
  if (quasi == T) {
  m<- glmer(formula=
    noquote(paste(response2[j],"~ (1| ",dam2[j],") + (1| ",sire2[j],") + (1|",dam2[j],":",sire2[j],") + (1|dispersion) +",mdl[j,],sep="")),
    family=fam_link,data=resamp) }
   if(paste(fam_link)[2]== "log") {
    #Get random effects names to generate null model
    rand.formula <- reformulate(sapply(findbars(formula(m)),function(x) paste0("(", deparse(x), ")")),response=".")
    #Generate null model (intercept and random effects only, no fixed effects)
    null.m <- update(m, rand.formula)
    res[j,]<- log(1/exp(fixef(null.m)[1]) + 1) }
  var_fixed[j,] <- var(as.vector(fixef(m) %*% t(m@pp$X)))  #total variation of fixed effects
  var_rand[j,]<- colSums(diag(VarCorr(m)))
  col_names<- as.data.frame(VarCorr(m))$grp; rm(m)  }
if (sum(var_fixed[,1]) != 0) {
  if(paste(fam_link)[2]!= "log") { comp<- as.data.frame(cbind(var_rand,m_err,var_fixed)) }
  if(paste(fam_link)[2]== "log") { comp<- as.data.frame(cbind(var_rand,res,var_fixed)) }
  colnames(comp)<- c(col_names,"Residual","Fixed") }
if (sum(var_fixed[,1]) == 0) {
  if(paste(fam_link)[2]!= "log") { comp<- as.data.frame(cbind(var_rand,m_err)) }
  if(paste(fam_link)[2]== "log") { comp<- as.data.frame(cbind(var_rand,res)) }
  colnames(comp)<- c(col_names,"Residual") }
comp$Total<- rowSums(comp)
colnames(comp)<- gsub(end,'', colnames(comp))
  colnames(comp)[which(colnames(comp)==dam)]<- "dam"
  colnames(comp)[which(colnames(comp)==sire)]<- "sire"
  colnames(comp)[which(colnames(comp)==noquote(paste(dam,":",sire,sep="")))]<- "dam:sire"
  colnames(comp)[which(colnames(comp)==noquote(paste(dam,".",sire,sep="")))]<- "dam:sire" #weird period sometimes
  comp$additive<- 4*comp$sire
  comp$nonadd<- 4*comp$'dam:sire'
  comp$maternal<- comp$dam- comp$sire
   print(Sys.time()- time1) #end time, keep no quote in one line
   invisible(comp)  #after time
}
