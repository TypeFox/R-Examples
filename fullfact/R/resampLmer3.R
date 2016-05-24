resampLmer3 <-
function(resamp,dam,sire,response,start,end,remain,ml=F) {
  if (missing(resamp)) stop("Need the resampled data frame")
  if (missing(dam)) stop("Need the dam column name")
  if (missing(sire)) stop("Need the sire column name")
  if (missing(response)) stop("Need the response column name")
  if (missing(start)) stop("Need the starting model number")
  if (missing(end)) stop("Need the ending model number")
  if (missing(remain)) stop("Need the remaining model formula with #")
  print(time1<- Sys.time()) #start time
  response2<- colnames(resamp)[grep(paste(response), colnames(resamp))]
  dam2<- colnames(resamp)[grep(paste(dam), colnames(resamp))]
  sire2<- colnames(resamp)[grep(paste(sire), colnames(resamp))]
  mdl<- matrix(0,ncol=1,nrow=length(response2))
    for (i in 1:length(response2)) { mdl[i,]<- gsub("#",i,remain)  }
  var_fixed<- matrix(0,ncol=1,nrow=length(response2))   #variance of fixed effects
    rand_l<- nchar(remain)- nchar(gsub(")","",remain))  #number of random effects
    var_rand<- matrix(0,ncol=rand_l+4,nrow=length(response2))   #variance of random effects + 4 constant
  for (j in start:end) {
    print(paste("Working on model: ", j, sep=""))
 if (ml == F) {
  m<- lmer(formula= noquote(paste(response2[j],"~ (1| ",dam2[j],") + (1| ",sire2[j],") + (1|",dam2[j],":",sire2[j],") +",mdl[j,],sep="")),
    data=resamp) }
 if (ml == T) {
  m<- lmer(formula= noquote(paste(response2[j],"~ (1| ",dam2[j],") + (1| ",sire2[j],") + (1|",dam2[j],":",sire2[j],") +",mdl[j,],sep="")),
    data=resamp, REML=F) }
  var_fixed[j,] <- var(as.vector(fixef(m) %*% t(m@pp$X)))  #total variation of fixed effects
  var_rand[j,]<- c(colSums(diag(VarCorr(m))),attr(VarCorr(m),"sc")^2)
  col_names<- as.data.frame(VarCorr(m))$grp; rm(m)  }
if (sum(var_fixed[,1]) != 0) { comp<- as.data.frame(cbind(var_rand,var_fixed))
  colnames(comp)<- c(col_names,"Fixed") }
if (sum(var_fixed[,1]) == 0) { comp<- as.data.frame(var_rand)
  colnames(comp)<- col_names }
comp$Total<- rowSums(comp)
colnames(comp)<- gsub(end,'', colnames(comp))
  colnames(comp)[which(colnames(comp)==dam)]<- "dam"
  colnames(comp)[which(colnames(comp)==sire)]<- "sire"
  colnames(comp)[which(colnames(comp)==noquote(paste(dam,":",sire,sep="")))]<- "dam:sire"
  colnames(comp)[which(colnames(comp)==noquote(paste(dam,".",sire,sep="")))]<- "dam:sire" #weird period sometimes
  comp$additive<- 4*comp$sire
  comp$nonadd<- 4*comp$'dam:sire'
  comp$maternal<- comp$dam- comp$sire
   print(Sys.time()- time1) #end time
   invisible(comp)  #after time
}
