JackLmer3 <-
function(observ,dam,sire,response,remain,ml=F,size=1,first=NULL) {
  if (missing(observ)) stop("Need the observed data frame")
  if (missing(dam)) stop("Need the dam column name")
  if (missing(sire)) stop("Need the sire column name")
  if (missing(response)) stop("Need the response column name")
  if (missing(remain)) stop("Need the remaining formula")
  if (size %% 1 != 0) stop("Size needs to be a positive integer")
  if (size %% 1 == 0 && size < 0) stop("Size needs to be a positive integer")
  if (!is.null(first) && first %% 1 != 0) stop("First needs to be a positive integer")
  if (!is.null(first) && first %% 1 == 0 && first < 0) stop("First needs to be a positive integer")
  print(time1<- Sys.time()) #start time
if (size == 1) {
  var_fixed<- matrix(0,ncol=1,nrow=nrow(observ))   #variance of fixed effects
  rand_l<- nchar(remain)- nchar(gsub(")","",remain))  #number of random effects
  var_rand<- matrix(0,ncol=rand_l+4,nrow=nrow(observ))   #variance of random effects (no total)
 if (is.null(first)) { iter<- nrow(observ) }
 if (!is.null(first)) { iter<- first }
  for (i in 1:iter) {
    print(paste0("Removing observation: ",i," of ",length(observ[,1])))
if (ml == F) {
  m<- lmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,") +",remain,sep="")),data=observ[-i,]) }
if (ml == T) {
  m<- lmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,") +",remain,sep="")),data=observ[-i,],REML=F) }
  var_fixed[i,] <- var(as.vector(fixef(m) %*% t(m@pp$X)))  #total variation of fixed effects
  var_rand[i,]<- c(colSums(diag(VarCorr(m))),attr(VarCorr(m),"sc")^2)
  col_names<- as.data.frame(VarCorr(m))$grp  } #end loop
if (sum(var_fixed[,1]) != 0) { jack<- as.data.frame(cbind(var_rand,var_fixed))
  colnames(jack)<- c(col_names,"Fixed") }
if (sum(var_fixed[,1]) == 0) { jack- as.data.frame(var_rand)
  colnames(jack)<- col_names }
jack$Total<- rowSums(jack)
} #end single Jack
if (size > 1) {
 observS1 <- observ[sample(nrow(observ)),];rm(observ) #shuffle once
 observS2 <- observS1[sample(nrow(observS1)),];rm(observS1) #shuffle twice
 n_mod<- floor(nrow(observS2)/size)  #number of models in data frame
  var_fixed<- matrix(0,ncol=1,nrow=n_mod)   #variance of fixed effects
  rand_l<- nchar(remain)- nchar(gsub(")","",remain))  #number of random effects
  var_rand<- matrix(0,ncol=rand_l+4,nrow=n_mod)   #variance of random effects (no total)
 if (is.null(first)) { iter<- n_mod }
 if (!is.null(first)) { iter<- first }
 beg<- matrix(0,ncol=1,nrow=n_mod)
 beg[1]<- 1 #always
 for (i in 2:n_mod) { beg[i]<- beg[i-1] + size }
 finn<- matrix(0,ncol=1,nrow=n_mod)
 finn[1]<- size #always
 for (i in 2:n_mod) { finn[i]<- finn[i-1] + size }
 for (i in 1:iter) {
  print(paste0("Removing block: ",i," of ",n_mod))
if (ml == F) {
  m<- lmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,") +",remain,sep="")),data=observS2[-c(beg[i]:finn[i]),]) }
if (ml == T) {
  m<- lmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,") +",remain,sep="")),data=observS2[-c(beg[i]:finn[i]),],REML=F) }
  var_fixed[i,] <- var(as.vector(fixef(m) %*% t(m@pp$X)))  #total variation of fixed effects
  var_rand[i,]<- c(colSums(diag(VarCorr(m))),attr(VarCorr(m),"sc")^2)
  col_names<- as.data.frame(VarCorr(m))$grp  } #end loop
if (sum(var_fixed[,1]) != 0) { jack<- as.data.frame(cbind(var_rand,var_fixed))
  colnames(jack)<- c(col_names,"Fixed") }
if (sum(var_fixed[,1]) == 0) { jack- as.data.frame(var_rand)
  colnames(jack)<- col_names }
jack$Total<- rowSums(jack)
}  #end block Jack
  colnames(jack)[which(colnames(jack)==dam)]<- "dam"
  colnames(jack)[which(colnames(jack)==sire)]<- "sire"
  colnames(jack)[which(colnames(jack)==noquote(paste(dam,":",sire,sep="")))]<- "dam:sire"
  jack$additive<- 4*jack$sire
  jack$nonadd<- 4*jack$'dam:sire'
  jack$maternal<- jack$dam- jack$sire
 print(Sys.time()- time1) #end time
 invisible(jack)  #after time
}
