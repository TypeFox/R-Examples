JackGlmer <-
function(observ,dam,sire,response,fam_link,quasi=F,size=1,first=NULL) {
  if (missing(observ)) stop("Need the observed data frame")
  if (missing(dam)) stop("Need the dam column name")
  if (missing(sire)) stop("Need the sire column name")
  if (missing(response)) stop("Need the response column name")
  if (missing(fam_link)) stop("Need the family(link) for the glmer")
  if (size %% 1 != 0) stop("Size needs to be a positive integer")
  if (size %% 1 == 0 && size < 0) stop("Size needs to be a positive integer")
  if (!is.null(first) && first %% 1 != 0) stop("First needs to be a positive integer")
  if (!is.null(first) && first %% 1 == 0 && first < 0) stop("First needs to be a positive integer")
  print(time1<- Sys.time()) #start time
  if(paste(fam_link)[2]== "logit") { m_err<- (pi^2)/3 }
  if(paste(fam_link)[2]== "probit") { m_err<- 1 }
  if(paste(fam_link)[2]== "sqrt") { m_err<- 0.25 }
if (size == 1) {
  if (quasi == F)  { jack<- matrix(0,ncol=5,nrow=length(observ[,1])) } #dam,sire,inter,residual,tot
  if (quasi == T)  { jack<- matrix(0,ncol=6,nrow=length(observ[,1])) } #add dispersion
  if (is.null(first)) { iter<- nrow(observ) }
  if (!is.null(first)) { iter<- first }
   for (i in 1:iter) {
   print(paste0("Removing observation: ",i," of ",length(observ[,1])))
  if (quasi == F)  {
  m<- glmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,")",sep="")),
    family=fam_link,data=observ[-i,])  }
  if (quasi == T)  {
  observ$dispersion<- as.factor(1:length(observ[,1]))
  m<- glmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,") + (1|dispersion)",sep="")),
    family=fam_link,data=observ[-i,])  }
  if(paste(fam_link)[2]== "log") { m_err<- log(1/exp(fixef(m)[1]) + 1) }
  tot<- sum(colSums(diag(VarCorr(m))),m_err)
  jack[i,]<- c(colSums(diag(VarCorr(m))),m_err,tot)
  col_names<- c(as.data.frame(VarCorr(m))$grp,"Residual","Total")  } #end loop
} #end single Jack
if (size > 1) {
 observS1 <- observ[sample(nrow(observ)),];rm(observ) #shuffle once
 observS2 <- observS1[sample(nrow(observS1)),];rm(observS1) #shuffle twice
 n_mod<- floor(nrow(observS2)/size)  #number of models in data frame
 if (quasi == F)  { jack<- matrix(0,ncol=5,nrow=n_mod) } #dam,sire,inter,residual,tot
 if (quasi == T)  { jack<- matrix(0,ncol=6,nrow=n_mod) } #add dispersion
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
  if (quasi == F)  {
  m<- glmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,")",sep="")),
    family=fam_link,data=observS2[-c(beg[i]:finn[i]),])  }
  if (quasi == T)  {
  observS2$dispersion<- as.factor(1:length(observS2[,1]))
  m<- glmer(formula= noquote(paste(response,"~ (1|",dam,") + (1|",sire,") + (1|",dam,":",sire,") + (1|dispersion)",sep="")),
    family=fam_link,data=observS2[-c(beg[i]:finn[i]),])  }
  if(paste(fam_link)[2]== "log") { m_err<- log(1/exp(fixef(m)[1]) + 1) }
  tot<- sum(colSums(diag(VarCorr(m))),m_err)
  jack[i,]<- c(colSums(diag(VarCorr(m))),m_err,tot)
  col_names<- c(as.data.frame(VarCorr(m))$grp,"Residual","Total")  } #end loop
}  #end block Jack
 jack<- as.data.frame(jack); colnames(jack)<- col_names
  colnames(jack)[which(colnames(jack)==dam)]<- "dam"
  colnames(jack)[which(colnames(jack)==sire)]<- "sire"
  colnames(jack)[which(colnames(jack)==noquote(paste(dam,":",sire,sep="")))]<- "dam:sire"
  jack$additive<- 4*jack$sire
  jack$nonadd<- 4*jack$'dam:sire'
  jack$maternal<- jack$dam- jack$sire
 print(Sys.time()- time1) #end time
 invisible(jack)  #after time
}
