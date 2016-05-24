cpg.perm <-
function(beta.values,indep,covariates=NULL,nperm,data=NULL,seed=NULL,
     logit.transform=FALSE,chip.id=NULL,subset=NULL,random=FALSE,fdr.cutoff=.05,fdr.method="BH",large.data=TRUE) {

name.holder<-list(deparse(substitute(beta.values)),deparse(substitute(chip.id)),cpg.everything(deparse(substitute(indep))))

if(is.null(ncol(beta.values))) {beta.values<-as.matrix(beta.values)}
beta.row<-nrow(beta.values)
beta.col<-ncol(beta.values)
if(class(covariates)=="formula") {
  variables<-gsub("[[:blank:]]","",strsplit(as.character(covariates)[2],"+",fixed=TRUE)[[1]])
  covariates<-data.frame(eval(parse(text=variables[1])))
  names(covariates)=variables[1]
  if(length(variables)>1) {
    for(i in 2:length(variables)) {
       covariates<-cbind(covariates,eval(parse(text=variables[i])))
       names(covariates)=variables[1:i]
      }
     } 
      }

cpg.length(indep,beta.col,covariates, chip.id)
if(is.character(indep)) {warnings("\nindep is a character class, converting to factor\n")
        indep<-as.factor(indep)
        }



if(!is.null(data)){
  nameofdata<-deparse(substitute(data))
  thecheck<-nameofdata %in% search()
  if(!thecheck) stop("\nMust attach data before using data option in CpGassoc package",
                     "\nPlease attach and resubmit command\n")  
}
if(is.matrix(covariates)| length(covariates)==length(indep) ) {
  if(is.character(covariates) &  is.matrix(covariates)) {
    stop("\nCan not analyze data with covariates given.\nNo characters allowed within",
          " a matrix")
          }
  else {
      covariates<-data.frame(covariates)
      }}
levin<-is.factor(indep)
Problems<-which(beta.values<0 |beta.values >1)

ob.data<-cpg.assoc(beta.values,indep,covariates,data,logit.transform,chip.id,subset,random,fdr.cutoff,fdr.method=fdr.method,large.data=large.data)
ob.data$info$Phenotype<-name.holder[[3]]
Min.P.Observed<-ob.data$info[1,1]


fdr <- beta.row>= 100
if(fdr.method=="qvalue" & !fdr) {
  fdr.method="BH"
  warning("\nCan not perform qvalue method with less than a 100 CpG sites\n")
  }
if(nperm>=100) {
  perm.pval<-matrix(NA,beta.row,nperm)
  if(!levin) {perm.tstat<-matrix(NA,beta.row,nperm)}
 }

cpg.everything(complex(1),first=TRUE,logit.transform,Problems,beta.values)
if(logit.transform) {
  beta.values<-as.matrix(beta.values)
  if (length(Problems)!=0) {
    beta.values[Problems]<-NA
            }
    onevalues<-which(beta.values==1)
    zerovalues<-which(beta.values==0)
    if(length(onevalues)>0 | length(zerovalues)>0) {
      if(length(onevalues)>0) {
         beta.values[onevalues]<-NA
         beta.values[onevalues]<-max(beta.values,na.rm=T)
            }
      if(length(zerovalues)>0) {
        beta.values[zerovalues]<-NA
        beta.values[zerovalues]<-min(beta.values,na.rm=T)
          }
          }
  beta.values=log(beta.values/(1-beta.values))
  beta.values<-data.frame(beta.values)
          }
     
if(!is.null(subset)){
  if(!is.null(chip.id)) {
    chip.id<-chip.id[subset]
        }
  beta.values<-beta.values[,subset]
  indep<-indep[subset]
  if(!is.null(covariates)) {
      covariates<-covariates[subset,]
      }
  subset=NULL
          }





Permutation <- matrix(nrow = nperm,ncol = 3)
compleval<-design(covariates,indep,chip.id,random)[[3]]
gc.Permutation<-matrix(nrow=nperm,ncol=3)
indep<-indep[compleval]
covariates<-data.frame(covariates[compleval,])
if(ncol(covariates)==0 & nrow(covariates)==0) {covariates=NULL}
beta.values<-beta.values[,compleval]
if(is.null(dim(beta.values))) {beta.values<-as.matrix(beta.values)}
chip.id<-chip.id[compleval]

for(i in 1:nperm) {
  if (!is.null(seed)) {
    set.seed(i*22424-seed)
        }
  Perm.var <- sample(indep);
  
  answers<-cpg.assoc(beta.values,Perm.var,covariates,data,logit.transform=FALSE
                ,chip.id,subset,random,fdr.cutoff,fdr.method=fdr.method,logitperm=TRUE,large.data=large.data)
    if(random) {
      problems<-sum(is.na(answers$results$P.value))
      if (problems>0){
        cpg.everything(complex(1),first=FALSE,logit.transform,problems)
                   }
               }
  if(nperm>=100 ) {
    perm.pval[,i]<-answers$results$P.value
    if(!levin) {perm.tstat[,i]<-answers$results$T.statistic}
      }
  Permutation[i,1:3] <- c(answers$info$Min.P.Observed,nrow(answers$Holm.sig),nrow(answers$FDR.sig))
  if (fdr.method=="qvalue") {
      fdr.adj<-tryCatch(qvalue::qvalue(answers$results$gc.p.value), error = function(e) NULL)
        if(is.null(fdr.adj)) {
          fdr.adj <- tryCatch(qvalue::qvalue(answers$results$gc.p.value, pi0.method = "bootstrap"), 
                              error = function(e) NULL)
          if(is.null(fdr.adj)) {
              fdr.method="BH"
        }}}
  if(fdr.method!="qvalue") {
          fdr.adj<-p.adjust(answers$results$gc.p.value,fdr.method)
          }
  if(fdr.method=="qvalue"){
          fdr.adj<-fdr.adj$qvalue
  }
  gc.Permutation[i,1:3]<-c(min(answers$results$gc.p.value,na.rm=TRUE),
                            sum(p.adjust(answers$results$gc.p.value,"holm")<.05),
                            sum(fdr.adj<.05))
  rm(answers)
  gc()
           }
      
Permutation<-data.frame(Permutation)

p.value.p <- sum(Permutation[,1]<=Min.P.Observed)/nperm
p.value.holm <- sum(Permutation[,2]>=nrow(ob.data$Holm.sig))/nperm
p.value.FDR <- sum(Permutation[,3]>=nrow(ob.data$FDR.sig))/nperm 

if(is.null(seed)) {seed<-"NULL"}
p.value.matrix <- data.frame(p.value.p,p.value.holm,p.value.FDR,nperm,seed)
names(Permutation)<-cpg.everything(fdr,perm=TRUE)
colnames(gc.Permutation)<-names(Permutation)
perm.data<-list(permutation.matrix=Permutation,perm.p.values=p.value.matrix,gc.permutation.matrix=gc.Permutation)
perm.data<-append(perm.data,ob.data)
names(perm.data)<-c("permutation.matrix","perm.p.values","gc.permutation.matrix",names(ob.data))
if(nperm>=100) {
  perm.pval<-apply(perm.pval,2,sort)
  perm.pval<- -log(perm.pval,base=10)
  perm.data$perm.pval<-t(apply(perm.pval,1,quantile,probs=c(.025,.975)))
  if(!levin ) {
  perm.tstat<-apply(perm.tstat,2,sort)
  perm.data$perm.tstat<-t(apply(perm.tstat,1,quantile,probs=c(.025,.975)))
    }}
perm.data$info$betainfo<-name.holder[[1]]
rm(Permutation,ob.data,p.value.matrix)
gc()

class(perm.data)<-"cpg.perm"
perm.data
     }
