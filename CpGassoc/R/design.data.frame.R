design.data.frame <-
function(covariates,indep,chip.id,random) {
    rowmiss<-function(x) {sum(is.na(x))}
      allzero<- function(x) nlevels(as.factor(x))==1
    Var<-paste(names(covariates),collapse="+")
    whomiss<-data.frame(covariates,indep)
    if (!random & !is.null(chip.id)){
      Var<-paste(Var,"factor(chip.id)",sep="+")
      whomiss<-data.frame(whomiss,chip.id)
          }
    missingval<-t(apply(whomiss,1,rowmiss))
    gu<-vector(length=length(missingval))
    gu[which(missingval==0)]=TRUE
    gu[which(missingval!=0)]=FALSE
    
    Varin<-paste("indep",Var,sep="+")   
    Varin<-paste("~",Varin,sep="")
    full<-model.matrix(formula(Varin),data=covariates)
    if(is.factor(indep)) {
     reduced<-full[,-c(2:(nlevels(indep)))]
     }
    else {
      reduced<-full[,-2]
        }
    problem<-which(apply(full[,2:ncol(full)],2,allzero))
    if(length(problem)>0) {
     full<-full[,-(problem+1)]
     reduced<-reduced[,-(which(apply(reduced[,2:ncol(reduced)],2,allzero))+1)]
     }
    list(full=full,reduced=reduced,miss=gu)
          }
