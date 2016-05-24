design.NULL <-
function(covariates,indep,chip.id,random) {
 whomiss<-data.frame(indep)
   allzero<- function(x) nlevels(as.factor(x))==1
 rowmiss<-function(x) {sum(is.na(x))}
 if (!random & !is.null(chip.id))  {
                 full<-model.matrix(~indep + factor(chip.id))
                 if(is.factor(indep)) {
                  reduced<-full[,-(2:(nlevels(indep)))]
                    }
                 else { reduced<-full[,-2]}
                 whomiss<-data.frame(whomiss,chip.id)  
                   }
 
 if(random | is.null(chip.id)) {
                   full<-model.matrix(~indep)
                   reduced<-matrix(1,nrow(full))
                 }

 
 problem<-which(apply(as.matrix(full[,2:ncol(full)]),2,allzero))
 if(length(problem)>0) {
     full<-full[,-(problem+1)]
     reduced<-reduced[,-which(apply(reduced,2,allzero))]
     }
 missingval<-t(apply(whomiss,1,rowmiss))
 gu<-vector(length=length(missingval))
 gu[which(missingval==0)]=TRUE
 gu[which(missingval!=0)]=FALSE
  list(full=full,reduced=reduced,miss=gu) }
