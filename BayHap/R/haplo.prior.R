haplo.prior<-
function(res.freq=NULL,haplos.name=NULL,coeff=NULL,sd.coeff=NULL,OR=NULL,CI.OR=NULL,sign.OR=NULL){
     
     
     if(length(haplos.name)==0){
       haplo.mat<-NULL
     }else{
       prior.val<-1
       if (any(!(haplos.name%in%print.freq(res.freq)$Haplotypes))) stop("Incorrect haplotype names.\n")
       if (is.null(OR)){    
           if ((length(haplos.name)!=length(coeff))|(length(haplos.name)!=length(sd.coeff))| (length(coeff)!=length(sd.coeff))) stop("Different long variable found for haplos.name, coeff and sd.coeff \n")
       
       }else{
            if ((length(haplos.name)!=length(OR))|(length(haplos.name)!=nrow(CI.OR))| (length(OR)!=nrow(CI.OR))) stop("Different long variable found for haplos.name, OR and CI.OR \n")
            if (length(sign.OR)!=length(haplos.name)) stop("Different long found for haplos.name and sign.")
            if (is.null(sign.OR)) stop("The signification of the CI for each OR must be specified in the sign argument.")
            ci.coef<-log(CI.OR)
            i<-1
            sd<-NULL
            for(i in 1:length(OR)){
               sd<-c(sd,(ci.coef[i,1]-log(OR[i]))/-(qnorm(1-sign.OR/2)))
            }
            coeff<-log(OR)
            sd.coeff<-sd
       }
  
           haplo.info.prior<-data.frame(haplos.name,mean.haplo=coeff,sd.haplo=sd.coeff)
           names(haplo.info.prior)[1]<-"Haplotypes"
           names(haplo.info.prior)[2]<-"mean.haplo"
           names(haplo.info.prior)[3]<-"sd.haplo"
           haplo.mat<-merge(
print.freq(res.freq),haplo.info.prior,by="Haplotypes",all.x=TRUE,all.y=TRUE)
              
     }
haplo.mat

}

