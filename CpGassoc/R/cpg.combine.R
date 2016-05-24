cpg.combine <-
function(allvalues,fdr.method="BH",fdr.cutoff=.05)    {
  if(length(allvalues)==1) {
    return(allvalues[[1]])
      }
  else {
  correctval<-list()
  right.val<-which(sapply(allvalues,is,"cpg"))
  for(i in right.val) {
        j=length(correctval)+1
        correctval[[j]]<-allvalues[[i]]
    }
  if(length(correctval)==1) {
    return(correctval)
    }
  levin<-is.factor(correctval[[1]]$indep)
  



  tnamepval<-function(x) {x$results[,1:3]}
  betainf<-function(x) {x$info$betainfo}
  temp<-sapply(correctval,tnamepval,USE.NAMES = FALSE)

  test.stat<-data.frame(unlist(temp[1,]),unlist(temp[2,]),unlist(temp[3,]))
  betainfo<-unlist(lapply(correctval,betainf))
  temp2<-t(sapply(correctval,coef))
  nonfactorinfo<-data.frame(df.top=unlist(temp2[,1]),df.bottom=unlist(temp2[,2]))
  if(!levin){
      nonfactorinfo<-data.frame(nonfactorinfo,cbind(unlist(temp2[,3]),unlist(temp2[,4]),unlist(temp2[,5])))

      names(nonfactorinfo)[3:5]<-c("adj.intercept","effect.size","std.error")
     }
  row.names(nonfactorinfo)<-test.stat[,1]
  
    

beta.col<-nrow(test.stat)
  
gcvalue<-median(nonfactorinfo[1,1]*ifelse(rep(levin,beta.col),test.stat[,2],test.stat[,2]**2),na.rm=TRUE)/qchisq(.5,nonfactorinfo[1,1])
gcvalue<-ifelse(gcvalue<1,1,gcvalue)
gc.p.val<-pf(ifelse(rep(levin,beta.col),test.stat[,2],test.stat[,2]**2)/gcvalue,
              nonfactorinfo[,1],nonfactorinfo[,2],lower.tail=FALSE)

if(correctval[[1]]$info$random & gcvalue<1) {
        gc.p.val<-test.stat[,3]
       }
  

betainfo<-gsub("[[:digit:]]","",betainfo)
betainfo<-gsub("[:,:]:]","",betainfo)
betainfo<-gsub("[[]","",betainfo)
betainfo<-gsub("[]:]","",betainfo)
betainfo<-gsub(",","",betainfo)
betainfo<-gsub("[[:blank:]]","",betainfo)
betainfo<-as.character(unique(betainfo))
if(length(betainfo)!=1) {betainfo="Could not get original source of beta values"}
test.stat[,4]<-NA
holmadjust<-p.adjust(test.stat[,3],"holm")
test.stat[,4]<-ifelse(holmadjust>.05,FALSE,TRUE)
names(test.stat)[4]<-"Holm.sig"

fdr<-TRUE
FDR<-matrix(NA,nrow(test.stat))
holder<-which(!is.na(test.stat[,3]) | !is.nan(test.stat[,3]))
if(fdr.method=="qvalue")  {
       if (!requireNamespace("qvalue", quietly = TRUE)) {
            stop("qvalue needed for this to work. Please install it.",
                      call. = FALSE)
            }
qv<-tryCatch(qvalue::qvalue(test.stat[holder,3]), error = function(e) NULL)
if(is.null(qv)) {
    qv <- tryCatch(qvalue::qvalue(test.stat[holder,3], pi0.method = "bootstrap"), error = function(e) NULL)
    if(is.null(qv)) {
      fdr<-FALSE
      warning("qvalue package failed to process p-values.\nUsing Benjamini & Hochberg\n")
      fdr.method="BH"
      }}

if(fdr) {
      FDR[holder,1]<-qv$qvalue
      test.stat<-cbind(test.stat,FDR)
      }
      }
if(!fdr | fdr.method!="qvalue") {
   FDR<-p.adjust(test.stat[,3],fdr.method)
   test.stat<-cbind(test.stat,FDR)
    }
names(test.stat)<-cpg.everything(fdr,perm=FALSE,levin)
test.stat<-data.frame(test.stat,gc.p.value=gc.p.val,stringsAsFactors=FALSE)
fdr.sites<-test.stat[which(test.stat$FDR<fdr.cutoff),]

holm.sites<-test.stat[which(holmadjust<.05),]
INFO<-data.frame(Min.P.Observed=min(test.stat$P.value,na.rm = TRUE),Num.Cov=correctval[[1]]$info$Num.Cov,fdr.cutoff,FDR.method=fdr.method,
              Phenotype=correctval[[1]]$info$Phenotype,betainfo,chipinfo=correctval[[1]]$info$chipinfo,random=correctval[[1]]$info$random,logittran=correctval[[1]]$info$logittran
              , stringsAsFactors=FALSE)
info.data<-list(results=test.stat,Holm.sig=holm.sites,FDR.sig=fdr.sites,
                  info=INFO,indep=correctval[[1]]$indep,covariates=correctval[[1]]$covariates,
                  chip=correctval[[1]]$chip,coefficients=nonfactorinfo)


rm(allvalues,test.stat,correctval)
gc()
class(info.data)<-"cpg"
info.data
      }

}
