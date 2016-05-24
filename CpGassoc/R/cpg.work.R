cpg.work <-
function(beta.values,indep,covariates=NULL,data=NULL,logit.transform=FALSE,
          chip.id=NULL,subset=NULL,random=FALSE,fdr.cutoff=.05,callarge=FALSE,fdr.method="BH",logitperm=FALSE,big.split=FALSE) {
gc()
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

if(callarge) {gc()}
nameholder<-list(deparse(substitute(beta.values)),deparse(substitute(chip.id)),
                  cpg.everything(deparse(substitute(indep))))
beta.values<-t(beta.values)
if(nrow(beta.values)==1) {beta.values<-t(beta.values)}

if(is.character(indep)) {warnings("\nindep is a character class, converting to factor\n")
        indep<-as.factor(indep)
        }

cpg.length(indep,nrow(beta.values),covariates, chip.id)
if(!is.null(data)){
  nameofdata<-deparse(substitute(data))
  thecheck<-nameofdata %in% search()
  if(!thecheck) stop("\nMust attach data before using data option in CpGassoc package",
                     "\nPlease attach and resubmit command\n")  
}

if(is.matrix(covariates)| length(covariates)==length(indep))  {
  if(is.character(covariates) &  is.matrix(covariates)) {
    stop("\nCan not analyze data with covariates given.\nNo characters allowed within",
          " a matrix")
          }
  covariates<-data.frame(covariates)
      }

if(!is.matrix(beta.values)) {beta.values<-as.matrix(beta.values)}
theholder<-list(covariates,chip.id)

levin<-is.factor(indep)
Problems<-which(beta.values<0 |beta.values >1)

if(is.character(beta.values[1,1])) {
  stop("beta.values contains characters, matrix/data.frame must be all numeric\n")
  }

beta.col<-ncol(beta.values)
fdr <- beta.col >= 100

if(fdr.method=="qvalue" & !fdr) {
    fdr.method="BH"
    warning("\nCan not perform qvalue method with less than a 100 CpG sites\n")
    }

if(callarge & big.split) {fdr=FALSE}

cpg.everything(complex(1),first=TRUE,logit.transform,Problems,beta.values,logitperm)

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
          }

if(!is.null(subset)){
  beta.values<-as.matrix(beta.values[subset,])
  indep<-indep[subset]
  if(!is.null(covariates)){
      covariates<-data.frame(covariates[subset,]) }
  if(!is.null(chip.id)) {
       chip.id<-chip.id[subset]
        }
         }

designmatrix<-design(covariates,indep,chip.id,random)
f.design<-designmatrix[[1]]
r.design<-designmatrix[[2]]

n=nrow(f.design) 
beta.values<-beta.values[designmatrix[[3]],]
if(!is.null(chip.id)){
	chip.id<-chip.id[designmatrix[[3]]]
	}
if(callarge) {rm(designmatrix)
  gc()
}
test.stat <- matrix(NA,beta.col,3)

if (!levin) {
  e.s<-matrix(ncol=2,nrow=beta.col)
  stderror<-matrix(NA,beta.col)                  
  }
              
else{ 
    e.s<-stderror<-NULL
   }

df.gc<-matrix(nrow=beta.col,ncol=2)

df.gc[,1]<-ifelse(levin,nlevels(indep)-1,1)

if (random) {
  if(is.null(chip.id)) {
    stop("No chip.id was input with the random call\n")
      }
if(!levin | is.null(covariates)) { 
    i.p<-NULL
    }
if(!is.null(covariates) | !levin) {
    p<-f.design[,-c(1:nlevels(indep))]
     }
if(levin){
 if(!is.null(covariates)) {i.p<-model.matrix(~indep)[,-1]}
 else { p<-model.matrix(~indep)[,-1]}
      }
  ran.function<-cpg.everything(p,i.p,levin,chip.id)
  ran.info<-t(apply(beta.values,2,ran.function))

  test.stat[,1:2]<-ran.info[,1:2]
  df.gc[,2]<-ran.info[,3]
  if(!levin) {
    e.s<-ran.info[,4:5]
    stderror[,1]<-abs(e.s[,2]/ran.info[,1])
     
    }
    
  if (sum(is.na(ran.info[,2]))>0){
    cpg.everything(complex(1),first=FALSE,logit.transform,sum(is.na(ran.info[,2])))
              }
             
  rm(ran.info)
  gc()}
    
else {
  if(!is.matrix(beta.values)) {beta.values<-as.matrix(beta.values)}
  nonmissing<-which(apply(is.na(beta.values),2,sum)==0) 
  non.m.beta<-as.matrix(beta.values[,nonmissing])
  if(callarge) gc()

  beta<-tryCatch(solve(t(f.design) %*% f.design)%*% t(f.design) %*% non.m.beta,
                 error = function(e) NULL)
  dfl<-list(df=n-ncol(f.design),df0=n-ncol(r.design))
  df.gc[nonmissing,2]<-dfl$df

  if(!is.null(beta)) {
  
    ressq<-(non.m.beta-f.design %*% beta)^2
    ssef<-t(ressq)%*% matrix(1,dim(ressq)[1],1)
    if(callarge) {
        rm(ressq)
        gc()}
    beta0<-solve(t(r.design) %*% r.design)%*% t(r.design) %*% non.m.beta
    r.ressq<-(non.m.beta-r.design %*% beta0)^2
    sser<-t(r.ressq)%*% matrix(1,dim(r.ressq)[1],1)   
    test.stat[nonmissing,1]<-(sser-ssef)/ssef*(dfl$df/(dfl$df0-dfl$df))

    test.stat[nonmissing,2]<-pf(test.stat[nonmissing,1],df1=(dfl$df0-dfl$df),df2=dfl$df,lower.tail=FALSE)

    #This rounding is for when cpgwork gets called during a permutation run
   
  test.stat[which(test.stat[,1]<0 & test.stat[,1] > -.001),1]<-0
     }
  if(is.null(beta)) {

     if(!levin) {
        mono.lm2<-tryCatch(lm(non.m.beta~f.design[,-1]), error = function(e) NULL)
        mono.results<-cpgassocsummary(mono.lm2)
        test.stat[nonmissing,1:2]<-mono.results[,c(1,3)]
        stderror[nonmissing,]<-mono.results[,2]
        if(is.null(ncol(mono.lm2$results))) {
               coeffic<-as.matrix(coef(mono.lm2))
                }
        else {coeffic<-coef(mono.lm2)}
        e.s[nonmissing,]<-cbind(t(as.matrix(colMeans(r.design)))%*% coeffic[-2,],
                      t(coeffic[2,]))
                      }
    else {
       if(ncol(r.design)==1) {
          mono.results<-aov(non.m.beta~f.design[,2:(nlevels(indep))])
                      }
       if(ncol(r.design)>1) {
          mono.results<-aov(non.m.beta~r.design[,-1]+f.design[,2:(nlevels(indep))])
                   }
       test.stat[nonmissing,1:2]<-cpgassocsummary(mono.results)
               }                                            
            }
  if(callarge) {
      rm(non.m.beta,sser,ssef,beta0,r.ressq)
      gc()
      }
  if(length(nonmissing)!=beta.col) {
    probsites<-which(is.na(test.stat[,2]))
    if(levin) {
      i.p<-f.design[,2:(nlevels(indep))]
      }


      
     for(i in probsites) {
        miss.check<-table(indep[!is.na(beta.values[,i])])
        miss.check<-ifelse(sum(miss.check!=0)>1,1,0)
        if(!levin) {
          temp<-tryCatch(summary(lm(beta.values[,i]~f.design[,-1],x=TRUE)), error = function(e) NULL)
          
          if(is.null(temp) | miss.check==0)    {
            stderror[i,]<-test.stat[i,1:2]<-df.gc[i,1:2]<-e.s[1,]<-NA

            }
          else{
       
            df.gc[i,2]<-c(temp$df[2])
            stderror[i,]<-coef(temp)[2,2]
            test.stat[i,1:2]<-coef(temp)[2,3:4]
            
            if(ncol(r.design)==1) {
              e.s[i,]<-coef(temp)[1:2,1]
                }
            else{
      
              altdesign<-f.design[!is.na(beta.values[,i]),-2]
              thecolmeans<-colMeans(altdesign[,colSums(altdesign)!=0])
              if( length(thecolmeans)!=length(coef(temp)[-2,1]))  {
                almost<-gsub(", -1","",gsub("f.design","",names(coef(temp)[-2,1])))
                almost<-gsub("[[]]","",almost)
                thecolmeans<-thecolmeans[which(names(f.design)[,-2] %in% almost)]
                }
              newint<-sum(thecolmeans*(coef(temp)[-2,1]))
              e.s[i,]<-c(newint,coef(temp)[2,1])
            }
            }
            }
         else{ 
   
           temp <- tryCatch(as.matrix(anova(lm(beta.values[,i]~r.design+i.p))), error = function(e) NULL)
        if(is.null(temp) | miss.check==0) {
            test.stat[i,1:2]<-NA
            }
          else{
            df.gc[i,2]<-tail(temp[,1],1)
            test.stat[i,1:2] <- temp[nrow(temp)-1,4:5]
          }
          } 
         }

        }
 
  if(!is.null(beta)) {

  if(!levin) {
    test.stat[nonmissing,1]<-sqrt(test.stat[nonmissing,1])*sign(beta[2,])
    if(beta.col>1) { 
    e.s[nonmissing,]<-cbind(colSums(colMeans(r.design)*matrix(beta[-2,],ncol(r.design))),beta[2,])
       }
    else {
      e.s[nonmissing,]<-cbind(colSums(as.matrix(colMeans(r.design)*matrix(beta[-2,],ncol(r.design)))),beta[2,])
      }
    stderror[nonmissing,]<-abs(e.s[nonmissing,2]/test.stat[nonmissing,1])
       }}
   if(callarge) {
   rm(beta)
   gc()}
   if(sum(is.na(test.stat[,1]))>0) {
      warning(sum(is.na(test.stat[,1]))," sites were not able to be analysed.\n",
          "Test statistics for these were set to NA")
      }  
     
     }


  gcvalue<-df.gc[1,1]*median(ifelse(rep(levin,beta.col),test.stat[,1],test.stat[,1]**2),na.rm=TRUE)/qchisq(.5,df.gc[1,1])
  gcvalue<-ifelse(gcvalue<1,1,gcvalue)
  gc.p.val<-pf(ifelse(rep(levin,beta.col),test.stat[,1],test.stat[,1]**2)/gcvalue,df.gc[,1],df.gc[,2],lower.tail=FALSE)
  if(random & gcvalue==1) {
        gc.p.val<-test.stat[,2]
        }
  




if (fdr & fdr.method=="qvalue") {
  FDR<-matrix(NA,beta.col)
  if (!requireNamespace("qvalue", quietly = TRUE)) {
       stop("qvalue needed for this to work. Please install it.",
                      call. = FALSE)
            }
  holder<-which(!is.na(test.stat[,2]) | !is.nan(test.stat[,2]))
  qv<-tryCatch(qvalue::qvalue(test.stat[holder,2]), error = function(e) NULL)
  if(is.null(qv)) {
    qv <- tryCatch(qvalue::qvalue(test.stat[holder,2], pi0.method = "bootstrap"), error = function(e) NULL)
    if(is.null(qv)) {
      fdr<-FALSE
      warning("qvalue package failed to process p-values.\nUsing Benjamini & Hochberg \n")
      fdr.method="BH"
      }}

  if(fdr) {
    FDR[holder,1]<-qv$qvalue
    test.stat<-cbind(test.stat,FDR)
      }}

if(!fdr | fdr.method!="qvalue") {
        if(fdr.method=="qvalue"){
                 test.stat<-cbind(test.stat,p.adjust(test.stat[,2],"BH"))
              }
        if(fdr.method!="qvalue"){
                 test.stat<-cbind(test.stat,p.adjust(test.stat[,2],fdr.method))
              }
        }
if(beta.col==1) {callarge<-FALSE}

if(is.null(dimnames(beta.values))& !callarge & beta.col>1) {colnames(beta.values)<-paste("X",1:beta.col,sep="") }
if(is.null(dimnames(beta.values))& !callarge & beta.col==1) {dimnames(beta.values)<- list(seq(1,length(beta.values)),paste("X",1:beta.col,sep=""))}
         
if(is.null(colnames(beta.values)) & !callarge) { colnames(beta.values)<-paste("X",1:beta.col,sep="")}

if(beta.col==1){callarge<-TRUE}

test.stat<-data.frame(colnames(beta.values),test.stat,stringsAsFactors=FALSE)
holmadjust<-p.adjust(test.stat[,3],"holm")
test.stat[,4]<-ifelse(holmadjust>.05,FALSE,TRUE)
if(sum(is.na(test.stat[,3]))>0) {
test.stat[which(is.na(test.stat[,3])),4]<-FALSE
      }

nonfactorinfo<-data.frame(df.top=df.gc[,1],df.bottom=df.gc[,2],stringsAsFactors=FALSE)
if(!levin) {
    nonfactorinfo<-data.frame(nonfactorinfo,adj.intercept=e.s[,1],
                effect.size=e.s[,2],std.error=stderror,stringsAsFactors=FALSE)
          }
row.names(nonfactorinfo)<-test.stat[,1]

names(test.stat)<-cpg.everything(fdr,perm=FALSE,levin)
test.stat<-data.frame(test.stat,gc.p.value=gc.p.val,stringsAsFactors=FALSE)
Num.Cov<-ncol(data.frame(covariates))
if(!is.null(chip.id) & !random) {Num.Cov<-Num.Cov+1}          
INFO<-data.frame(Min.P.Observed=min(test.stat$P.value,na.rm = TRUE),Num.Cov,fdr.cutoff,FDR.method=fdr.method,Phenotype=nameholder[[3]],
                  betainfo=nameholder[[1]],chipinfo=nameholder[[2]],random,logittran=logit.transform,stringsAsFactors=FALSE)

holm.sites<-subset(test.stat,test.stat[,4]==TRUE)
FDR.sites<- subset(test.stat,FDR<=fdr.cutoff)

info.data<-list(results=test.stat,Holm.sig=holm.sites,FDR.sig=FDR.sites,
        info=INFO,indep=indep,covariates=theholder[[1]],chip=theholder[[2]],coefficients=nonfactorinfo)
rm(test.stat,holm.sites,INFO,FDR.sites,f.design,r.design,Problems,beta.values,
            e.s,indep,theholder)
gc()
class(info.data)<-"cpg"
info.data

  }
