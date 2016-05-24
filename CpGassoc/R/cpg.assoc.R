cpg.assoc <-
function(beta.val,indep,covariates=NULL,data=NULL,logit.transform=FALSE,
              chip.id=NULL,subset=NULL,random=FALSE,fdr.cutoff=.05,large.data=TRUE,fdr.method="BH",logitperm=FALSE)  {
p.name.holder<-list(deparse(substitute(beta.val)),deparse(substitute(chip.id)),cpg.everything(deparse(substitute(indep))))
beta.val<-as.matrix(beta.val)
gc()

if(is.null(ncol(beta.val))) {beta.val<-as.matrix(beta.val)}
if(ncol(beta.val)==1) beta.val<-t(beta.val)

beta.row<-nrow(beta.val)
beta.col<-ncol(beta.val)

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
cpg.length(indep,beta.col,covariates,chip.id)
if(is.character(indep)) {warnings("\nindep was character class, converted to factor\n")
        indep<-as.factor(indep)
        }

gc()


if(!large.data) {
    results<-cpg.work(beta.val,indep,covariates,data,logit.transform,
                      chip.id,subset,random,fdr.cutoff,callarge=FALSE,fdr.method,logitperm)
    }
              
           
else  {

if(is.null(row.names(beta.val))) {
      row.names(beta.val)<-paste("X",1:beta.row,sep="")
      }
    if(.Platform$OS.type=="windows") {
          maxmemo<-memory.limit()/15
    
      }
else {
    whatsystem<-system('uname',intern=TRUE)
    if(whatsystem %in% c("Linux","Solaris")){
      info<-system('free',intern=TRUE)
      maxmemo<-as.numeric(strsplit(info[2],"\\s+")[[1]][4])/(1024*15)
            }
     else{
      top_output<-system('top -l1 -n 20 | grep -Ei "mem|vm"',intern=TRUE)
      there<-gsub("[[:alpha:]]","",top_output[2])
      mborgb<-gsub("[[:digit:]]","",top_output[2])
      mborgb<-gsub("free","", mborgb)
      locat<-gregexpr("[[:alpha:]]",mborgb)[[1]]
      mbgb<-substr(mborgb,locat[length(locat)],locat[length(locat)])
      positi<-gregexpr("[[:blank:]]",there)[[1]]
      maxmemo<-as.numeric(substr(there,(positi[length(positi)-1]+1),(positi[length(positi)]-1)))/15
      if(toupper(mbgb)=="G") {maxmemo=maxmemo*1024}
        } }
    mainobjectsize<-object.size(beta.val)/1048576
    if(!is.matrix(beta.val)) {beta.val<-as.matrix(beta.val)}
    allresults<-list()
    gc()
    i<-as.numeric(ceiling(mainobjectsize/maxmemo))
    bigobjectsize<-mainobjectsize
    bigobjectsize<-mainobjectsize/i 
    while(bigobjectsize >maxmemo) {
      i<-i+1
      bigobjectsize<-mainobjectsize/i
      }
  
    div<-trunc(beta.row/i)
    if(i>1) {big.split=TRUE}
    else{ big.split=FALSE}
    if(logit.transform) {
        onevalues<-which(beta.val==1)
        zerovalues<-which(beta.val==0)
      if(length(onevalues)>0 | length(zerovalues)>0) {
         
        if(length(onevalues)>0) {
           beta.val[onevalues]<-NA
           beta.val[onevalues]<-max(beta.val,na.rm=T)
            }
        if(length(zerovalues)>0) {
           beta.val[zerovalues]<-NA
           beta.val[zerovalues]<-min(beta.val,na.rm=T)
          }
        }
    
          }
    for(j in 1:i){
        gc()
        if(j<i) {
           allresults[[j]]<-cpg.work(beta.val[1:div,],indep,covariates,data,logit.transform,
                                    chip.id,subset,random,fdr.cutoff,callarge=TRUE,fdr.method,logitperm,big.split=big.split)
           gc()
           beta.val<-beta.val[(div+1):nrow(beta.val),]
           gc()
              }
        else {
           allresults[[j]]<-cpg.work(beta.val,indep,covariates,data,logit.transform,
                                    chip.id,subset,random,fdr.cutoff,callarge=TRUE,fdr.method,logitperm,big.split=big.split)
           gc()
              }
             }

  results<-cpg.combine(allresults,fdr.method)

  rm(allresults,fdr.method)
  gc()
      }
  rm(beta.val)
  gc()
  results$info$betainfo<-p.name.holder[[1]]
  results$info$Phenotype<-p.name.holder[[3]]
  results$info$chipinfo<- p.name.holder[[2]]
  results   
       
  }
