
mbmdr <-
function(y,data,order,covar=NULL,exclude=NA,risk.threshold=0.1,output=NULL,adjust=c("none","covariates","main effects","both"),first.model=NULL,list.models=NULL,use.logistf=TRUE,printStep1=FALSE,...)
{
 called <- match.call()
 adjust <- match.arg(adjust)
 if(use.logistf && !(class(called$family)=="call" && called$family[[1]]=="binomial" && called$family$link=="logit")) use.logistf <- FALSE
 	
## if(nullDist==TRUE){
##   data(NullDist)
##   null.dist <- NullDist
## }
## if(nullDist==FALSE) null.dist <- NULL

 #Initialization
 data <- as.data.frame(data)
 for (j in 1:ncol(data)) data[,j] <- factor(data[,j], exclude=exclude)
 
 EXPO <- list()
 EXPO[["0"]] <- (data==0)
 EXPO[["1"]] <- (data==1)
 EXPO[["2"]] <- (data==2)
 
 #Functions
 NextModel <- function(models=NULL){
   if(is.null(models)){
     switch(class(list.models),
       NULL = { ifelse(is.null(first.model), models <- list(c(ncol(data):(ncol(data)-order+1)),1,0,0), models <- list(first.model,1,0,0)) },
       integer = { models <- list(list.models,1,1,list.models) },
       numeric = { models <- list(list.models,1,1,list.models) },
       matrix = {models <- list(list.models[1,],1,nrow(list.models),list.models) },
       character = {if(length(list.models)==1){
       	              aux <- as.matrix(read.table(list.models))
       	              models <- list(aux[1,],1,nrow(aux),aux)
       	            }
       	            else models <- list(list.models,1,1,list.models) 
       	           }
     )
   }
   else{
     if(models[[3]]>0){
       if(models[[2]]==models[[3]]) return(NULL)
       models[[2]] <- models[[2]] + 1
       models[[1]] <- models[[4]][models[[2]],]
     }
     else{
       n <- length(models[[1]])
       model <- models[[1]][n:1]
       for(i in 1:n){
         if (model[i]>i){
           k <- model[i]-1
           model[1:i] <- c((k-i+1):k)
           models[[1]] <- model[n:1]
           lastmodel <- FALSE
           break
         }
         lastmodel <- TRUE
       }
       if(lastmodel) models <- NULL
     }
   }
   return(models)
 }
 
 OffsetFit <- function(resp,datos,covar,adjust,...){
  switch(adjust,
    "none"={ return(NULL) },
    "covariates"={ offsetgl <- glm(resp~.,data=data.frame(resp,covar),...) },
    "main effects"={ offsetgl <- glm(resp~.,data=data.frame(resp,datos),...) },
    "both"={ offsetgl <- glm(resp~.,data=data.frame(resp,datos,covar),...) }
  )
  offsaux <- rep(NA,nrow(datos))
  names(offsaux) <- rownames(datos)
  offsaux[names(offsetgl$linear.predictors)] <- offsetgl$linear.predictors
  return(offsaux)
 }
 
## EmpiricalPval <- function(null.dist,n,W,maf){
##   mafs <- sapply(maf,function(x) null.dist$mafs[order(abs(null.dist$mafs-x))[1]])
##   distname <- paste(sort(mafs,decreasing=TRUE),collapse="_")
##   if(null.dist[[distname]]$prob[null.dist[[distname]]$prob[,1]==n,2]==0) return(0)
##   m <- sum(null.dist[[distname]]$prob[,1]>0 & null.dist[[distname]]$prob[,2]>=null.threshold)
##   refer <- null.dist[[distname]]$dist[null.dist[[distname]]$dist[,1]==n,2]
##   pval  <- null.dist[[distname]]$prob[null.dist[[distname]]$prob[,1]==n,2]*m*(length(refer[refer>W]) / length(refer))
##   return(max(pval,1))
## }
 
## MAF <- function(snp,exclude){
##   auxtab <- table(snp,exclude=exclude)
##   freq <- c("0"=0,"1"=0,"2"=0)
##   for(f in names(auxtab)) freq[f] <- auxtab[f]
##   fa <- (2*freq["0"]+freq["1"])/(2*sum(auxtab))
##   return(min(fa,1-fa))
## }
 
 Regresion <- function(resp,exposed,offsetfit,adjust,use.logistf,...){
 	if(use.logistf){
    tab <- table(resp,exposed,exclude=exclude)
    if(any(tab==0) || sum(dim(tab))<4){
      if(adjust=="none"){
      	if(class(resp)=="factor") resp <- as.numeric(resp) - 1
        auxgl <- logistf(resp~.,data=data.frame(resp,exposed),pl=FALSE)
        return(auxgl)
      }
      #else cat("\n model: ",model, " - CORRECTION NEEDED")
    }
  }
  auxgl <- glm(resp~exposed,offset=offsetfit,data=data.frame(resp,exposed),...)
  return(auxgl)
 }

 Exposed <- function(snps,genotip){
  k <- 1
  genotip <- as.character(genotip)
  e <- EXPO[[genotip[k]]][,snps[k]]
  while(k<length(snps)){
    k <- k+1
    e <- as.logical(e * EXPO[[genotip[k]]][,snps[k]])
  }
  return(e)
 }
 
 Exposition <- function(snps,expolist){
  switch(class(expolist),
    matrix = { es <- apply(expolist, 1, function(x){Exposed(snps,x)})
               ifelse(ncol(es)>1, e <- as.logical(rowSums(es)),  e <- es)},
    numeric = { e <- Exposed(snps,expolist)}
  )
  return(e)
 }

 GenotipReg <- function(x,...){
 e <- 1*Exposition(model,x[-1])
  if(sum(e==0,na.rm=TRUE)==0 | sum(e==1,na.rm=TRUE)==0) return(c(0,0,0,1))
 r <- Regresion(resp=y,exposed=e,offsetfit=offsetfit,adjust=adjust,use.logistf=use.logistf,...)
 if(is.na(r$coefficients["exposed"])) return(c(0,0,0,1))
 if (class(r)[1]=="logistf"){
   coefic <- r$coefficients[2]
   sd.coef <- diag(r$var)[2]^0.5
   return(c(coefic, sd.coef, coefic/sd.coef, r$prob[2]))
 }
 return(summary(r)$coef["exposed",])
 }
  
 FirstStep <- function(model,exclude,...){
  tab <- table(data[,model],exclude=exclude)
  dimens <- dim(tab)
  n <- length(dimens)
  part <- matrix(,prod(dimens),(n + 5))
  part[,1] <- 1:prod(dimens)
  aux1 <- c(dimens,1)
  aux2 <- c(1,dimens)
  for (i in 1:n){
    part[,i+1] <- rep(rep(as.numeric(dimnames(tab)[[i]]),rep(prod(aux2[1:i]),length(dimnames(tab)[[i]]))),prod(aux1[(i+1):(n+1)]))
  }
  reg <- apply(part,1,function(x,...) return(GenotipReg(x[1:(n+1)],...)),...)
  part[,(n+2):(n+5)] <- t(reg)
  return(part)
 }


#Main program
 
 models <- NextModel()
 model <- models[[1]]
 offsetfit <- OffsetFit(resp=y,datos=data[,model],covar=covar,adjust=adjust,...)
 result <- data.frame()
 
 ## if(!is.null(null.dist)) maf <- apply(data[,model],2,MAF,exclude=exclude)
 
 if(!is.null(output) && !file.exists(output)) write.table(t(c(paste("SNP",1:order,sep=""),"NH","WH","PH","NL","WL","PL","MIN.P")),file=output,append=FALSE,row.names=FALSE,col.names=FALSE,sep=";")

 repeat{
   part <- FirstStep(model,exclude,...)
   
   if(printStep1){
   	 printStep1Out <- "default"
   	 if(!is.null(called$family) && ((class(called$family)=="call" && called$family[[1]]=="binomial") || called$family=="binomial")) printStep1Out <- "binomial"
   	 switch(printStep1Out,
   	   "binomial" = { ans <- data.frame(part[,c(2:(order+1))],NA,NA,part[,c(2+order,5+order)],factor("0",levels=c("0","H","L")))
                      colnames(ans) <- c(colnames(data[,model]),"cases","controls","beta","p.value","category")
                      for(i in 1:nrow(ans)){
                        e <- 1*Exposition(model,as.numeric(ans[i,1:order]))
                        auxtab <- table(e,y,exclude=exclude)       
                        ans[i,"cases"] <- ifelse(class(aux<-try(auxtab["1","1"],TRUE))=="try-error",0,aux)
                        ans[i,"controls"] <- ifelse(class(aux<-try(auxtab["1","0"],TRUE))=="try-error",0,aux)
                        ans[i,"category"] <- ifelse(ans[i,"p.value"]<=risk.threshold, ifelse(ans[i,"beta"]>0,"H","L"),"0")
                      }
                      cat("\n\nModel: ", model, "\n")
                      print(ans,digits=4,row.names=FALSE)
                    },
        "default" = { ans <- data.frame(part[,c(2:(order+1))],part[,c(2+order,5+order)],factor("0",levels=c("0","H","L")))
                      colnames(ans) <- c(colnames(data[,model]),"beta","p.value","category")
                      for(i in 1:nrow(ans)) ans[i,"category"] <- ifelse(ans[i,"p.value"]<=risk.threshold, ifelse(ans[i,"beta"]>0,"H","L"),"0")
                      cat("\n\nModel: ", model, "\n")
                      print(ans,digits=4,row.names=FALSE)
                    })
   }

   
   #second step
   h.list <- part[(part[,(order+2)]>0 & part[,(order+5)]<=risk.threshold),(2:(1+order))]
   l.list <- part[(part[,(order+2)]<0 & part[,(order+5)]<=risk.threshold),(2:(1+order))]
   
   if(!is.na(l.list[1])){
     switch(class(l.list), matrix={NL <- nrow(l.list)}, numeric={NL <- 1})
     l.e <- 1*Exposition(model,l.list)
     l.r <- Regresion(resp=y,exposed=l.e,offsetfit=offsetfit,adjust=adjust,use.logistf=use.logistf,...)
     if (class(l.r)[1]=="logistf"){
     coefic <- l.r$coefficients[2]
     sd.coef <- diag(l.r$var)[2]^0.5
     l.regout <- c(coefic, sd.coef, coefic/sd.coef, l.r$prob[2])
   }  
     else{
       l.regout <- summary(l.r)$coef["exposed",]
     }
     WL <- l.regout[3]^2
     PL <- l.regout[4]
##     PL <- ifelse(!is.null(null.dist), EmpiricalPval(null.dist,NL,WL,maf), NA)
     betaL <- l.regout[1]
   }
   else{
     NL <- 0
     betaL <- NA
     WL <- NA
     PL <- NA
   }

   if(!is.na(h.list[1])){
     switch(class(h.list), matrix={NH <- nrow(h.list)}, numeric={NH <- 1})
     h.e <- 1*Exposition(model,h.list)
     h.r <- Regresion(resp=y,exposed=h.e,offsetfit=offsetfit,adjust=adjust,use.logistf=use.logistf,...)
     if (class(h.r)[1]=="logistf"){
     coefic <- h.r$coefficients[2]
     sd.coef <- diag(h.r$var)[2]^0.5
     h.regout <- c(coefic, sd.coef, coefic/sd.coef, h.r$prob[2])
   }
     else{
       h.regout <- summary(h.r)$coef["exposed",]
     }
     WH <- h.regout[3]^2
     PH <- h.regout[4]
##     PH <- ifelse(!is.null(null.dist), EmpiricalPval(null.dist,NH,WH,maf), NA)
     betaH <- h.regout[1]
   }
   else{
     NH <- 0
     betaH <- NA
     WH <- NA
     PH <- NA
   }
   
   if (NL+NH>0){
   	  MIN.P <- min(PH,PL,na.rm=TRUE)
##     aux <- c(PH,PL)
##     ADJ.P <- ifelse(is.null(null.dist), NA, 2*min(aux[!is.na(aux)]))
     if(is.null(output)){
       result <- rbind(result,data.frame(t(colnames(data[,model])),as.numeric(NH),betaH,WH,PH,as.numeric(NL),betaL,WL,PL,MIN.P,row.names=NULL))
##       result <- rbind(result,data.frame(t(colnames(data[,model])),as.numeric(NH),WH,PH,as.numeric(NL),WL,PL,ADJ.P,row.names=NULL))
     }
     else
       write.table(data.frame(t(colnames(data[,model])),as.numeric(NH),betaH,WH,PH,as.numeric(NL),betaL,WL,PL,MIN.P,row.names=NULL),file=output,append=TRUE,row.names=FALSE,col.names=FALSE,sep=";")
##       write.table(data.frame(t(colnames(data[,model])),as.numeric(NH),WH,PH,as.numeric(NL),WL,PL,ADJ.P,row.names=NULL),file=output,append=TRUE,row.names=FALSE,col.names=FALSE,sep=";")
   }

   #next model
   models <- NextModel(models)
   if(is.null(models)) break
   model <- models[[1]]
   if(adjust=="main effects" | adjust=="both") 
   	  offsetfit <- OffsetFit(resp=y,datos=data[,model],covar=covar,adjust=adjust,...)
 }
 
## if(is.null(output) & nrow(result)>0) colnames(result) <- c(paste("SNP",1:order,sep=""),"NH","WH","PH","NL","WL","PL","ADJ.P")
 if(is.null(output) & nrow(result)>0) colnames(result) <- c(paste("SNP",1:order,sep=""),"NH","betaH","WH","PH","NL","betaL","WL","PL","MIN.P")
 if(!is.null(output)) result <- paste(c("file: ",output),collapse="")
 object <- list()
 object$call <- called
 object$y <- y
 object$data <- data
## object$exposition <- EXPO
 object$covar <- covar
 object$result <- result
 class(object) <- "mbmdr" 
 return(object)
}

