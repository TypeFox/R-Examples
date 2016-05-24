#summary.output.x12 <- function(object,fullSummary=FALSE,spectra.detail=FALSE,almostout=FALSE,rsd.autocorr=NULL,q2=FALSE,likelihood.stat=FALSE,aape=FALSE,id.rsdseas=FALSE,slidingspans=FALSE,...){
#	summaryworkhorse(object$dg,fullSummary=fullSummary,spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,identify=identify) 
#}
summary.output.workhorse <- function(x,fullSummary=FALSE,spectra.detail=FALSE,almostout=FALSE,rsd.autocorr=NULL,quality.stat=FALSE,likelihood.stat=FALSE,aape=FALSE,id.rsdseas=FALSE,slidingspans=FALSE,history=FALSE,identify=FALSE){
  #cat("File: \"",x$file,"\"",sep="","\n")	
  if(length(nchar(unlist(x$outlier)))==0)
    x$outlier<-"-"	
  if(fullSummary){
    spectra.detail=TRUE
    almostout=TRUE
    rsd.autocorr=c("acf","pacf","acf2")
    quality.stat=TRUE
    likelihood.stat=TRUE
    aape=TRUE
    id.rsdseas=TRUE 
    slidingspans=TRUE
    history=TRUE
    identify=TRUE
  }
  summary.output<-data.frame()
  summary.output[dim(summary.output)[1]+1,1]<-"Frequency"
  summary.output[dim(summary.output)[1],2]<-x$frequency	
  if(length(x$span)>1){
    span <- str_trim(unlist(strsplit(x$span,":")))
    span.index <- which(span=="span")
    summary.output[dim(summary.output)[1]+1,1]<-"Span"
    summary.output[dim(summary.output)[1],2]<-span[span.index+1]		
    span.index <- which(span=="modelspan")
    if(length(span.index)>0)
      modelspan <- c("Model Span",span[span.index+1])
    span.index <- which(span=="outlierspan")
    if(length(span.index)>0)
      outlierspan <- c("Outlier Span",span[span.index+1])
  }else{
    summary.output[dim(summary.output)[1]+1,1]<-"Span"
    summary.output[dim(summary.output)[1],2]<-x$span		
  }	
  if(x$x11regress=="no"){			
#		colnames(summary.output)<-c("Diagnostic","Series")
    summary.output[dim(summary.output)[1]+1,1]<-"X11 Regression"
    summary.output[dim(summary.output)[1],2]<-"FALSE"		
    summary.output[dim(summary.output)[1]+1,1]<-"Model Definition"
    if(x$automdl!="-"){
      summary.output[dim(summary.output)[1],2]<-paste("ARIMA Model:",unlist(x$arimamdl),"(Automatic Model Choice)")			
    }else{
      summary.output[dim(summary.output)[1],2]<- paste("ARIMA Model:",unlist(x$arimamdl))		
    }
    if(exists("modelspan"))
      summary.output[dim(summary.output)[1]+1,] <- modelspan
    if(x$transform=="Automatic selection"){
      summary.output[dim(summary.output)[1]+1,1]<-"Transformation"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$transform),":",unlist(x$autotransform))
    }else{
      summary.output[dim(summary.output)[1]+1,1]<-"Transformation"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$transform))
    }	
    
    summary.output[dim(summary.output)[1]+1,1]<-"Regression Model"
    summary.output[dim(summary.output)[1],2] <- paste(unlist(x$regmdl),collapse=" ")
#		cat("\n\tOutlier Detection\n")
    if(x$ifout=="Outlier detection performed"){
      summary.output[dim(summary.output)[1]+1,1]<-"Outlier detection performed"
      summary.output[dim(summary.output)[1],2]<-paste("TRUE")
      if(exists("outlierspan"))
        summary.output[dim(summary.output)[1]+1,] <- outlierspan
      
#			cat("Critical |t| for outliers:\t\n")
      for(i in 1:length(names(x$crit))){
        summary.output[dim(summary.output)[1]+1,1]<-names(x$crit)[[i]]
        summary.output[dim(summary.output)[1],2]<-paste(x$crit[[i]],collapse=" ")
      }
      summary.output[dim(summary.output)[1]+1,1]<-"Total Number of Outliers"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$nout))
      summary.output[dim(summary.output)[1]+1,1]<-"Nr of Automatically Identified Outliers"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$nautoout))
      if(almostout){
        summary.output[dim(summary.output)[1]+1,1]<-"Nr of Almost Outliers"
        summary.output[dim(summary.output)[1],2]<-paste(unlist(x$nalmostout))
      }
    }	
    else{
      summary.output[dim(summary.output)[1]+1,1]<-"Outlier detection performed"
      summary.output[dim(summary.output)[1],2]<-paste("FALSE")
    }
    
    rest<-unlist(lapply(strsplit(as.character(x$regmdl),"+",fixed=TRUE),function(x)gsub("^\\s+|\\s+$", "",x)))
    rest<-names(x)[which(names(x)%in%rest)]
    if(almostout){
      liste <- c("outlier","userdefined","leapyear","td","autooutlier",rest,"almostoutlier")
    }else
      liste <- c("outlier","userdefined","leapyear","td",rest,"autooutlier")
    liste<-liste[which(liste%in%names(x))]
    empty <- which(unlist(lapply(1:length(x),function(y)any(x[[y]]=="-"))))
    res <- as.data.frame(do.call(rbind,lapply(which(!liste %in% names(x[empty])),function(j){
                  if(!any(grepl(names(x[liste[j]]),names(x[[liste[j]]])))){
                    names(x[[liste[j]]])<-paste(names(x[liste[j]]),"_",names(x[[liste[j]]]),sep="")	
                  }
                  do.call(rbind,lapply(1:length(x[[liste[j]]]),function(i){
                            c(names(x[[liste[j]]][i]),unlist(x[[liste[j]]][[i]]))}))})))
    if(all(dim(res))>0){
      res[,2:4] <- apply(res[,2:4],2,function(x)as.numeric(formatC(as.numeric(as.character(x)),digits=3,format="f")))
      colnames(res)[1]<-"variable"
      res2 <- cbind(paste(1:length(res[,1]),"variable, coef, stderr, tval"),apply(res,1,paste,collapse=", "))
      summary.output<-rbind(summary.output,res2)
      if(!is.null(x[["derived.coef"]])){
        summary.output[dim(summary.output)[1]+1,1]<-"* Derived parameter estimates"
        summary.output[dim(summary.output)[1],2]<-paste(x[["derived.coef"]])
      }
    }
    if(likelihood.stat){
#			cat("\n\tLikelihood Statistics\n")
      lstat<-matrix(c("AIC","AICC","BIC","HQ ","Log Likelihood",x$aic,x$aicc,x$bic,x$hq,x$loglikelihood),ncol=2)	  
      summary.output <- rbind(summary.output,lstat)
    }	  
    
    if(aape && length(x$aape)>1){
#			cat("\nAverage absolute percentage error\n")
      mode<-ifelse(x$aape$aape.mode=="outofsample","out of sample","within sample")
      summary.output[dim(summary.output)[1]+1,1]<-"AAPE mode"
      summary.output[dim(summary.output)[1],2]<-paste(mode)
      
      aape.mat<-matrix(c("AAPE Last year","AAPE Last-1 year","AAPE Last-2 year","AAPE Last 3 years",x$aape$aape.0,x$aape$aape.1,x$aape$aape.2,x$aape$aape.3),ncol=2)	  
      summary.output <- rbind(summary.output,aape.mat)
    }
#		cat("\n\tSeasonal Adjustment\n\n")
    summary.output[dim(summary.output)[1]+1,1]<-"Identifiable Seasonality"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$id.seas))
    if(id.rsdseas){
      summary.output[dim(summary.output)[1]+1,1]<-"Residual Seasonality"
      if(x$id.rsdseas=="none")
        summary.output[dim(summary.output)[1],2]<-"none"
      else
        summary.output[dim(summary.output)[1],2]<-"yes"
    }		
    summary.output[dim(summary.output)[1]+1,1]<-"Seasonal Peaks"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$peaks.seas))
    summary.output[dim(summary.output)[1]+1,1]<-"Trading Day Peaks"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$peaks.td))
    summary.output[dim(summary.output)[1]+1,1]<-"Q Statistic"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$q))
    if(quality.stat){
      summary.output[dim(summary.output)[1]+1,1]<-"Q2 Statistic"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$q2))
      summary.output[dim(summary.output)[1]+1,1]<-"M1"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m1))
      summary.output[dim(summary.output)[1]+1,1]<-"M2"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m2))
      summary.output[dim(summary.output)[1]+1,1]<-"M3"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m3))
      summary.output[dim(summary.output)[1]+1,1]<-"M4"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m4))
      summary.output[dim(summary.output)[1]+1,1]<-"M5"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m5))
      summary.output[dim(summary.output)[1]+1,1]<-"M6"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m6))
      summary.output[dim(summary.output)[1]+1,1]<-"M7"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m7))
      summary.output[dim(summary.output)[1]+1,1]<-"M8"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m8))
      summary.output[dim(summary.output)[1]+1,1]<-"M9"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m9))
      summary.output[dim(summary.output)[1]+1,1]<-"M10"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m10))
      summary.output[dim(summary.output)[1]+1,1]<-"M11"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m11))
      
    }
    summary.output[dim(summary.output)[1]+1,1]<-"Nr of M stats outside limits"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$nmfail))
    if(spectra.detail){
      new.names<-function(z){
        for(i in 1:length(x[[z]])){
          if(length(x[[z]][[i]])>1){
            x[[z]][[i]] <- do.call(paste,x[[z]][[i]])	}}
        y<-as.data.frame(x[[z]],row.names="")			
        colnames(y)<-gsub(paste(z,".",sep=""),replacement="",names(unlist(x[[z]])),fixed=TRUE)
        return(y)}
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of the original series"
      summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcori")),":",new.names("spcori"),collapse=",",sep="")
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of the regARIMA model residuals"
      if(any(x$spcrsd!="-"))				
        summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcrsd")),":",new.names("spcrsd"),collapse=",",sep="")
      else
        summary.output[dim(summary.output)[1],2]<-"-"			
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of differenced seasonally adjusted series"
      summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcsa")),":",new.names("spcsa"),collapse=",",sep="")
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of modified irregular series"
      summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcirr")),":",new.names("spcirr"),collapse=",",sep="")
    }
    #		cat("\n\tSeasonal and Trend Moving Averages\n\n")
    summary.output[dim(summary.output)[1]+1,1]<-"SA decomposition"
    summary.output[dim(summary.output)[1],2]<-paste(x$samode[[length(x$samode)]])
    
    if(x$seasonalma[[1]]=="M.S.R."){
      summary.output[dim(summary.output)[1]+1,1]<-"Seasonal moving average"
      summary.output[dim(summary.output)[1],2]<-paste(x$seasonalma[[length(x$seasonalma)]],"(Based on msr size)",sep=" ")
    }else{
      summary.output[dim(summary.output)[1]+1,1]<-"Seasonal moving average"
      summary.output[dim(summary.output)[1],2]<-paste(x$seasonalma[[length(x$seasonalma)]])
    }
    summary.output[dim(summary.output)[1]+1,1]<-"Henderson filter"
    summary.output[dim(summary.output)[1],2]<-paste(x$trendma[[length(x$trendma)]],"-term",sep="")
    if(!is.null(rsd.autocorr)){
      #rsd.autocorr=c("rsd.acf","rsd.pacf","rsd.acf2")
      if("acf"%in%rsd.autocorr && !is.null(x$rsd.acf)){
        sig<-rep("",dim(x$rsd.acf)[1])
        sig[which(x$rsd.acf$pval<0.05 & x$rsd.acf$df.q>0)]<-"*"
        rsd.acf<-cbind(round(x$rsd.acf,digits=3),sig)
        summary.output<-rbind(summary.output,cbind(rep("acf: lag, sample.acf, stderr.acf, Ljung-Box.q, df.q, pval",dim(rsd.acf)[1]),paste(apply(rsd.acf[,-dim(rsd.acf)[2]],1,paste,collapse=", "),rsd.acf[,dim(rsd.acf)[2]])))
      }
      if("pacf"%in%rsd.autocorr && !is.null(x$rsd.pacf)){
        summary.output<-rbind(summary.output,cbind(rep("pacf: lag, sample.pacf, stderr.pacf",dim(x$rsd.pacf)[1]),apply(round(x$rsd.pacf,digits=3),1,paste,collapse=", ")))	
      }
      if("acf2"%in%rsd.autocorr && !is.null(x$rsd.acf2)){
        summary.output<-rbind(summary.output,cbind(rep("acf2: lag, sample.acf2, stderr.acf2",dim(x$rsd.acf2)[1]),apply(round(x$rsd.acf2,digits=3),1,paste,collapse=", ")))	
      }
    }
    if(slidingspans){
      if(length(grep("ss.",names(x)))>0){
        summary.output[dim(summary.output)[1]+1,1]<-"Sliding spans analysis performed"
        summary.output[dim(summary.output)[1],2]<-"TRUE"
        if(all(x$ss.options!="-")){
          ss.options<-cbind(c("Nr of spans","Length of spans","First period in first span","First year in first span"),as.character(x$ss.options))
          summary.output<-rbind(summary.output,ss.options)	
        }
        if(any(x$ss.seasTests!="-")){
          seasTests <- x$ss.seasTests
          seasTests <- cbind(paste(rownames(seasTests),": M7, Identifiable Seasonality",sep=""), 
              as.character(apply(seasTests[which(colnames(seasTests)%in%c("m7","idseas"))],1,function(x)paste(x,collapse=", "))))
          summary.output <- rbind(summary.output,seasTests)					
        }	
        if(all(x$ss.S1!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"S1 Table generated (Period means of seasonal factors)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          ss.S1 <- x$ss.S1[[1]]
          ss.S1 <- cbind(paste("S1: ",paste(colnames(ss.S1),collapse=", "),sep=""),
              apply(ss.S1,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,ss.S1)
          
          ss.S1 <- x$ss.S1[[2]]
          ss.S1 <- cbind(paste("S1.summary",paste(row.names(ss.S1),paste(colnames(ss.S1),collapse=", "),sep=": "),sep="."),
              as.character(apply(ss.S1,1,function(x)paste(x,collapse=", "))))
          summary.output <- rbind(summary.output,ss.S1)					
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"S1 Table generated (Period means of seasonal factors)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        if(all(x$ss.S2!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"S2 Table generated (Percentage of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          ss.S2 <- x$ss.S2
          rn.S2 <- which(c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")%in%row.names(ss.S2))
          ss.S2<-cbind(paste(paste("S2.",row.names(ss.S2)[rn.S2],sep=""),": nUnstable, nPeriods, percUnstable",sep=""),	
              as.character(apply(ss.S2,1,function(x)paste(x,collapse=", "))))
          summary.output <- rbind(summary.output,ss.S2)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"S2 Table generated (Percentage of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"
        }
        if(all(x$ss.S3!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"S3 Table(s) generated (Breakdown of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          ss.S3 <- x$ss.S3
          rn.S3 <- c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")
          rn.S3.index <- which(rn.S3%in%names(ss.S3))
          for(i in 1:length(rn.S3.index)){
            table.ss.S3 <- cbind(paste("S3.",rn.S3[rn.S3.index[i]],": ",paste(colnames(ss.S3[[i]]),collapse=", "),sep=""),
                as.character(apply(ss.S3[[i]],1,function(x)paste(x,collapse=", "))))
            summary.output <- rbind(summary.output,table.ss.S3)
          }
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"S3 Table(s) generated (Breakdown of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"	
        }
        
      }else{
        summary.output[dim(summary.output)[1]+1,1]<-"Sliding spans analysis performed"
        summary.output[dim(summary.output)[1],2]<-"FALSE"
      }
    }
    if(history){
      if(length(grep("h.",names(x),fixed=TRUE))>0){
        summary.output[dim(summary.output)[1]+1,1]<-"History analysis performed"
        summary.output[dim(summary.output)[1],2]<-"TRUE"
        ## R1
        if(all(x$h.R1!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R1 Summary table generated (SA series)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R1 <- x$h.R1[[1]]
          h.R1 <- cbind(paste("R1: ",paste(colnames(h.R1),collapse=", "),sep=""),
              apply(h.R1,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R1)
          h.R1 <- x$h.R1[[2]]
          h.R1 <- cbind(paste("R1: total, ",paste(colnames(h.R1),collapse=", "),sep=""),
              paste("total",paste(h.R1,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R1)
          h.R1 <- x$h.R1[[3]]
          h.R1 <- cbind(paste("R1: ",paste(colnames(h.R1),collapse=", "),sep=""),
              apply(h.R1,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R1)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R1 Summary table generated (SA series)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R2
        if(all(x$h.R2!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R2 Summary table generated (period-period changes in SA)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R2 <- x$h.R2[[1]]
          h.R2 <- cbind(paste("R2: ",paste(colnames(h.R2),collapse=", "),sep=""),
              apply(h.R2,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R2)
          h.R2 <- x$h.R2[[2]]
          h.R2 <- cbind(paste("R2: total, ",paste(colnames(h.R2),collapse=", "),sep=""),
              paste("total",paste(h.R2,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R2)
          h.R2 <- x$h.R2[[3]]
          h.R2 <- cbind(paste("R2: ",paste(colnames(h.R2),collapse=", "),sep=""),
              apply(h.R2,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R2)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R2 Summary table generated (period-period changes in SA)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R4
        if(all(x$h.R4!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R4 Summary table generated (Henderson trend component)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R4 <- x$h.R4[[1]]
          h.R4 <- cbind(paste("R4: ",paste(colnames(h.R4),collapse=", "),sep=""),
              apply(h.R4,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R4)
          h.R4 <- x$h.R4[[2]]
          h.R4 <- cbind(paste("R4: total, ",paste(colnames(h.R4),collapse=", "),sep=""),
              paste("total",paste(h.R4,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R4)
          h.R4 <- x$h.R4[[3]]
          h.R4 <- cbind(paste("R4: ",paste(colnames(h.R4),collapse=", "),sep=""),
              apply(h.R4,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R4)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R4 Summary table generated (Henderson trend component))"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R5
        if(all(x$h.R5!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R5 Summary table generated (period-period changes in trend)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R5 <- x$h.R5[[1]]
          h.R5 <- cbind(paste("R5: ",paste(colnames(h.R5),collapse=", "),sep=""),
              apply(h.R5,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R5)
          h.R5 <- x$h.R5[[2]]
          h.R5 <- cbind(paste("R5: total, ",paste(colnames(h.R5),collapse=", "),sep=""),
              paste("total",paste(h.R5,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R5)
          h.R5 <- x$h.R5[[3]]
          h.R5 <- cbind(paste("R5: ",paste(colnames(h.R5),collapse=", "),sep=""),
              apply(h.R5,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R5)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R5 Summary table generated (period-period changes in trend)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R6
        if(all(x$h.R6!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R6 Summary table generated (conc. and proj. seasonal component)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R6 <- x$h.R6[[1]]
          h.R6 <- cbind(paste("R6: ",paste(colnames(h.R6),collapse=", "),sep=""),
              apply(h.R6,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R6)
          h.R6 <- x$h.R6[[2]]
          h.R6 <- cbind(paste("R6: total, ",paste(colnames(h.R6),collapse=", "),sep=""),
              paste("total",paste(h.R6,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R6)
          h.R6 <- x$h.R6[[3]]
          h.R6 <- cbind(paste("R6: ",paste(colnames(h.R6),collapse=", "),sep=""),
              apply(h.R6,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R6)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R6 Summary table generated (conc. and proj. seasonal component)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }		
        ## R7
        if(!is.null(x$h.R7)){
          summary.output[dim(summary.output)[1]+1,1]<-"R7 Table generated (Likelihood stats from estimating model)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R7 <- x$h.R7
          h.R7 <- cbind(paste("R7: ",paste(colnames(h.R7),collapse=", "),sep=""),
              apply(h.R7,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R7)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R7 Table generated (Likelihood stats from estimating model)))"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R8
        if(!is.null(x$h.R8)){
          summary.output[dim(summary.output)[1]+1,1]<-"R8 Table generated (Cum SumSq Fcst Errors at spec leads)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R8 <- x$h.R8[[1]]
          h.R8 <- cbind(paste("R8: ",paste(colnames(h.R8),collapse=", "),sep=""),
              apply(h.R8,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R8)
          h.R8 <- x$h.R8[[2]]
          h.R8 <- cbind(paste("R8: mean, ",paste(colnames(h.R8),collapse=", "),sep=""),
              paste("mean",paste(h.R8,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R8)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R8 Table generated (SumSq Fcst Errors at spec leads)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        
        
      }else{
        summary.output[dim(summary.output)[1]+1,1]<-"History analysis performed"
        summary.output[dim(summary.output)[1],2]<-"FALSE"
      }
    }	
# identify
    if(identify){
      if(!is.null(x$rsd.iac)){
        for(i in 1:length(x$rsd.iac)){
          sig<-rep("",dim(x$rsd.iac[[i]])[1])
          sig[which(x$rsd.iac[[i]]$pval<0.05 & x$rsd.iac[[i]]$df.q>0)]<-"*"
          rsd.iac<-cbind(round(x$rsd.iac[[i]],digits=3),sig)
          summary.output<-rbind(summary.output,cbind(rep(paste("iac_",names(x$rsd.iac)[[i]],": lag, sample.iac, stderr.iac, Ljung-Box.q, df.q, pval",sep=""),dim(rsd.iac)[1]),paste(apply(rsd.iac[,-dim(rsd.iac)[2]],1,paste,collapse=", "),rsd.iac[,dim(rsd.iac)[2]])))
        }
      }
      if(!is.null(x$rsd.ipc)){
        for(i in 1:length(x$rsd.ipc)){
          sig<-rep("",dim(x$rsd.ipc[[i]])[1])
          sig[which(x$rsd.ipc[[i]]$pval<0.05 & x$rsd.ipc[[i]]$df.q>0)]<-"*"
          rsd.ipc<-cbind(round(x$rsd.ipc[[i]],digits=3),sig)
          summary.output<-rbind(summary.output,cbind(rep(paste("ipc_",names(x$rsd.ipc)[[i]],": lag, sample.ipc, stderr.ipc, Ljung-Box.q, df.q, pval",sep=""),dim(rsd.ipc)[1]),paste(apply(rsd.ipc[,-dim(rsd.ipc)[2]],1,paste,collapse=", "),rsd.ipc[,dim(rsd.ipc)[2]])))
        }
      }	
    }	
#		names.sumout<-unique(summary.output[,1])
  }else{
#		cat("\n\tX11 Regression\n\n")
    summary.output[dim(summary.output)[1]+1,1]<-"X11 Regression"
    summary.output[dim(summary.output)[1],2]<-"TRUE"		
    
    summary.output[dim(summary.output)[1]+1,1]<-"Regression Model"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$regmdl))
#		cat("\n\tOutlier Detection\n")
#	if(x$ifout=="Outlier detection performed"){
#		summary.output[dim(summary.output)[1]+1,1]<-"Outlier detection performed"
#		summary.output[dim(summary.output)[1],2]<-paste("TRUE")
    
#			cat("Critical |t| for outliers:\t\n")
    summary.output[dim(summary.output)[1]+1,1]<-"aocrit"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$crit))
    summary.output[dim(summary.output)[1]+1,1]<-"Total Number of Outliers"
    summary.output[dim(summary.output)[1],2]<-paste(length(x$out)+length(x$autooutlier)-length(which(x$out=="-"))-length(which(x$autooutlier=="-")))
    summary.output[dim(summary.output)[1]+1,1]<-"Nr of Automatically Identified Outliers"
    summary.output[dim(summary.output)[1],2]<-paste(length(x$autooutlier)-length(which(x$autooutlier=="-")))
#	}	
#	else{
#		summary.output[dim(summary.output)[1]+1,1]<-"Outlier detection performed"
#		summary.output[dim(summary.output)[1],2]<-paste("FALSE")
#	}
#		cat("\n\tRegression Model\n")
    rest<-unlist(lapply(strsplit(as.character(x$regmdl),"+",fixed=TRUE),function(x)gsub("^\\s+|\\s+$", "",x)))
    rest<-names(x)[which(names(x)%in%rest)]
    liste <- c("outlier","userdefined","leapyear","td",rest,"autooutlier")#,"almostoutlier")
    liste<-liste[which(liste%in%names(x))]
    empty <- which(unlist(lapply(1:length(x),function(y)any(x[[y]]=="-"))))
    res <- as.data.frame(do.call(rbind,lapply(which(!liste %in% names(x[empty])),function(j){
                  if(!any(grepl(names(x[liste[j]]),names(x[[liste[j]]])))){
                    names(x[[liste[j]]])<-paste(names(x[liste[j]]),"_",names(x[[liste[j]]]),sep="")	
                  }
                  do.call(rbind,lapply(1:length(x[[liste[j]]]),function(i){
                            c(names(x[[liste[j]]][i]),unlist(x[[liste[j]]][[i]]))}))})))
    if(all(dim(res))>0){
      res[,2:4] <- apply(res[,2:4],2,function(x)as.numeric(formatC(as.numeric(as.character(x)),digits=3,format="f")))
      colnames(res)[1]<-"variable"
      res2 <- cbind(paste(1:length(res[,1]),"variable, coef, stderr, tval"),apply(res,1,paste,collapse=", "))
      summary.output<-rbind(summary.output,res2)
      if(!is.null(x[["derived.coef"]])){
        summary.output[dim(summary.output)[1]+1,1]<-"* Derived parameter estimates"
        summary.output[dim(summary.output)[1],2]<-paste(x[["derived.coef"]])
      }
    }
#		cat("\n\tSeasonal Adjustment\n\n")
    
    summary.output[dim(summary.output)[1]+1,1]<-"Identifiable Seasonality"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$id.seas))
    if(id.rsdseas){
      summary.output[dim(summary.output)[1]+1,1]<-"Residual Seasonality"
      if(x$id.rsdseas=="none")
        summary.output[dim(summary.output)[1],2]<-"none"
      else
        summary.output[dim(summary.output)[1],2]<-"yes"
    }		
    summary.output[dim(summary.output)[1]+1,1]<-"Seasonal Peaks"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$peaks.seas))
    summary.output[dim(summary.output)[1]+1,1]<-"Trading Day Peaks"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$peaks.td))
    summary.output[dim(summary.output)[1]+1,1]<-"Q Statistic"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$q))
    if(quality.stat){
      summary.output[dim(summary.output)[1]+1,1]<-"Q2 Statistic"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$q2))
      summary.output[dim(summary.output)[1]+1,1]<-"M1"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m1))
      summary.output[dim(summary.output)[1]+1,1]<-"M2"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m2))
      summary.output[dim(summary.output)[1]+1,1]<-"M3"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m3))
      summary.output[dim(summary.output)[1]+1,1]<-"M4"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m4))
      summary.output[dim(summary.output)[1]+1,1]<-"M5"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m5))
      summary.output[dim(summary.output)[1]+1,1]<-"M6"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m6))
      summary.output[dim(summary.output)[1]+1,1]<-"M7"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m7))
      summary.output[dim(summary.output)[1]+1,1]<-"M8"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m8))
      summary.output[dim(summary.output)[1]+1,1]<-"M9"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m9))
      summary.output[dim(summary.output)[1]+1,1]<-"M10"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m10))
      summary.output[dim(summary.output)[1]+1,1]<-"M11"
      summary.output[dim(summary.output)[1],2]<-paste(unlist(x$m11))
      
    }
    summary.output[dim(summary.output)[1]+1,1]<-"Nr of M stats outside limits"
    summary.output[dim(summary.output)[1],2]<-paste(unlist(x$nmfail))
    
    
    if(spectra.detail){
      new.names<-function(z){
        for(i in 1:length(x[[z]])){
          if(length(x[[z]][[i]])>1){
            x[[z]][[i]] <- do.call(paste,x[[z]][[i]])	}}
        y<-as.data.frame(x[[z]],row.names="")			
        colnames(y)<-gsub(paste(z,".",sep=""),replacement="",names(unlist(x[[z]])),fixed=TRUE)
        return(y)}
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of the original series"
      summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcori")),":",new.names("spcori"),collapse=",",sep="")
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of differenced seasonally adjusted series"
      summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcsa")),":",new.names("spcsa"),collapse=",",sep="")
      summary.output[dim(summary.output)[1]+1,1]<-"Spectrum of modified irregular series"
      summary.output[dim(summary.output)[1],2]<-paste(names(new.names("spcirr")),":",new.names("spcirr"),collapse=",",sep="")
    }
    
    
#		cat("\n\tSeasonal and Trend Moving Averages\n\n")
    
    summary.output[dim(summary.output)[1]+1,1]<-"SA decomposition"
    summary.output[dim(summary.output)[1],2]<-paste(x$samode[[length(x$samode)]])
    
    if(x$seasonalma[[1]]=="M.S.R."){
      summary.output[dim(summary.output)[1]+1,1]<-"Seasonal moving average"
      summary.output[dim(summary.output)[1],2]<-paste(x$seasonalma[[length(x$seasonalma)]],"(Based on msr size)",sep=" ")
    }else{
      summary.output[dim(summary.output)[1]+1,1]<-"Seasonal moving average"
      summary.output[dim(summary.output)[1],2]<-paste(x$seasonalma[[length(x$seasonalma)]])
    }
    summary.output[dim(summary.output)[1]+1,1]<-"Henderson filter"
    summary.output[dim(summary.output)[1],2]<-paste(x$trendma[[length(x$trendma)]],"-term",sep="")
    if(slidingspans){
      if(length(grep("ss.",names(x)))>0){
        summary.output[dim(summary.output)[1]+1,1]<-"Sliding spans analysis performed"
        summary.output[dim(summary.output)[1],2]<-"TRUE"
        if(all(x$ss.options!="-")){
          ss.options<-cbind(c("Nr of spans","Length of spans","First period in first span","First year in first span"),as.character(x$ss.options))
          summary.output<-rbind(summary.output,ss.options)	
        }
        if(any(x$ss.seasTests!="-")){
          seasTests <- x$ss.seasTests
          seasTests <- cbind(paste(rownames(seasTests),": M7, Identifiable Seasonality",sep=""), 
              as.character(apply(seasTests[which(colnames(seasTests)%in%c("m7","idseas"))],1,function(x)paste(x,collapse=", "))))
          summary.output <- rbind(summary.output,seasTests)					
        }	
        if(all(x$ss.S1!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"S1 Table generated (Period means of seasonal factors)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          ss.S1 <- x$ss.S1[[1]]
          ss.S1 <- cbind(paste("S1: ",paste(colnames(ss.S1),collapse=", "),sep=""),
              apply(ss.S1,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,ss.S1)
          
          ss.S1 <- x$ss.S1[[2]]
          ss.S1 <- cbind(paste("S1.summary",paste(row.names(ss.S1),paste(colnames(ss.S1),collapse=", "),sep=": "),sep="."),
              as.character(apply(ss.S1,1,function(x)paste(x,collapse=", "))))
          summary.output <- rbind(summary.output,ss.S1)					
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"S1 Table generated (Period means of seasonal factors)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        if(all(x$ss.S2!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"S2 Table generated (Percentage of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          ss.S2 <- x$ss.S2
          rn.S2 <- which(c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")%in%row.names(ss.S2))
          ss.S2<-cbind(paste(paste("S2.",row.names(ss.S2)[rn.S2],sep=""),": nUnstable, nPeriods, percUnstable",sep=""),	
              as.character(apply(ss.S2,1,function(x)paste(x,collapse=", "))))
          summary.output <- rbind(summary.output,ss.S2)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"S2 Table generated (Percentage of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"
        }
        if(all(x$ss.S3!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"S3 Table(s) generated (Breakdown of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          ss.S3 <- x$ss.S3
          rn.S3 <- c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")
          rn.S3.index <- which(rn.S3%in%names(ss.S3))
          for(i in 1:length(rn.S3.index)){
            table.ss.S3 <- cbind(paste("S3.",rn.S3[rn.S3.index[i]],": ",paste(colnames(ss.S3[[i]]),collapse=", "),sep=""),
                as.character(apply(ss.S3[[i]],1,function(x)paste(x,collapse=", "))))
            summary.output <- rbind(summary.output,table.ss.S3)
          }
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"S3 Table(s) generated (Breakdown of unstable periods)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"	
        }
        
      }else{
        summary.output[dim(summary.output)[1]+1,1]<-"Sliding spans analysis performed"
        summary.output[dim(summary.output)[1],2]<-"FALSE"
      }
    }
    if(history){
      if(length(grep("h.",names(x),fixed=TRUE))>0){
        summary.output[dim(summary.output)[1]+1,1]<-"History analysis performed"
        summary.output[dim(summary.output)[1],2]<-"TRUE"
        ## R1
        if(all(x$h.R1!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R1 Summary table generated (SA series)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R1 <- x$h.R1[[1]]
          h.R1 <- cbind(paste("R1: ",paste(colnames(h.R1),collapse=", "),sep=""),
              apply(h.R1,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R1)
          h.R1 <- x$h.R1[[2]]
          h.R1 <- cbind(paste("R1: total, ",paste(colnames(h.R1),collapse=", "),sep=""),
              paste("total",paste(h.R1,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R1)
          h.R1 <- x$h.R1[[3]]
          h.R1 <- cbind(paste("R1: ",paste(colnames(h.R1),collapse=", "),sep=""),
              apply(h.R1,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R1)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R1 Summary table generated (SA series)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R2
        if(all(x$h.R2!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R2 Summary table generated (period-period changes in SA)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R2 <- x$h.R2[[1]]
          h.R2 <- cbind(paste("R2: ",paste(colnames(h.R2),collapse=", "),sep=""),
              apply(h.R2,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R2)
          h.R2 <- x$h.R2[[2]]
          h.R2 <- cbind(paste("R2: total, ",paste(colnames(h.R2),collapse=", "),sep=""),
              paste("total",paste(h.R2,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R2)
          h.R2 <- x$h.R2[[3]]
          h.R2 <- cbind(paste("R2: ",paste(colnames(h.R2),collapse=", "),sep=""),
              apply(h.R2,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R2)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R2 Summary table generated (period-period changes in SA)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R4
        if(all(x$h.R4!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R4 Summary table generated (Henderson trend component)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R4 <- x$h.R4[[1]]
          h.R4 <- cbind(paste("R4: ",paste(colnames(h.R4),collapse=", "),sep=""),
              apply(h.R4,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R4)
          h.R4 <- x$h.R4[[2]]
          h.R4 <- cbind(paste("R4: total, ",paste(colnames(h.R4),collapse=", "),sep=""),
              paste("total",paste(h.R4,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R4)
          h.R4 <- x$h.R4[[3]]
          h.R4 <- cbind(paste("R4: ",paste(colnames(h.R4),collapse=", "),sep=""),
              apply(h.R4,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R4)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R4 Summary table generated (Henderson trend component))"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R5
        if(all(x$h.R5!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R5 Summary table generated (period-period changes in trend)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R5 <- x$h.R5[[1]]
          h.R5 <- cbind(paste("R5: ",paste(colnames(h.R5),collapse=", "),sep=""),
              apply(h.R5,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R5)
          h.R5 <- x$h.R5[[2]]
          h.R5 <- cbind(paste("R5: total, ",paste(colnames(h.R5),collapse=", "),sep=""),
              paste("total",paste(h.R5,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R5)
          h.R5 <- x$h.R5[[3]]
          h.R5 <- cbind(paste("R5: ",paste(colnames(h.R5),collapse=", "),sep=""),
              apply(h.R5,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R5)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R5 Summary table generated (period-period changes in trend)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R6
        if(all(x$h.R6!="-")){
          summary.output[dim(summary.output)[1]+1,1]<-"R6 Summary table generated (conc. and proj. seasonal component)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R6 <- x$h.R6[[1]]
          h.R6 <- cbind(paste("R6: ",paste(colnames(h.R6),collapse=", "),sep=""),
              apply(h.R6,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R6)
          h.R6 <- x$h.R6[[2]]
          h.R6 <- cbind(paste("R6: total, ",paste(colnames(h.R6),collapse=", "),sep=""),
              paste("total",paste(h.R6,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R6)
          h.R6 <- x$h.R6[[3]]
          h.R6 <- cbind(paste("R6: ",paste(colnames(h.R6),collapse=", "),sep=""),
              apply(h.R6,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R6)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R6 Summary table generated (conc. and proj. seasonal component)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }		
        ## R7
        if(!is.null(x$h.R7)){
          summary.output[dim(summary.output)[1]+1,1]<-"R7 Table generated (Likelihood stats from estimating model)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R7 <- x$h.R7
          h.R7 <- cbind(paste("R7: ",paste(colnames(h.R7),collapse=", "),sep=""),
              apply(h.R7,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R7)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R7 Table generated (Likelihood stats from estimating model)))"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        ## R8
        if(!is.null(x$h.R8)){
          summary.output[dim(summary.output)[1]+1,1]<-"R8 Table generated (Cum SumSq Fcst Errors at spec leads)"
          summary.output[dim(summary.output)[1],2]<-"TRUE"
          
          h.R8 <- x$h.R8[[1]]
          h.R8 <- cbind(paste("R8: ",paste(colnames(h.R8),collapse=", "),sep=""),
              apply(h.R8,1,function(x)paste(x,collapse=", ")))
          summary.output <- rbind(summary.output,h.R8)
          h.R8 <- x$h.R8[[2]]
          h.R8 <- cbind(paste("R8: mean, ",paste(colnames(h.R8),collapse=", "),sep=""),
              paste("mean",paste(h.R8,collapse=", "),sep=", "))
          summary.output <- rbind(summary.output,h.R8)
        }else{
          summary.output[dim(summary.output)[1]+1,1]<-"R8 Table generated (SumSq Fcst Errors at spec leads)"
          summary.output[dim(summary.output)[1],2]<-"FALSE"					
        }
        
        
      }else{
        summary.output[dim(summary.output)[1]+1,1]<-"History analysis performed"
        summary.output[dim(summary.output)[1],2]<-"FALSE"
      }
    }	
# identify
    if(identify){
      if(!is.null(x$rsd.iac)){
        for(i in 1:length(x$rsd.iac)){
          sig<-rep("",dim(x$rsd.iac[[i]])[1])
          sig[which(x$rsd.iac[[i]]$pval<0.05 & x$rsd.iac[[i]]$df.q>0)]<-"*"
          rsd.iac<-cbind(round(x$rsd.iac[[i]],digits=3),sig)
          summary.output<-rbind(summary.output,cbind(rep(paste("iac_",names(x$rsd.iac)[[i]],": lag, sample.iac, stderr.iac, Ljung-Box.q, df.q, pval",sep=""),dim(rsd.iac)[1]),paste(apply(rsd.iac[,-dim(rsd.iac)[2]],1,paste,collapse=", "),rsd.iac[,dim(rsd.iac)[2]])))
        }
      }
      if(!is.null(x$rsd.ipc)){
        for(i in 1:length(x$rsd.ipc)){
          sig<-rep("",dim(x$rsd.ipc[[i]])[1])
          sig[which(x$rsd.ipc[[i]]$pval<0.05 & x$rsd.ipc[[i]]$df.q>0)]<-"*"
          rsd.ipc<-cbind(round(x$rsd.ipc[[i]],digits=3),sig)
          summary.output<-rbind(summary.output,cbind(rep(paste("ipc_",names(x$rsd.ipc)[[i]],": lag, sample.ipc, stderr.ipc, Ljung-Box.q, df.q, pval",sep=""),dim(rsd.ipc)[1]),paste(apply(rsd.ipc[,-dim(rsd.ipc)[2]],1,paste,collapse=", "),rsd.ipc[,dim(rsd.ipc)[2]])))
        }
      }	
    }	
    
  }
#		sumout<-unlist(summary.output[,1])
#		names.sumout<-unique(sumout)
#		length.names.sumout<-unlist(lapply(names.sumout,function(x)length(grep(x,sumout))))
  names(summary.output)<-c("DIAGNOSTICS","--- Rout ---")
  spl<-split(summary.output[,1,drop=FALSE],factor(summary.output[,1],levels=unique(summary.output[,1])))		
  ind <- which(sapply(spl, nrow) > 1)		
  v <- lapply(ind, function(x) {data.frame(DIAGNOSTICS=paste(1:nrow(spl[[x]]), names(spl)[x]),stringsAsFactors=FALSE)} )	
  spl[ind] <- v
  new.col<-do.call("rbind",spl)
  summary.output[,1]<-new.col	
  return(summary.output)
#		list(summary.output=summary.output,names.sumout=names.sumout,length.names.sumout=length.names.sumout)
}


