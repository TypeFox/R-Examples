summary.x12work <- function(object,fullSummary=FALSE,spectra.detail=FALSE,almostout=FALSE,rsd.autocorr=NULL,quality.stat=FALSE,likelihood.stat=FALSE,aape=FALSE,id.rsdseas=FALSE,slidingspans=FALSE,history=FALSE,identify=FALSE,...){
  summaryworkhorse(object$dg,fullSummary=fullSummary,spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,identify=identify) 
}
summaryworkhorse <- function(x,fullSummary=FALSE,spectra.detail=FALSE,almostout=FALSE,rsd.autocorr=NULL,quality.stat=FALSE,likelihood.stat=FALSE,aape=FALSE,id.rsdseas=FALSE,slidingspans=FALSE,history=FALSE,identify=FALSE){
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
	cat("\n\tTime Series\n\n")
	cat("Frequency:",x$frequency,"\n")	
	if(length(x$span)>1){
	span <- str_trim(unlist(strsplit(x$span,":")))
	span.index <- which(span=="span")
		cat("Span:",span[span.index+1],"\n")
  }else{
		cat("Span:",x$span,"\n")}	
  
  span.index <- which(span=="modelspan")
  if(length(span.index)>0)
    modelspan <- span[span.index+1]
  span.index <- which(span=="outlierspan")
  if(length(span.index)>0)
    outlierspan <- span[span.index+1]
### RegARIMA Option:
	
	if(x$x11regress=="no"){	
    cat("\n\tModel Definition\n\n")
    if(x$automdl!="-"){
      cat("ARIMA Model:",unlist(x$arimamdl),"(Automatic Model Choice)\n")	
    }else{
      cat("ARIMA Model:",unlist(x$arimamdl),"\n")#automdl erwaehnen
    }
    #zz <- data.frame("(row names)"= c("aaaaa", "b"), check.names=FALSE)
    if(exists("modelspan"))
      cat("Model Span:",span[span.index+1],"\n")
    if(x$transform=="Automatic selection"){
      cat("Transformation:",unlist(x$transform),":",unlist(x$autotransform),"\n")
    }else{
      cat("Transformation:",unlist(x$transform),"\n")
    }	
    cat("Regression Model:",unlist(x$regmdl),"\n")
    cat("\n\tOutlier Detection\n\n")
    if(x$ifout=="Outlier detection performed"){
      if(exists("outlierspan"))
        cat("Outlier Span:",span[span.index+1],"\n")
      cat("Critical |t| for outliers:\t\n")
      print(unlist(x$crit))
      cat("Total Number of Outliers:",unlist(x$nout),"\n")
      cat("Automatically Identified Outliers:",unlist(x$nautoout),"\n")
      if(almostout){
        cat("Number of ts values that were almost identified as outliers:",unlist(x$nalmostout),"\n")
      }}	
    else{
      cat("No outlier detection performed\n")	
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
      cat("\n\tRegression Model\n")
      print(res)
      if(!is.null(x[["derived.coef"]])){
        cat("* Derived parameter estimates: ",x[["derived.coef"]],"\n")	
      }
      }
	  if(likelihood.stat){
		  cat("\n\tLikelihood Statistics\n")
#		  cat("\n")
#		  lstat<-as.data.frame(matrix(c("Log Likelihood","AIC","AICC","BIC","HQ",x$loglikelihood,x$aic,x$aicc,x$bic,x$hq),ncol=2))
#		  colnames(lstat)<-c("Likelihood Statistics"," ")
#		  print(lstat)
		lstat<-matrix(c("AIC:","AICC:","BIC: ","HQ: ","Log Likelihood:",x$aic,x$aicc,x$bic,x$hq,x$loglikelihood),ncol=2)	  
		  write.table(lstat,quote=FALSE,row.names=FALSE, col.names=FALSE,sep="\t")
#write.table(lstat,quote=FALSE,row.names=FALSE, col.names=FALSE)		  
		  
	  }	  
	  if(aape && length(x$aape)>1){
		cat("\nAverage absolute percentage error\n")
		mode<-ifelse(x$aape$aape.mode=="outofsample","out of sample","within sample")
		cat("\tin ",mode," forecasts",sep="","\n")
		aape.mat<-matrix(c("Last year:  ","Last-1 year:","Last-2 year:","Last 3 years:",x$aape$aape.0,x$aape$aape.1,x$aape$aape.2,x$aape$aape.3),ncol=2)	  
		write.table(aape.mat,quote=FALSE,row.names=FALSE, col.names=FALSE,sep="\t")
		}
	cat("\n\tSeasonal Adjustment\n\n")
    cat("Identifiable Seasonality:",unlist(x$id.seas),"\n")
	if(id.rsdseas){
	if(x$id.rsdseas=="none")
	cat("Residual Seasonality: none\n")		
	else
	cat("Residual Seasonality: yes\n")	
	}		
    cat("Seasonal Peaks:",unlist(x$peaks.seas),"\n")
    cat("Trading Day Peaks:",unlist(x$peaks.td),"\n")
    cat("Overall Index of Quality of SA\n(Acceptance Region from 0 to 1)\nQ:",unlist(x$q),"\n")
    if(quality.stat){
      cat("Q2:",unlist(x$q2),"\n")#(Q Statistic computed w/o the M2 Quality Control Statistic)\n")	
	  cat("Quality Control Statistics\n")
	  cat("M1:",unlist(x$m1),"\n")
	  cat("M2:",unlist(x$m2),"\n")
	  cat("M3:",unlist(x$m3),"\n")
	  cat("M4:",unlist(x$m4),"\n")
	  cat("M5:",unlist(x$m5),"\n")
	  cat("M6:",unlist(x$m6),"\n")
	  cat("M7:",unlist(x$m7),"\n")
	  cat("M8:",unlist(x$m8),"\n")
	  cat("M9:",unlist(x$m9),"\n")
	  cat("M10:",unlist(x$m10),"\n")
	  cat("M11:",unlist(x$m11),"\n")
    }
    cat("Number of M statistics outside the limits:",unlist(x$nmfail),"\n")
	if(spectra.detail){
      new.names<-function(z){
        for(i in 1:length(x[[z]])){
          if(length(x[[z]][[i]])>1){
            x[[z]][[i]] <- do.call(paste,x[[z]][[i]])	}}
        y<-as.data.frame(x[[z]],row.names="")			
        colnames(y)<-gsub(paste(z,".",sep=""),replacement="",names(unlist(x[[z]])),fixed=TRUE)
        return(y)}
      cat("Spectrum of the original series\n")
      print(new.names("spcori"))	
      cat("Spectrum of the regARIMA model residuals\n")
	  if(any(x$spcrsd!="-"))
	  print(new.names("spcrsd"))	
      else
	  cat("\t- Not available -\n")  
	  cat("Spectrum of differenced seasonally adjusted series\n")
      print(new.names("spcsa"))
      cat("Spectrum of modified irregular series\n")
      print(new.names("spcirr"))
    }
    #		cat("\n\tSeasonal and Trend Moving Averages\n\n")
    cat("\nSA decomposition:",x$samode[[length(x$samode)]],"\n")		
    if(x$seasonalma[[1]]=="M.S.R."){
      cat("Seasonal moving average used for the final iteration: \n",x$seasonalma[[length(x$seasonalma)]],
          " (Based on the size of the global moving seasonality ratio (msr))\n",sep="")
    }else{
      cat("Moving average used to estimate the seasonal factors:",x$seasonalma[[length(x$seasonalma)]],"\n")	
    }
    cat("Moving average used to estimate the final trend-cycle: ",x$trendma[[length(x$trendma)]],"-term Henderson filter\n",sep="")
    if(!is.null(rsd.autocorr)){
      #rsd.autocorr=c("rsd.acf","rsd.pacf","rsd.acf2")
      if("acf"%in%rsd.autocorr){
        cat("\n\tSample Autocorrelations of the Residuals\n")
		if(!is.null(x$rsd.acf)){			
        #cat("p-values approximate the probability of observing a q-value at least this
        #large when the model fitted is correct.")
        #!!Diese Schranke (<0.05) kann sich noch aendern falls pickmodel spec implementiert wird  
        cat("(Small p-values (<0.05) indicate model inadequacy (for df.q >0))\n\n")
        sig<-rep("",dim(x$rsd.acf)[1])
        sig[which(x$rsd.acf$pval<0.05 & x$rsd.acf$df.q>0)]<-"*"
        rsd.acf<-cbind(round(x$rsd.acf,digits=3),sig)
        colnames(rsd.acf)<-c(colnames(x$rsd.acf),"")
        print(rsd.acf)}else{
		cat("\n\t- Not available -\n\n")
		}
      }
      if("pacf"%in%rsd.autocorr){
        cat("\n\tSample Partial Autocorrelations of the Residuals\n\n")
		if(!is.null(x$rsd.pacf))
		print(round(x$rsd.pacf,digits=3))	
		else
		cat("\t- Not available -\n\n")
      }
      if("acf2"%in%rsd.autocorr){
        cat("\n\tSample Autocorrelations of the Squared Residuals\n\n")
        if(!is.null(x$rsd.acf2))
		print(round(x$rsd.acf2,digits=3))
		else
		cat("\t- Not available -\n\n")
      }
    }
	
#	if(slidingspans){
#	if(length(grep("ss.",names(x)))>0){
#		if(!is.null(x$ss.))
#		
##		if(length(grep("rsd.acf",names(x)))>0)
#			
##		rsd.autocorr=c("rsd.acf","rsd.pacf","rsd.acf2")
##		if("acf"%in%rsd.autocorr){
##			cat("\n\tSample Autocorrelations of the Residuals\n")
##			if(!is.null(x$rsd.acf)){
#write.table(ss.options,quote=FALSE,row.names=FALSE,sep="\t")
#	}}
	if(slidingspans){
		cat("\n\tSlidingspans\n\n")	
		if(length(grep("ss.",names(x),fixed=TRUE))>0){
## S0
			cat("Summary of options selected for this run:")	
			if(all(x$ss.options!="-")){
			ss.options<-as.data.frame((t(as.data.frame(x$ss.options))))
			row.names(ss.options)<-c("Number of spans:","Length of spans:","First period in first span:","First year in first span:")
			colnames(ss.options)<-""
			print(ss.options)
			}else{cat("\t- Not available -\n")}
			
			cat("\nSeasonality:\n")	
			if(any(x$ss.seasTests!="-")){
#			F statistics not printed (might change this if requested)			
#			cat("stabseas (F-Value of test for the presence of seasonality assuming stability)\n",
#			"movseas (F-Value of moving seasonality test)\n")		
			seasTests <- x$ss.seasTests
			seasTests <- seasTests[which(colnames(seasTests)%in%c("m7","idseas"))]
			colnames(seasTests)<-gsub("m7","M7",colnames(seasTests))
			colnames(seasTests)<-gsub("idseas","Identifiable Seasonality",colnames(seasTests))
			print(seasTests)}else{cat("\n\t- Not available -\n")}
## S1
			cat("\nS1 Period means of Seasonal Factors\n(Movements within a period should be small)\n")
			if(all(x$ss.S1!="-")){
			print(x$ss.S1[[1]])
		
			cat("\nSummary statistics for mean seasonal factor\n")	
			print(x$ss.S1[[2]])
			}else{cat("\t- Not available -\n")}
## S2		
			cat("\nS2 Percentage of unstable periods\n")	
			if(all(x$ss.S2!="-")){
			ss.S2 <- x$ss.S2
			rn.S2 <- which(c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")%in%row.names(ss.S2))
			row.names(ss.S2) <- c("a.Seasonal Factors","b.Trading Days","c.Final SA series","d.Period to period changes","e.Year on year changes")[rn.S2]		
			print(ss.S2)
			}else{cat("\t- Not available -\n")}
## S3			
			cat("\nS3 Breakdown of unstable periods\n")	
			if(all(x$ss.S3!="-")){
			ss.S3 <- x$ss.S3
			rn.S3 <- which(c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")%in%names(ss.S3))
			new.rn.S3 <- c("a.Seasonal Factors","b.Trading Days","c.Final SA series","d.Period to period changes","e.Year on year changes")
			for(i in 1:length(rn.S3)){
				cat("\n",new.rn.S3[rn.S3[i]],"\n")
				print(x$ss.S3[[i]])
			}
			}else{cat("\t- Not available -\n\n")}
		}else{cat("\t- Not available -\n\n")}
	}
if(history){
	cat("\n\tHistory\n")	
	if(length(grep("h.",names(x),fixed=TRUE))>0){
## R1
		cat("\nR1 Average absolute percent revisions of the seasonal adjusments\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
		if(all(x$h.R1!="-")){	
			x$h.R1[[1]][,-1] <- round(x$h.R1[[1]][,-1],digits=3)
			x$h.R1[[2]] <- round(x$h.R1[[2]],digits=3)						
			x$h.R1[[3]][,-1] <- round(x$h.R1[[3]][,-1],digits=3)						
			print(rbind(x$h.R1[[1]],
					rep("",dim(x$h.R1[[1]])[2]),		
					c("total",unlist(x$h.R1[[2]])),rep("",dim(x$h.R1[[1]])[2])))
			print(x$h.R1[[3]])
		}else{cat("\t- Not available -\n")}
## R2
		cat("\nR2 Average absolute percent revisions of the period-to-period percent change of the ajusments\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
		if(all(x$h.R2!="-")){	
			x$h.R2[[1]][,-1] <- round(x$h.R2[[1]][,-1],digits=3)
			x$h.R2[[2]] <- round(x$h.R2[[2]],digits=3)						
			x$h.R2[[3]][,-1] <- round(x$h.R2[[3]][,-1],digits=3)						
			print(rbind(x$h.R2[[1]],
							rep("",dim(x$h.R2[[1]])[2]),		
							c("total",unlist(x$h.R2[[2]])),rep("",dim(x$h.R2[[1]])[2])))
			print(x$h.R2[[3]])
		}else{cat("\t- Not available -\n")}
## R4
		cat("\nR4 Average absolute percent revisions of the final Henderson trend component\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
		if(all(x$h.R4!="-")){	
			x$h.R4[[1]][,-1] <- round(x$h.R4[[1]][,-1],digits=3)
			x$h.R4[[2]] <- round(x$h.R4[[2]],digits=3)						
			x$h.R4[[3]][,-1] <- round(x$h.R4[[3]][,-1],digits=3)						
			print(rbind(x$h.R4[[1]],
							rep("",dim(x$h.R4[[1]])[2]),		
							c("total",unlist(x$h.R4[[2]])),rep("",dim(x$h.R4[[1]])[2])))
			print(x$h.R4[[3]])
		}else{cat("\t- Not available -\n")}
## R5
		cat("\nR5 Average absolute percent revisions of period-to-period percent change of the trend-cycle\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
		if(all(x$h.R5!="-")){
			x$h.R5[[1]][,-1] <- round(x$h.R5[[1]][,-1],digits=3)
			x$h.R5[[2]] <- round(x$h.R5[[2]],digits=3)						
			x$h.R5[[3]][,-1] <- round(x$h.R5[[3]][,-1],digits=3)						
			print(rbind(x$h.R5[[1]],
							rep("",dim(x$h.R5[[1]])[2]),		
							c("total",unlist(x$h.R5[[2]])),rep("",dim(x$h.R5[[1]])[2])))
			print(x$h.R5[[3]])
		}else{cat("\t- Not available -\n")}
## R6
		cat("\nR6 Average absolute percent revisions of the concurrent and projected seasonal component\n")
		if(all(x$h.R6!="-")){	
			x$h.R6[[1]][,-1] <- round(x$h.R6[[1]][,-1],digits=3)
			x$h.R6[[2]] <- round(x$h.R6[[2]],digits=3)						
			x$h.R6[[3]][,-1] <- round(x$h.R6[[3]][,-1],digits=3)			
			print(rbind(x$h.R6[[1]],
							rep("",dim(x$h.R6[[1]])[2]),		
							c("total",unlist(x$h.R6[[2]]),rep("",dim(x$h.R6[[1]])[2]))))
			print(x$h.R6[[3]])
		}else{cat("\t- Not available -\n")}		
## R7
		cat("\nR7 Likelihood statistics from estimating regARIMA model")
		if(!is.null(x$h.R7)){
		cat(" over spans with ending dates",paste(substr(x$h.R7[1,1],1,4),".",substr(x$h.R7[1,1],5,6),sep=""),"to",paste(substr(x$h.R7[dim(x$h.R7)[1],1],1,4),".",substr(x$h.R7[dim(x$h.R7)[1],1],5,6),sep=""),"\n")
			x$h.R7[dim(x$h.R7)[1],1]
			print(x$h.R7)
		}else{cat("\n\t- Not available -\n")}
## R8
		cat("\nR8 Sum of squared forecast errors at specified leads from the end of each span")
		if(!is.null(x$h.R8)){					
			R8 <- x$h.R8[[1]]
			colnames(R8) <-gsub("sumSqFcstError","lead",colnames(R8))
			colnames(R8) <- gsub("date","forecast.date",colnames(R8))
			cat("\nCumulative sum of squared forecast errors\n")			
			print(R8)
			cat("\nMean sum of squared forecast errors\n")			
			print(x$h.R8[[2]])
		}else{cat("\n\t- Not available -\n")}
	}else{cat("\n\t- Not available -\n\n")}
}

if(identify){
	cat("\n\tModel Identification\n")	
#iac	
	cat("\nSample Autocorrelations of the Residuals\n")
	if(!is.null(x$rsd.iac)){	
		cat("(Small p-values (<0.05) indicate model inadequacy (for df.q >0))\n\n")
		for(i in 1:length(x$rsd.iac)){
			sig<-rep("",dim(x$rsd.iac[[i]])[1])
			sig[which(x$rsd.iac[[i]]$pval<0.05 & x$rsd.iac[[i]]$df.q>0)]<-"*"
			rsd.iac<-cbind(round(x$rsd.iac[[i]],digits=3),sig)
			colnames(rsd.iac)<-c(colnames(x$rsd.iac[[i]]),"")
			cat("\n",names(x$rsd.iac)[[i]],"\n")
			print(rsd.iac)
		}
	}else{cat("\t- Not available -\n")}	
#ipc
	cat("\nSample Partial Autocorrelations of the Residuals\n")
	if(!is.null(x$rsd.ipc)){	
		cat("(Small p-values (<0.05) indicate model inadequacy (for df.q >0))\n\n")
		for(i in 1:length(x$rsd.ipc)){
			sig<-rep("",dim(x$rsd.ipc[[i]])[1])
			sig[which(x$rsd.ipc[[i]]$pval<0.05 & x$rsd.ipc[[i]]$df.q>0)]<-"*"
			rsd.ipc<-cbind(round(x$rsd.ipc[[i]],digits=3),sig)
			colnames(rsd.ipc)<-c(colnames(x$rsd.ipc[[i]]),"")
			cat("\n",names(x$rsd.ipc)[[i]],"\n")
			print(rsd.ipc)
		}
	}else{cat("\t- Not available -\n")}		
}

}else{
#End RegARIMA Option
# x11Regression Option:	  
    cat("\n\tX11 Regression\n\n")
    cat("Regression Model:",unlist(x$regmdl),"\n")
    cat("\n\tOutlier Detection\n")
    cat("Critical |t| for outliers:",unlist(x$crit),"\t\n")
    cat("Total Number of Outliers:",length(x$out)+length(x$autooutlier)-length(which(x$out=="-"))-length(which(x$autooutlier=="-")),"\n")
    cat("Automatically Identified Outliers:",length(x$autooutlier)-length(which(x$autooutlier=="-")),"\n")
    cat("\n\tRegression Model\n")
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
    res[,2:4] <- apply(res[,2:4],2,function(x)as.numeric(formatC(as.numeric(as.character(x)),digits=3,format="f")))
    colnames(res)[1]<-"variable"
    print(res)
    if(!is.null(x[["derived.coef"]])){
      cat("* Derived parameter estimates: ",x[["derived.coef"]],"\n")	
    }
    cat("\n\tSeasonal Adjustment\n\n")
    cat("Identifiable Seasonality:",unlist(x$id.seas),"\n")
	if(id.rsdseas){
		if(x$id.rsdseas=="none")
			cat("Residual Seasonality: none\n")		
		else
			cat("Residual Seasonality: yes\n")	
	}		
    cat("Seasonal Peaks:",unlist(x$peaks.seas),"\n")
    cat("Trading Day Peaks:",unlist(x$peaks.td),"\n")
    cat("Overall Index of Quality of SA\n(Acceptance region from 0 to 1)\nQ:",unlist(x$q),"\n")
    if(quality.stat){
      cat("Q2:",unlist(x$q2),"(Q statistic computed w/o the M2 quality control statistic)\n")	
	  cat("Quality Control Statistics\n")
	  cat("M1:",unlist(x$m1),"\n")
	  cat("M2:",unlist(x$m2),"\n")
	  cat("M3:",unlist(x$m3),"\n")
	  cat("M4:",unlist(x$m4),"\n")
	  cat("M5:",unlist(x$m5),"\n")
	  cat("M6:",unlist(x$m6),"\n")
	  cat("M7:",unlist(x$m7),"\n")
	  cat("M8:",unlist(x$m8),"\n")
	  cat("M9:",unlist(x$m9),"\n")
	  cat("M10:",unlist(x$m10),"\n")
	  cat("M11:",unlist(x$m11),"\n")
	  
  }
    cat("Number of M statistics outside the limits:",unlist(x$nmfail),"\n")
    if(spectra.detail){
      new.names<-function(z){
        for(i in 1:length(x[[z]])){
          if(length(x[[z]][[i]])>1){
            x[[z]][[i]] <- do.call(paste,x[[z]][[i]])	}}
        y<-as.data.frame(x[[z]],row.names="")			
        colnames(y)<-gsub(paste(z,".",sep=""),replacement="",names(unlist(x[[z]])),fixed=TRUE)
        return(y)}
      cat("Spectrum of the original series\n")
      print(new.names("spcori"))	
      cat("Spectrum of differenced seasonally adjusted series\n")
      print(new.names("spcsa"))
      cat("Spectrum of modified irregular series\n")
      print(new.names("spcirr"))
    }
    #		cat("\n\tSeasonal and Trend Moving Averages\n\n")
    cat("\nSA decomposition:",x$samode[[length(x$samode)]],"\n")		
    if(x$seasonalma[[1]]=="M.S.R."){
      cat("Seasonal moving average used for the final iteration: \n",x$seasonalma[[length(x$seasonalma)]],
          " (Based on the size of the global moving seasonality ratio (msr))\n",sep="")
    }else{
      cat("Moving average used to estimate the seasonal factors:",x$seasonalma[[length(x$seasonalma)]],"\n")	
    }
    cat("Moving average used to estimate the final trend-cycle: ",x$trendma[[length(x$trendma)]],"-term Henderson filter\n",sep="")
	if(slidingspans){
		cat("\n\tSlidingspans\n\n")	
		if(length(grep("ss.",names(x),fixed=TRUE))>0){
			## S0
			cat("Summary of options selected for this run:")	
			if(all(x$ss.options!="-")){
				ss.options<-as.data.frame((t(as.data.frame(x$ss.options))))
				row.names(ss.options)<-c("Number of spans:","Length of spans:","First period in first span:","First year in first span:")
				colnames(ss.options)<-""
				print(ss.options)
			}else{cat("\t- Not available -\n")}
			
			cat("\nSeasonality:\n")	
			if(any(x$ss.seasTests!="-")){
#			F statistics not printed (might change this if requested)			
#			cat("stabseas (F-Value of test for the presence of seasonality assuming stability)\n",
#			"movseas (F-Value of moving seasonality test)\n")		
				seasTests <- x$ss.seasTests
				seasTests <- seasTests[which(colnames(seasTests)%in%c("m7","idseas"))]
				colnames(seasTests)<-gsub("m7","M7",colnames(seasTests))
				colnames(seasTests)<-gsub("idseas","Identifiable Seasonality",colnames(seasTests))
				print(seasTests)}else{cat("\n\t- Not available -\n")}
			## S1
			cat("\nS1 Period means of Seasonal Factors\n(Movements within a period should be small)\n")
			if(all(x$ss.S1!="-")){
				print(x$ss.S1[[1]])
				
				cat("\nSummary statistics for mean seasonal factor\n")	
				print(x$ss.S1[[2]])
			}else{cat("\t- Not available -\n")}
			## S2		
			cat("\nS2 Percentage of unstable periods\n")	
			if(all(x$ss.S2!="-")){
				ss.S2 <- x$ss.S2
				rn.S2 <- which(c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")%in%row.names(ss.S2))
				row.names(ss.S2) <- c("a.Seasonal Factors","b.Trading Days","c.Final SA series","d.Period to period changes","e.Year on year changes")[rn.S2]		
				print(ss.S2)
			}else{cat("\t- Not available -\n")}
			## S3			
			cat("\nS3 Breakdown of unstable periods\n")	
			if(all(x$ss.S3!="-")){
				ss.S3 <- x$ss.S3
				rn.S3 <- which(c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")%in%names(ss.S3))
				new.rn.S3 <- c("a.Seasonal Factors","b.Trading Days","c.Final SA series","d.Period to period changes","e.Year on year changes")
				for(i in 1:length(rn.S3)){
					cat("\n",new.rn.S3[rn.S3[i]],"\n")
					print(x$ss.S3[[i]])
				}
			}else{cat("\t- Not available -\n\n")}
		}else{cat("\t- Not available -\n\n")}
	}
	if(history){
		cat("\n\tHistory\n")	
		if(length(grep("h.",names(x),fixed=TRUE))>0){
			## R1
			cat("\nR1 Average absolute percent revisions of the seasonal adjusments\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
			if(all(x$h.R1!="-")){	
				x$h.R1[[1]][,-1] <- round(x$h.R1[[1]][,-1],digits=3)
				x$h.R1[[2]] <- round(x$h.R1[[2]],digits=3)						
				x$h.R1[[3]][,-1] <- round(x$h.R1[[3]][,-1],digits=3)			
				
				print(rbind(x$h.R1[[1]],
								rep("",dim(x$h.R1[[1]])[2]),		
								c("total",unlist(x$h.R1[[2]])),rep("",dim(x$h.R1[[1]])[2])))
				print(x$h.R1[[3]])
			}else{cat("\t- Not available -\n")}
			## R2
			cat("\nR2 Average absolute percent revisions of the period-to-period percent change of the ajusments\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
			if(all(x$h.R2!="-")){	
				x$h.R2[[1]][,-1] <- round(x$h.R2[[1]][,-1],digits=3)
				x$h.R2[[2]] <- round(x$h.R2[[2]],digits=3)						
				x$h.R2[[3]][,-1] <- round(x$h.R2[[3]][,-1],digits=3)			
				
				print(rbind(x$h.R2[[1]],
								rep("",dim(x$h.R2[[1]])[2]),		
								c("total",unlist(x$h.R2[[2]])),rep("",dim(x$h.R2[[1]])[2])))
				print(x$h.R2[[3]])
			}else{cat("\t- Not available -\n")}
			## R4
			cat("\nR4 Average absolute percent revisions of the final Henderson trend component\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
			if(all(x$h.R4!="-")){	
				x$h.R4[[1]][,-1] <- round(x$h.R4[[1]][,-1],digits=3)
				x$h.R4[[2]] <- round(x$h.R4[[2]],digits=3)						
				x$h.R4[[3]][,-1] <- round(x$h.R4[[3]][,-1],digits=3)			
				
				print(rbind(x$h.R4[[1]],
								rep("",dim(x$h.R4[[1]])[2]),		
								c("total",unlist(x$h.R4[[2]])),rep("",dim(x$h.R4[[1]])[2])))
				print(x$h.R4[[3]])
			}else{cat("\t- Not available -\n")}
			## R5
			cat("\nR5 Average absolute percent revisions of period-to-period percent change of the trend-cycle\n(Deviation from",unlist(x$h.target),"estimate defines revisions of seasonal adjustments calculated at lags)\n")
			if(all(x$h.R5!="-")){
				x$h.R5[[1]][,-1] <- round(x$h.R5[[1]][,-1],digits=3)
				x$h.R5[[2]] <- round(x$h.R5[[2]],digits=3)						
				x$h.R5[[3]][,-1] <- round(x$h.R5[[3]][,-1],digits=3)			
				
				print(rbind(x$h.R5[[1]],
								rep("",dim(x$h.R5[[1]])[2]),		
								c("total",unlist(x$h.R5[[2]])),rep("",dim(x$h.R5[[1]])[2])))
				print(x$h.R5[[3]])
			}else{cat("\t- Not available -\n")}
			## R6
			cat("\nR6 Average absolute percent revisions of the concurrent and projected seasonal component\n")
			if(all(x$h.R6!="-")){	
				x$h.R6[[1]][,-1] <- round(x$h.R6[[1]][,-1],digits=3)
				x$h.R6[[2]] <- round(x$h.R6[[2]],digits=3)						
				x$h.R6[[3]][,-1] <- round(x$h.R6[[3]][,-1],digits=3)			
				print(rbind(x$h.R6[[1]],
								rep("",dim(x$h.R6[[1]])[2]),		
								c("total",unlist(x$h.R6[[2]]),rep("",dim(x$h.R6[[1]])[2]))))
				print(x$h.R6[[3]])
			}else{cat("\t- Not available -\n")}		
			## R7
			cat("\nR7 Likelihood statistics from estimating regARIMA model")
			if(!is.null(x$h.R7)){
				cat(" over spans with ending dates",paste(substr(x$h.R7[1,1],1,4),".",substr(x$h.R7[1,1],5,6),sep=""),"to",paste(substr(x$h.R7[dim(x$h.R7)[1],1],1,4),".",substr(x$h.R7[dim(x$h.R7)[1],1],5,6),sep=""),"\n")
				x$h.R7[dim(x$h.R7)[1],1]
				print(x$h.R7)
			}else{cat("\n\t- Not available -\n")}
			## R8
			cat("\nR8 Sum of squared forecast errors at specified leads from the end of each span")
			if(!is.null(x$h.R8)){					
				R8 <- x$h.R8[[1]]
				colnames(R8) <-gsub("sumSqFcstError","lead",colnames(R8))
				colnames(R8) <- gsub("date","forecast.date",colnames(R8))
				cat("\nCumulative sum of squared forecast errors\n")			
				print(R8)
				cat("\nMean sum of squared forecast errors\n")			
				print(x$h.R8[[2]])
			}else{cat("\n\t- Not available -\n")}
		}else{cat("\n\t- Not available -\n\n")}
	}
	if(identify){
		cat("\n\tModel Identification\n")	
#iac	
		cat("\nSample Autocorrelations of the Residuals\n")
		if(!is.null(x$rsd.iac)){	
			cat("(Small p-values (<0.05) indicate model inadequacy (for df.q >0))\n\n")
			for(i in 1:length(x$rsd.iac)){
				sig<-rep("",dim(x$rsd.iac[[i]])[1])
				sig[which(x$rsd.iac[[i]]$pval<0.05 & x$rsd.iac[[i]]$df.q>0)]<-"*"
				rsd.iac<-cbind(round(x$rsd.iac[[i]],digits=3),sig)
				colnames(rsd.iac)<-c(colnames(x$rsd.iac[[i]]),"")
				cat("\n",names(x$rsd.iac)[[i]],"\n")
				print(rsd.iac)
			}
		}else{cat("\t- Not available -\n")}	
#ipc
		cat("\nSample Partial Autocorrelations of the Residuals\n")
		if(!is.null(x$rsd.ipc)){	
			cat("(Small p-values (<0.05) indicate model inadequacy (for df.q >0))\n\n")
			for(i in 1:length(x$rsd.ipc)){
				sig<-rep("",dim(x$rsd.ipc[[i]])[1])
				sig[which(x$rsd.ipc[[i]]$pval<0.05 & x$rsd.ipc[[i]]$df.q>0)]<-"*"
				rsd.ipc<-cbind(round(x$rsd.ipc[[i]],digits=3),sig)
				colnames(rsd.ipc)<-c(colnames(x$rsd.ipc[[i]]),"")
				cat("\n",names(x$rsd.ipc)[[i]],"\n")
				print(rsd.ipc)
			}
		}else{cat("\t- Not available -\n")}		
}
	
	}}


