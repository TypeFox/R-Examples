# Underlying S3 function
x12work <- function(tso,period=frequency(tso),file="Rout",
    series.span=NULL,series.modelspan=NULL,
    transform.function="auto",transform.power=NULL,transform.adjust=NULL,
    regression.variables=NULL,regression.user=NULL,regression.file=NULL,
    regression.usertype=NULL,regression.centeruser=NULL,regression.start=NULL,
    regression.aictest=NULL,
    outlier.types=NULL,outlier.critical=NULL,outlier.span=NULL,outlier.method=NULL,
    identify=FALSE,identify.diff=NULL,identify.sdiff=NULL,identify.maxlag=NULL,
    arima.model=NULL,arima.smodel=NULL,arima.ar=NULL,arima.ma=NULL,
    automdl=FALSE,automdl.acceptdefault=FALSE,automdl.balanced=TRUE,
    automdl.maxorder=c(3,2),automdl.maxdiff=c(1,1),
    forecast_years=NULL,backcast_years=NULL,forecast_conf=.95,
    estimate=FALSE,estimate.outofsample=TRUE,
    check=TRUE,check.maxlag=NULL,
    slidingspans=FALSE,
    slidingspans.fixmdl=NULL,slidingspans.fixreg=NULL,
    slidingspans.length=NULL,slidingspans.numspans=NULL,
    slidingspans.outlier=NULL,
    slidingspans.additivesa=NULL,slidingspans.start=NULL,
    history=FALSE,
    history.estimates=NULL,history.fixmdl=FALSE,
    history.fixreg=NULL,history.outlier=NULL,
    history.sadjlags=NULL,history.trendlags=NULL,history.start=NULL,history.target=NULL,
    x11.sigmalim=c(1.5,2.5),x11.type=NULL,x11.sfshort=FALSE,x11.samode=NULL,
    x11.seasonalma=NULL,x11.trendma=NULL,
    x11.appendfcst=TRUE,x11.appendbcst=FALSE,x11.calendarsigma=NULL,
    x11.excludefcst=TRUE,x11.final="user",
    x11regression=FALSE,
    tblnames=NULL,Rtblnames=NULL,
    x12path=NULL,x13path=NULL,use="x12",
    keep_x12out=TRUE,showWarnings=TRUE
){
  ### Quick Fix: Rename the parameters to previous version:
  seats=FALSE 
  seatsparameter=NULL
  span <-	series.span
  modelspan <-series.modelspan
  transform<-transform.function	
  regvariables <- regression.variables
  reguser <- regression.user
  regfile <- regression.file
  usertype <- regression.usertype
  centeruser <- regression.centeruser
  regfilestart <- regression.start
  aictest<-regression.aictest
#	outlier.detection <- outlier
  outlier <- outlier.types
  critical <- outlier.critical
  outlier_span <- outlier.span
  outlier_method <- outlier.method	
  arima <- arima.model
  sarima <- arima.smodel
  acceptdefault <- automdl.acceptdefault
  balanced <- automdl.balanced
  maxorder <- automdl.maxorder
  maxdiff	<- automdl.maxdiff
  estOutofsample <- estimate.outofsample
  sigmalim <- x11.sigmalim
  onlytd <- x11.type 
  sfshort <- x11.sfshort
  samode <- x11.samode
  seasonalma <- x11.seasonalma
  trendma <- x11.trendma
  x11appendfcst <- x11.appendfcst
  x11appendbcst <- x11.appendbcst
  x11calendarsigma <- x11.calendarsigma
  x11excludefcst <- x11.excludefcst
  x11final <- x11.final
  x11regress <- x11regression
  basename(file)
  
  ext <- c("out","err","spc","otl","dat","a1","b1","d10","d11","d12","ftr","log")
  for(e in ext){
    f <- paste(basename(file),".",e,sep="")
    if(file.exists(f))
      file.remove(f)
  }
  dirgra <- paste("gra_",gsub("\\.","_",basename(file)),sep="")
  unlink(paste(dirname(file),"/",dirgra,sep=""),recursive=TRUE)
  if((length(tso)/period)>15 && !is.null(backcast_years) && !showWarnings){
    cat("\nWarning: x12 cannot produce backcasts for time series that are more than 15 years long!\n")
  }
  
  
  header <- vector()
  header[length(header)+1] <- "series{"
  header[length(header)+1] <- "save=(a1 b1)"
  header[length(header)+1] <- 'title="R Output for x12a"'
  header[length(header)+1] <- paste("start=",paste(start(tso),collapse="."),sep="")
  if(!is.null(span)){
    topaste<-span
    tocollapse<-c(".",".")
    if(any(is.na(span))){
      topaste[which(is.na(span))]<-""
      tocollapse[which(is.na(span))[2]/2]<-""
    }
    header[length(header)+1] <- paste("span=(",paste(topaste[1:2],collapse=tocollapse[1]),",",paste(topaste[3:4],collapse=tocollapse[2]),")",sep="")
  }
  if(!is.null(modelspan)){
    topaste<-modelspan
    tocollapse<-c(".",".")
    if(any(is.na(modelspan))){
      topaste[which(is.na(modelspan))]<-""
      tocollapse[which(is.na(modelspan))[2]/2]<-""
    }
    header[length(header)+1] <- paste("modelspan=(",paste(topaste[1:2],collapse=tocollapse[1]),",",paste(topaste[3:4],collapse=tocollapse[2]),")",sep="")
  }
#	  if(!is.null(series.comptype)){
#		  header[length(header)+1] <- paste("comptype=",series.comptype,sep="") 
#	  }
#	  if(!is.null(series.compwt)){
#		  header[length(header)+1] <- paste("compwt=",series.compwt,sep="") 
#	  }
  
  header[length(header)+1] <- paste("period=",period,sep="")
  
#	  if(!is.null(series.type)){
#		  header[length(header)+1] <- paste("type=",series.type,sep="") 
#	  }
#ERROR:  Argument name "type" not found 
  header[length(header)+1] <- "DECIMALS=5"
  header[length(header)+1] <- paste("file = \"",file,".dat\"",sep="")
  tsoout <- as.character(round(as.vector(tso),5))
  tsoout[is.na(tsoout)] <- "-99999"
  write(tsoout,file=paste(file,".dat",sep=""),ncolumns =1)
  #header[length(header)+1] <- "data=("
  #datarows<-as.vector(tso)
  #datarows[length(datarows)+1] <- ")"
  #datarows[length(datarows)+1] <- "}"
  header[length(header)+1] <- "}"
  addcommands <- vector()
  if(!x11regress){#transform ausschalten falls x11 Regression
    addcommands[length(addcommands)+1] <- paste("transform{") 
    if(is.null(transform.power))
      addcommands[length(addcommands)+1] <- paste("function=",transform,sep="") 
    else
      addcommands[length(addcommands)+1] <- paste("power=",transform.power,sep="") 
    if(!is.null(transform.adjust))
      addcommands[length(addcommands)+1] <- paste("adjust=",transform.adjust,sep="") 
    
    addcommands[length(addcommands)+1] <- "}"
  }
  if(!is.null(c(arima,sarima,arima.ar,arima.ma))&&!x11regress){
    arima <- paste("(",paste(arima,collapse=","),")",sep="")
    if(!is.null(sarima))  
      sarima <- paste("(",paste(sarima,collapse=","),")",sep="")
    addcommands[length(addcommands)+1] <- paste("arima{")
    addcommands[length(addcommands)+1] <- paste("model=",arima,sarima,sep="")
    if(!is.null(arima.ar)){
      arima.ar[is.na(arima.ar)]<-" "
      addcommands[length(addcommands)+1] <- paste("ar=",paste("(",paste(arima.ar,collapse=","),")",sep=""),sep="") 
    }
    if(!is.null(arima.ma)){
      arima.ma[is.na(arima.ma)]<-" "
      addcommands[length(addcommands)+1] <- paste("ma=",paste("(",paste(arima.ma,collapse=","),")",sep=""),sep="") 
    }
    addcommands[length(addcommands)+1] <- "}"
  }
  if(!is.null(c(arima,sarima,arima.ar,arima.ma))&&automdl&&!x11regress)
    cat("Warning: 'automdl' is ignored because an ARIMA model has been specified! \n")
  #cat("Arima and Sarima model specifications are ignored, because automdl is activated! \n")
  if(any(!is.null(c(regvariables,reguser,regfile,aictest,regfilestart,usertype,centeruser))) &&! x11regress){
    addcommands[length(addcommands)+1] <- "regression{"
    addcommands[length(addcommands)+1] <- "save=otl"
    if(!is.null(regvariables))
      addcommands[length(addcommands)+1] <- paste("variables=(",paste(regvariables,collapse=" "),")",sep="")
    if(!is.null(aictest))
      addcommands[length(addcommands)+1] <- paste("aictest=(",aictest,") savelog= aictest",sep="")
    if(!is.null(reguser)){
      forbidden.regression.user <- c("x11regress:",
          "samode:","finmode:","seasonalma:","trendma:","sfmsr:",
          "finalxreg","x11irrcrtval:",
          "$AO","User-defined$","Automatically Identified Outliers$",
          "peaks.seas:","peaks.td:","f2.idseasonal:","d11.f:",
          "spcori","spcsa","spcirr",
          "f3.m01:","f3.m02:","f3.m03:","f3.m04:","f3.m05:","f3.m06:",
          "f3.m07:","f3.m08:","f3.m09:","f3.m10:","f3.m11:",
          "f3.q:","f3.qm2:","f3.fail:",
          "ssa:","ssfstab:","ssfmov:","ssm7:","ssident:","ssran.","s2.","s3.",
          "historytarget","r01.lag","r02.lag","r04.lag","r05.lag","r06","meanssfe")
      if(!any(unlist(lapply(forbidden.regression.user,function(x)grepl(x,reguser))))){				  				  
        addcommands[length(addcommands)+1] <- paste("user=(",paste(reguser,collapse=" "),")",sep="")
      }else{
        bad.name.regression.user <- unlist(lapply(forbidden.regression.user,function(x)grep(x,reguser,value=TRUE)))
        bad.name.index <- which(reguser%in%bad.name.regression.user)
        reguser[bad.name.index]<-paste("user_",1:length(bad.name.regression.user),sep="")
        cat("Warning: the user paramter/s",bad.name.regression.user,"in the regression argument 'regression.user' has/have been renamed to",reguser[bad.name.index],"due to conflicts! \n")
        addcommands[length(addcommands)+1] <- paste("user=(",paste(reguser,collapse=" "),")",sep="")
      }
    }
    if(!is.null(regfile))
      addcommands[length(addcommands)+1] <- paste("file='",regfile,"'",sep="")
    if(!is.null(regfilestart))
      addcommands[length(addcommands)+1] <- paste("start=",paste(regfilestart,collapse="."),"",sep="")
    if(!is.null(usertype))
      addcommands[length(addcommands)+1] <- paste("usertype=(",paste(usertype,collapse=" "),")",sep="")
    if(!is.null(centeruser))
      addcommands[length(addcommands)+1] <- paste("centeruser=",centeruser,sep="")
    addcommands[length(addcommands)+1] <- "}"
  }
  if(!is.null(outlier) &&! x11regress){
#		  if((outlier.detection || !is.null(outlier)) &&! x11regress){			  
    addcommands[length(addcommands)+1] <- "outlier {"
#		  if(!is.null(outlier)){
#			  outlier.detection <- TRUE
    if(all(outlier=="all"))
      addcommands[length(addcommands)+1] <- "types=(all)"
    else
      addcommands[length(addcommands)+1] <- paste("types=(",paste(outlier,collapse=" "),")",sep="")
    
#	  		}
    if(!is.null(critical)){
      if(is.list(critical)){
        names(critical)<-toupper(names(critical))
        critval <- vector()
        ifelse(is.null(critical$AO),critval[1] <- "",critval[1] <- critical$AO)
        ifelse(is.null(critical$LS),critval[2] <- "",critval[2] <- critical$LS)
        ifelse(is.null(critical$TC),critval[3] <- "",critval[3] <- critical$TC)
        addcommands[length(addcommands)+1] <- paste("critical=(",paste(critval,collapse=","),")",sep="")
      }else{addcommands[length(addcommands)+1] <- paste("critical=(",paste(critical,collapse=","),")",sep="")	
      }
    }
    if(!is.null(outlier_span)){
      topaste<-outlier_span
      tocollapse<-c(".",".")
      if(any(is.na(outlier_span))){
        topaste[which(is.na(outlier_span))]<-""
        tocollapse[which(is.na(outlier_span))[2]/2]<-""
      }
      addcommands[length(addcommands)+1] <- paste("span=(",paste(topaste[1:2],collapse=tocollapse[1]),",",paste(topaste[3:4],collapse=tocollapse[2]),")",sep="")
    }
    
#		  if(!is.null(outlier_span))
#			  addcommands[length(addcommands)+1] <- paste("span=(",paste(outlier_span,collapse=","),")",sep="")	
    addcommands[length(addcommands)+1] <- "print=(default)"
    if(!is.null(outlier_method) &&! x11regress){
      addcommands[length(addcommands)+1] <- paste("method=",paste(outlier_method,collapse=","),sep="")	
    }
    addcommands[length(addcommands)+1] <- "}"
  }
  if(identify){
    addcommands[length(addcommands)+1] <- "identify {"
    
    if(!is.null(identify.diff))
      addcommands[length(addcommands)+1] <- paste("diff=",paste("(",paste(identify.diff,collapse=","),")",sep=""),sep="")	
    if(!is.null(identify.sdiff))
      addcommands[length(addcommands)+1] <- paste("sdiff=",paste("(",paste(identify.sdiff,collapse=","),")",sep=""),sep="")
    if(!is.null(identify.maxlag))
      addcommands[length(addcommands)+1] <- paste("maxlag=",identify.maxlag,sep="")	
    addcommands[length(addcommands)+1] <- "}"	
  }	
  if(slidingspans){
    addcommands[length(addcommands)+1] <- "slidingspans{" 
    if(!is.null(slidingspans.fixmdl))
      addcommands[length(addcommands)+1] <- paste("fixmdl=",slidingspans.fixmdl,sep="")	
    
    if(!is.null(slidingspans.fixreg))
      addcommands[length(addcommands)+1] <- paste("fixreg=(",paste(slidingspans.fixreg,collapse=" "),")",sep="")	
    
    if(!is.null(slidingspans.length))
      addcommands[length(addcommands)+1] <- paste("length=",slidingspans.length,sep="")	
    
    if(!is.null(slidingspans.numspans))
      addcommands[length(addcommands)+1] <- paste("numspans=",slidingspans.numspans,sep="")	
    
    if(!is.null(slidingspans.outlier))
      addcommands[length(addcommands)+1] <- paste("outlier=",slidingspans.outlier,sep="")	
    
    if(!is.null(slidingspans.start))
      addcommands[length(addcommands)+1] <-  paste("start=",paste(slidingspans.start,collapse="."),"",sep="")
    
    if(!is.null(slidingspans.additivesa))
      addcommands[length(addcommands)+1] <- paste("additivesa=",slidingspans.additivesa,sep="")	
    
    addcommands[length(addcommands)+1] <- "}" 
  }
  if(history){
    addcommands[length(addcommands)+1] <- "history{" 
    if(!is.null(history.estimates))
      addcommands[length(addcommands)+1] <- paste("estimates=(",paste(history.estimates,collapse=" "),")",sep="")	
    
    if(history.fixmdl)
      addcommands[length(addcommands)+1] <- "fixmdl=yes"
    
    if(!is.null(history.fixreg))
      addcommands[length(addcommands)+1] <- paste("fixreg=(",paste(history.fixreg,collapse=" "),")",sep="")
    
    if(!is.null(history.outlier))
      addcommands[length(addcommands)+1] <- paste("outlier=",history.outlier,sep="")	
    
    if(!is.null(history.sadjlags))
      addcommands[length(addcommands)+1] <- paste("sadjlags=",paste("(",paste(history.sadjlags,collapse=","),")",sep=""),sep="")	
    
    if(!is.null(history.trendlags))
      addcommands[length(addcommands)+1] <- paste("trendlags=",paste("(",paste(history.trendlags,collapse=","),")",sep=""),sep="")	
    
    if(!is.null(history.start))
      addcommands[length(addcommands)+1] <- paste("start=",paste(history.start,collapse="."),"",sep="")
    
    if(!is.null(history.target))
      addcommands[length(addcommands)+1] <- paste("target=",history.target,sep="")
    
    addcommands[length(addcommands)+1] <- "}" 
  }
  if(!x11regress){#nicht bei x11 Regression
    if(estimate){
      addcommands[length(addcommands)+1] <- "estimate {"
      if(estOutofsample){	
        addcommands[length(addcommands)+1] <- "outofsample=yes"}
      addcommands[length(addcommands)+1] <- "print=(default + rts)"
      addcommands[length(addcommands)+1] <- "savelog=(aic bic afc)"
      addcommands[length(addcommands)+1] <- "}"			  
      if(check){
        addcommands[length(addcommands)+1] <- "check{"
        if(!is.null(check.maxlag))
          addcommands[length(addcommands)+1] <- paste("maxlag=",check.maxlag,sep="")	
        #addcommands[length(addcommands)+1] <- "print=(default+specresidual+pacfplot)"
        addcommands[length(addcommands)+1] <- "savelog=(nrm lbq)"
        addcommands[length(addcommands)+1] <- "}"
      }
    }
    if(automdl && is.null(c(arima,sarima,arima.ar,arima.ma))){
      addcommands[length(addcommands)+1] <- "automdl{"
      if(acceptdefault)
        addcommands[length(addcommands)+1] <- "acceptdefault=yes"
      else
        addcommands[length(addcommands)+1] <- "acceptdefault=no"
      if(balanced)
        addcommands[length(addcommands)+1] <- "balanced=yes"
      else
        addcommands[length(addcommands)+1] <- "balanced=no"
      
      maxorder[is.na(maxorder)]<-" "
      addcommands[length(addcommands)+1] <- paste("maxorder=",paste("(",paste(maxorder,collapse=","),")",sep=""),sep="") 
      maxdiff[is.na(maxdiff)]<-" "
      addcommands[length(addcommands)+1] <- paste("maxdiff=",paste("(",paste(maxdiff,collapse=","),")",sep=""),sep="") 
      
      addcommands[length(addcommands)+1] <- "savelog=(adf amd b5m mu)"
      addcommands[length(addcommands)+1] <- "}" }
#Forecasts Backcasts
    if(!is.null(forecast_years) | !is.null(backcast_years)){
      addcommands[length(addcommands)+1] <- "forecast {"
      addcommands[length(addcommands)+1] <- "save=ftr"
      if(!is.null(forecast_years)){
        addcommands[length(addcommands)+1] <- paste("maxlead=",forecast_years*frequency(tso),sep="")
      }
      if(!is.null(backcast_years)){
        addcommands[length(addcommands)+1] <- paste("maxback=",backcast_years*frequency(tso),sep="")
      }
      addcommands[length(addcommands)+1] <- "}"
    }
  }#end nicht bei x11 Regression
  if(!seats){
    addcommands[length(addcommands)+1] <- "x11{"
    addcommands[length(addcommands)+1] <- "save=(d10 d11 d12)"
    if(!is.null(onlytd)){
      addcommands[length(addcommands)+1] <- paste("type=",onlytd,sep="")		
    }  			  
    if(sfshort)  
      addcommands[length(addcommands)+1] <- "sfshort=yes"
    if(!is.null(sigmalim)){
      sigmalim <- paste("(",sigmalim[1],",",sigmalim[2],")",sep="")
      addcommands[length(addcommands)+1] <- paste("sigmalim=",sigmalim,sep="")
    }
    if(!is.null(samode))
      addcommands[length(addcommands)+1] <- paste("mode=",samode,sep="")	
    if(!is.null(seasonalma)){
      addcommands[length(addcommands)+1] <- paste("seasonalma=(",paste(seasonalma,collapse=" "),")",sep="")}	
    if(!is.null(trendma)){
      addcommands[length(addcommands)+1] <- paste("trendma=",trendma,sep="")		
    }
    if(!is.null(x11calendarsigma))
      addcommands[length(addcommands)+1] <- paste("calendarsigma=",x11calendarsigma,sep="")
    if(x11excludefcst)
      addcommands[length(addcommands)+1] <- "excludefcst=yes"
    if(x11appendbcst)
      addcommands[length(addcommands)+1] <- "appendbcst=yes" ###backcast
    if(x11final!="none")
      addcommands[length(addcommands)+1] <- paste("final=(",paste(x11final,collapse=" "),")",sep="")
    if(x11appendfcst)
      addcommands[length(addcommands)+1] <- "appendfcst=yes" ###forecast		  
    addcommands[length(addcommands)+1] <- "savelog=all"
    addcommands[length(addcommands)+1] <- "}" 
  }else{
    addcommands[length(addcommands)+1] <- paste("seats{",seatsparameter,"}",sep="")
  }
  if(x11regress){
#start: The start date for the values of the user-defined regression variables.
# The default is the start date of the series. 
# Valid values are any date up to the start date of the series 
# (or up to the start date of the span specified by the span argument of the series spec, if present).
    addcommands[length(addcommands)+1] <- "x11regression{"
    if(!is.null(regfilestart))
      addcommands[length(addcommands)+1] <- paste("start=",paste(regfilestart,collapse="."),sep="")
    else
      addcommands[length(addcommands)+1] <- paste("start=",paste(start(tso),collapse="."),sep="")
    if(!is.null(critical)){
      if(is.list(critical) & length(critical)>1 &!"AO"%in%names(critical)){
        cat("X11 Regression only allows for the detection of Additive Outliers (AO)! \n")}
      else
        addcommands[length(addcommands)+1] <- paste("critical=",critical,sep="")	
    }		
    if(!is.null(outlier_method)){
      addcommands[length(addcommands)+1] <- paste("outliermethod=",paste(outlier_method,collapse=","),sep="")	
    }
    if(!is.null(outlier_span)){
      topaste<-outlier_span
      tocollapse<-c(".",".")
      if(any(is.na(outlier_span))){
        topaste[which(is.na(outlier_span))]<-""
        tocollapse[which(is.na(outlier_span))[2]/2]<-""
      }
      addcommands[length(addcommands)+1] <- paste("outlierspan=(",paste(topaste[1:2],collapse=tocollapse[1]),",",paste(topaste[3:4],collapse=tocollapse[2]),")",sep="")
    }		  
    if(!is.null(regvariables))
      addcommands[length(addcommands)+1] <- paste("variables=(",paste(regvariables,collapse=" "),")",sep="")
    if(!is.null(reguser)){
      forbidden.regression.user <- c("x11regress:",
          "samode:","finmode:","seasonalma:","trendma:","sfmsr:",
          "finalxreg","x11irrcrtval:",
          "$AO","User-defined$","Automatically Identified Outliers$",
          "peaks.seas:","peaks.td:","f2.idseasonal:","d11.f:",
          "spcori","spcsa","spcirr",
          "f3.m01:","f3.m02:","f3.m03:","f3.m04:","f3.m05:","f3.m06:",
          "f3.m07:","f3.m08:","f3.m09:","f3.m10:","f3.m11:",
          "f3.q:","f3.qm2:","f3.fail:",
          "ssa:","ssfstab:","ssfmov:","ssm7:","ssident:","ssran.","s2.","s3.",
          "historytarget","r01.lag","r02.lag","r04.lag","r05.lag","r06","meanssfe")
      if(!any(unlist(lapply(forbidden.regression.user,function(x)grepl(x,reguser))))){				  				  
        addcommands[length(addcommands)+1] <- paste("user=(",paste(reguser,collapse=" "),")",sep="")
      }else{
        bad.name.regression.user <- unlist(lapply(forbidden.regression.user,function(x)grep(x,reguser,value=TRUE)))
        bad.name.index <- which(reguser%in%bad.name.regression.user)
        reguser[bad.name.index]<-paste("user_",1:length(bad.name.regression.user),sep="")
        cat("Warning: the user paramter/s",bad.name.regression.user,"in the regression argument 'regression.user' has/have been renamed to",reguser[bad.name.index],"due to conflicts! \n")
        addcommands[length(addcommands)+1] <- paste("user=(",paste(reguser,collapse=" "),")",sep="")
      }
    }
    if(!is.null(regfile))
      addcommands[length(addcommands)+1] <- paste("file='",regfile,"'",sep="")				
    if(!is.null(centeruser))
      addcommands[length(addcommands)+1] <- paste("centeruser=",centeruser,sep="")
    if(!is.null(usertype))
      addcommands[length(addcommands)+1] <- paste("usertype=(",paste(usertype,collapse=" "),")",sep="")
    addcommands[length(addcommands)+1] <- "}"	
    
  }
  
  
  con <- file(paste(file,".spc",sep=""))
  
  #writeLines(c(header,datarows,addcommands),con)
  writeLines(c(header,addcommands),con)
  close(con)
  
  # Rewritten, for not using sh or bat files (Suggestions by Peter Ellis)
  if(Sys.info()[1]=="Windows"){
    #con1 <- file("run.bat")
    #mdcommand <- "md gra"
    file_1 <- gsub("/","\\\\",file)
    if((!is.null(x12path)) && use=="x12"){
      x12path_1 <- gsub("/","\\\\",x12path)
      command <- paste(x12path_1," ",file_1," -g ",dirgra,sep="")
    }else if((!is.null(x13path)) && use!="x12"){
      x13path_1 <- gsub("/","\\\\",x13path)
      command <- paste(x13path_1," ",file_1," -g ",dirgra,sep="")
    }else 
      stop("Please define the path to the X12 binaries!")
  }else{
    #con1 <- file("run.sh")
    #mdcommand <- "mkdir gra"
    if((!is.null(x12path)) && use=="x12"){
      command <- paste(x12path," ",file," -g ",dirgra,sep="")
    }else if((!is.null(x13path)) && use!="x12"){
      command <- paste(x13path," ",file," -g ",dirgra,sep="")
    }else
      stop("Please define the path to the X12 binaries!")
  }
  #writeLines(c(mdcommand,command),con1)
  #close(con1)
  if(!file.exists(dirgra))
    dir.create(dirgra) 
  system(command) 
#  if(Sys.info()[1]=="Windows"){
#    system("run.bat")
#  }else{
#    system("chmod 744 run.sh")
#    system("./run.sh")
#  }
  
#  out <- list()
  start_series <- start(tso)
  end_series <- end(tso)
  if(!is.null(span)){
    if(!any(is.na(span[1:2])))
      start_series <- as.numeric(span[1:2])
    if(!any(is.na(span[3:4])))
      end_series <- as.numeric(span[3:4])
  
  }
  out <- readx12Out(file,freq_series=frequency(tso),start_series=start_series,end_series=end_series,tblnames=tblnames,Rtblnames=Rtblnames,transform=transform,slidingspans=slidingspans,history=history,x11regress=x11regress,outlier=outlier,showWarnings=showWarnings,keep_x12out=keep_x12out)
#  Rtblnames <- c("Original series", "Final seasonal factors", "Final seasonally adjusted data", "Final trend cycle",
#		    "Final irregular components","Combined adjustment factors","Final weights for irregular component",
#			"Final replacements for SI ratios",
#			"Differenced, transformed, seasonally adjusted data",
#			Rtblnames)
#  if(seats==TRUE)
#    tblnames <- c("a1", "s10", "s11", "s12", "s13","s16","c17","s9","e2", tblnames)
#  else
#    tblnames <- c("a1", "d10", "d11", "d12", "d13","d16","c17","d9","e2", tblnames)
#  for(i in 1:length(tblnames)){
#    if(file.exists(paste("gra\\",file,".",tblnames[i],sep="")))
#      out[[tblnames[i]]] <- ts(read.table(paste("gra\\",file,".",tblnames[i],sep=""),header=FALSE,skip=2,sep="	",na.strings="-999")[,2],frequency=frequency(tso),start=start(tso))
#  }
#  spnames <- c("Spectrum_AdjOri","Spectrum_SA","Spectrum_Irr")
#  sptblnames <- c("sp0", "sp1", "sp2")
#  if(!seats){
#    for(i in 1:length(sptblnames)){
#      out[[spnames[i]]] <- read.table(paste("gra\\",file,".",sptblnames[i],sep=""),header=FALSE,skip=2,sep="	")[,2:3]
#      names(out[[spnames[i]]]) <- c("frequency","spectrum")
#    }
#  }
#  out[["d9"]][out[["d9"]]==-999]<-NA
#  out[["Forecast with CI"]] <- list()
#  fct <- read.table(paste("gra\\",file,".","fct",sep=""),header=FALSE,skip=2,sep="	")
#  out[["Forecast with CI"]]$estimate <-ts(fct[,2],frequency=frequency(tso),start=end(tso)) 
#  out[["Forecast with CI"]]$lower <-ts(fct[,3],frequency=frequency(tso),start=end(tso))
#  out[["Forecast with CI"]]$upper <-ts(fct[,4],frequency=frequency(tso),start=end(tso))
#  out$seats <- seats
#  out$file <- file
#  out$tblnames <- tblnames
#  out$Rtblnames <- Rtblnames
#  class(out) <- "x12"
  
  ext <- c("out","err","spc","otl","dat","a1","b1","d10","d11","d12","ftr","log")
  if(!keep_x12out){
    for(e in ext){
      f <- paste(basename(file),".",e,sep="")
      if(file.exists(f))
        file.remove(f)
    }
  }
  
  if(!keep_x12out)
    unlink(paste(dirname(file),"/",dirgra,sep=""),recursive=TRUE)
  out
}




#print.x12work <- function(x,editor=getOption("editor"),...){
#  if(!(x$file=="Example_for_X12"))
#    filename <- paste(x$file,".out",sep="")
#  else
#    filename <- paste(paste(searchpaths()[grep("x12",searchpaths())],"/doc/Rout",sep=""),".out",sep="")
#  edit(file=filename,editor=editor,...)
#}

