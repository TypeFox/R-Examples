readx12Out <- function(file,tblnames=NULL,Rtblnames=NULL,freq_series,start_series,end_series,
    seats=FALSE,transform,slidingspans,history,x11regress,outlier,showWarnings,keep_x12out){
#  cat("start_series: ")
#  print(start_series)
  dirgra <- paste("gra_",gsub("\\.","_",basename(file)),sep="")
  out<-list()
  Rtblnames <- c("Original series", "Final seasonal factors", "Final seasonally adjusted data", "Final trend cycle",
      "Final irregular components","Combined adjustment factors","Final weights for irregular component",
      "Final replacements for SI ratios",
      "Differenced, transformed, seasonally adjusted data","Final unmodified SI Ratios","Orig2","Trading day component",
      Rtblnames)
  if(seats==TRUE)
    tblnames <- c("a1", "s10", "s11", "s12", "s13","s16","c17","s9","e2","d8","b1","td", tblnames)
  else
    tblnames <- c("a1", "d10", "d11", "d12", "d13","d16","c17","d9","e2","d8","b1","td", tblnames)
  if(!(file=="Example_for_X1")){
    sp_file <- strsplit(file,"/")[[1]]
    filename <- paste(paste(sp_file[-length(sp_file)],collapse="/"),"/",dirgra,"/",sp_file[length(sp_file)],sep="")
    if(substring(filename,1,1)=="/")
      filename <- substring(filename,2)
  }else
    filename <- paste(searchpaths()[grep("x12",searchpaths())],"/doc/Rout",sep="")
  ############  
  
#if(file.exists(paste(file,".","err",sep="")) && file.exists(paste(file,".","spc",sep=""))){#evt unnoetige Abfrage
  
#	x1<-	file.info(paste(filename,".","udg",sep=""))
#	cat("udg:",str(x1),"\n")
#	x2<-	file.info(paste(file,".","err",sep=""))
#		cat("err:",str(x2),"\n")
#	ind<-	!identical(x1$mtime,x2$mtime)
#	cat("ind:",ind,"\n")	
#		udgtime<-file.info(paste(filename,".","udg",sep=""))$mtime
#		errtime<-file.info(paste(file,".","err",sep=""))$mtime
  
#		if(ind){
#			#cat("udg:",file.info(paste(filename,".","udg",sep=""))$mtime,"\n")
#			#cat("err:",file.info(paste(file,".","err",sep=""))$mtime,"\n")
#			errorfile <- readLines(con=paste(file,".","err",sep=""),n=-1)
#			for(i in 1:length(errorfile)){
#				cat(errorfile[i],"\n")	
#			}
#			stop("Error! No proper run of x12! Check your parameter settings.\n=> Beware that files in \"gra\" directory do not represent current x12 output")
#		}
  
  if(!file.exists(paste(filename,".","udg",sep=""))){
    
    errorfile <- readLines(con=paste(file,".","err",sep=""),n=-1)
    for(i in 1:length(errorfile)){
      cat(errorfile[i],"\n")	
    }
    if(!keep_x12out)
      unlink(paste(dirname(file),"/",dirgra,sep=""),recursive=TRUE)
    
    stop("Error! No proper run of x12! Check your parameter settings.")	
  }
  udg <- readLines(con=paste(filename,".","udg",sep=""),n=-1)
  
  
  if(showWarnings && file.exists(paste(file,".","err",sep=""))){
    errorfile <- readLines(con=paste(file,".","err",sep=""),n=-1)
    for(i in min(which(errorfile=="")):length(errorfile)){
      cat(errorfile[i],"\n")	
    }
  }
  if(any(any(grepl("errorstop: yes",udg)),length(udg)==0)){
    if(!showWarnings){
      errorfile <- readLines(con=paste(file,".","err",sep=""),n=-1)
      for(i in which(grepl("ERROR:",errorfile)):length(errorfile)){
        cat(errorfile[i],"\n")	
      }
      if(!keep_x12out)
        unlink(paste(dirname(file),"/",dirgra,sep=""),recursive=TRUE)
      
      stop("An error occured when running x12! Program halted!","\n") 	
    }else{
      if(!keep_x12out)
        unlink(paste(dirname(file),"/",dirgra,sep=""),recursive=TRUE)
      
      stop("An error occured when running x12! Program halted!","\n") 	
    }
#stop("Errorstop! (Check error file \"",file,".err\")",sep="")	
  }
  
  ############	
  for(i in 1:length(tblnames)){
    if(file.exists(paste(filename,".",tblnames[i],sep="")))
      out[[tblnames[i]]] <- ts(read.table(paste(filename,".",tblnames[i],sep=""),header=FALSE,skip=2,sep="	",na.strings="-999")[,2],frequency=freq_series,start=start_series)
  }
#  if(!x11regress){
#  spnames <- c("Spectrum_AdjOri","Spectrum_SA","Spectrum_Irr","Spectrum_Rsd")
#  sptblnames <- c("sp0", "sp1", "sp2","spr")
#	}else{
#	spnames <- c("Spectrum_AdjOri","Spectrum_SA","Spectrum_Irr")
#	sptblnames <- c("sp0", "sp1", "sp2")	
#	}
  
  sptblnames <- vector()
  for(i in c("sp0", "sp1", "sp2","spr")){
    if(file.exists(paste(filename,".",i,sep="")))
      sptblnames <- c(sptblnames,i)	
  }
  if(!seats){
    for(i in 1:length(sptblnames)){
      out[[sptblnames[i]]] <- read.table(paste(filename,".",sptblnames[i],sep=""),header=FALSE,skip=2,sep="	")[,2:3]
      names(out[[sptblnames[i]]]) <- c("frequency","spectrum")
    }
  }
  out[["d9"]][out[["d9"]]==-999]<-NA
  if(!x11regress){
    ### Forecasts:	  
    if(file.exists(paste(filename,".","fct",sep=""))){
      out[["forecast"]] <- list()
      fct <- read.table(paste(filename,".","fct",sep=""),header=FALSE,skip=2,sep="	")
      
      if((freq_series==12 && end_series[2]==12) | (freq_series==4 && end_series[2]==4)){
        start_forecast <- c(end_series[1]+1,1)
#  }else if(freq_series==4 && end_series[2]==4){
#	  start_forecast <- c(end_series[1]+1,1)}
      }else{
        start_forecast <- end_series
        start_forecast[2] <- start_forecast[2]+1 
      }
      out[["forecast"]]$estimate <-ts(fct[,2],frequency=freq_series,start=start_forecast) 
      out[["forecast"]]$lowerci <-ts(fct[,3],frequency=freq_series,start=start_forecast)
      out[["forecast"]]$upperci <-ts(fct[,4],frequency=freq_series,start=start_forecast)
    }
    ### Backcasts:
    if(file.exists(paste(filename,".","bct",sep=""))){
      out[["backcast"]] <- list()
      bct <- read.table(paste(filename,".","bct",sep=""),header=FALSE,skip=2,sep="	")
      if(start_series[2]==1){
        if(freq_series==12)
          #start_backcast <- c(start_series[1]-1,12)
          end_backcast <- c(start_series[1]-1,12)	
        if(freq_series==4)
          end_backcast <- c(start_series[1]-1,4)	
      }else{
        end_backcast <- start_series
        end_backcast[2] <- end_backcast[2]-1 
      }
      out[["backcast"]]$estimate <-ts(bct[,2],frequency=freq_series,end=end_backcast) 
      out[["backcast"]]$lowerci <-ts(bct[,3],frequency=freq_series,end=end_backcast)
      out[["backcast"]]$upperci <-ts(bct[,4],frequency=freq_series,end=end_backcast)
    }
    
  }
#Testbeispiele:
#filename <- "M:/Meraner/Workspace/Saisonbereinigung_Test/gra/Air"
#filename <- "M:/Meraner/Saisonbereinigung/x12probierfiles/loge_d"  
#filename <- "M:/Meraner/Saisonbereinigung/x12probierfiles/a"  
#filename <- "M:/Meraner/Saisonbereinigung/x12probierfiles/b05"
#file<-"M:/Meraner/Workspace/Saisonbereinigung_Test/Air"
#filename <- "M:/Meraner/Workspace/Saisonbereinigung_Test/gra/Rout"
  
#udg <- readLines(con=paste(filename,".","udg",sep=""),n=-1)
  
  
  if("x11regress: no" %in% udg){
#Informationen die aus udg File eingelesen werden sollen:
    dglist <- c("x11regress:","transform:","samode:","finmode:","seasonalma:","trendma:","sfmsr:",
        "arimamdl:","automdl:", 
        "finalreg",
        "outlier.total:","autoout:","nalmostout:",
        "almostoutlier$","crit:",
        "Outlier$","User-defined$","AutoOutlier$",
        #Autokorrelationen
        #"sfmsr:","finaltrendma:",
        "peaks.seas:","peaks.td:",  
        "f2.idseasonal:",
        "d11.f:",
        "spcrsd",
        "spcori",
        "spcsa",
        "spcirr",
        "f3.m01:","f3.m02:","f3.m03:","f3.m04:","f3.m05:","f3.m06:",
        "f3.m07:","f3.m08:","f3.m09:","f3.m10:","f3.m11:",
        "f3.q:","f3.qm2:","f3.fail:",
        "loglikelihood:","aic:","aicc:","bic:","hq:","aape")
    
#Hier die gewuenschten Variablennamen eintragen:
    dglistnames <-  c("x11regress","transform","samode","finalsamode","seasonalma","trendma","finalseasonalma",
        "arimamdl","automdl", 
        "regmdl",
        "nout","nautoout","nalmostout",
        "almostoutlier","crit",
        "outlier","userdefined","autooutlier",
        #Autokorrelationen
        #"sfmsr","finaltrendma",
        "peaks.seas","peaks.td",  
        "id.seas",
        "id.rsdseas",
        "spcrsd",
        "spcori",
        "spcsa",
        "spcirr",
        "m1","m2","m3","m4","m5","m6","m7","m8",
        "m9","m10","m11",
        "q","q2","nmfail",
        "loglikelihood","aic","aicc","bic","hq","aape") 
    
    if(transform=="auto"){# &&!x11regress --- x11regress Abfrage hier eigentl nicht mehr notw
      dglist[length(dglist)+1]<-"aictrans:"
      dglistnames[length(dglistnames)+1]<-"autotransform"
    }
    
    numvariables <- c("nout","nautoout","nalmostout",
        "crit","spcrsd",
        "spcori",
        "spcsa",
        "spcirr",
        "m1","m2","m3","m4","m5","m6","m7","m8",
        "m9","m10","m11",
        "q","q2","nmfail","seasonalma","trendma","finalseasonalma",
        "loglikelihood","aic","aicc","bic","hq","aape")
    
    if(slidingspans){
      dglist[(length(dglist)+1):(length(dglist)+8)]<-c("ssa:","ssfstab:","ssfmov:","ssm7:","ssident:",
          "ssran.","s2.","s3.")
      
      dglistnames[(length(dglistnames)+1):(length(dglistnames)+8)]<-c("ss.options","ss.stabseas","ss.movseas","ss.m7","ss.idseas",
          "ss.S1","ss.S2","ss.S3")
    }
    if(history){
      dglist[(length(dglist)+1):(length(dglist)+7)]<-c("historytarget",
          "r01.lag","r02.lag","r04.lag","r05.lag","r06","meanssfe")
      
      dglistnames[(length(dglistnames)+1):(length(dglistnames)+7)]<-c("h.target",
          "h.R1","h.R2","h.R4","h.R5","h.R6","h.meanssfe")
    }
    
    
    dg <- lapply(dglist,function(x)grep(x,udg,value=TRUE,fixed=TRUE))
    
#Extrawurst fuer Regression variables Teil 1:
    if(length(dg[[which(dglist=="finalreg")]])>1){
      dg[[which(dglist=="finalreg")]]<-dg[[which(dglist=="finalreg")]][-grep("nfinalreg",dg[[which(dglist=="finalreg")]])]	
    }
    
    regvar<-unlist(strsplit(dg[[which(dglist=="finalreg")]],": "))
    if(any(grepl("+",regvar,fixed=TRUE))){
      regvar <- strsplit(regvar,"+",fixed=TRUE)
      regvar<-unlist(lapply(regvar,function(x)gsub("^\\s+|\\s+$", "", x))[-(grep("finalreg",regvar))])
      if(""%in%regvar){
        regvar<-regvar[-which(regvar=="")]}  
    }
    othername<-which(regvar %in% c("User-defined","Automatically Identified Outliers"))
    if(length(othername)>0){
      regvar <- regvar[-othername]}
#End Extrawurst Teil 1
    empty <- which(lapply(dg,function(x)length(x))==0)  
    dglist[which(dglist%in%grep("$",grep(":",dglist[empty],invert=TRUE,fixed=TRUE,value=TRUE),invert=TRUE,fixed=TRUE,value=TRUE))]<-paste(grep("$",grep(":",dglist[empty],invert=TRUE,fixed=TRUE,value=TRUE),invert=TRUE,fixed=TRUE,value=TRUE),":",sep="")  
    dg[empty] <- paste(gsub("$",replacement=":",dglist[empty],fixed=TRUE),"-")
    dg <- lapply(dg,function(x)strsplit(x,": "))
    names(dg)<- dglistnames
    
    if(length(which(dg[["outlier"]]%in%dg[["autooutlier"]]))>0){
      dg[["outlier"]]<-dg[["outlier"]][-which(dg[["outlier"]]%in%dg[["autooutlier"]])]}
    
    if(length(dg[["outlier"]])==0){
      dg["outlier"] <- list(strsplit(paste("outlier","-")," ",fixed=TRUE))
      empty <- c(empty,which(names(dg)=="outlier"))
    }
    
    grone <- which(lapply(1:length(dg),function(x){length(dg[[x]])})>1)
    grone.nodoll <- unlist(lapply(grep("$",dglist[grone],fixed=TRUE,value=TRUE,invert=TRUE),function(x){
              grep(x,dglist,fixed=TRUE)}))	
    for(i in c(1:length(dg))[-grone.nodoll]){
#	  cat("i=",i,"\n")
#	  cat("name=",names(dg)[i],"\n")
      names(dg[[i]])<-lapply(1:length(dg[[i]]),function(x){
            if(grepl("$",dg[[i]][[x]][1],fixed=TRUE)){
              gsub(grep("$",dglist[i],fixed=TRUE,value=TRUE),replacement=paste(dglistnames[i],"_",sep=""),dg[[i]][[x]][1],fixed=TRUE)
            }else{
              gsub(dg[[i]][[x]][1],replacement=dglistnames[i],dg[[i]][[x]][1],fixed=TRUE)
            }})
      for(j in 1:length(dg[[i]])){
        dg[[i]][[j]] <- dg[[i]][[j]][-1]
        # returns string w/o leading or trailing whitespace:
        dg[[i]][[j]] <- gsub("^\\s+|\\s+$", "", dg[[i]][[j]])
      }}
    
#Extrawurst fuer Regression variables Teil 2:
    if(length(regvar)!=0 &&regvar!="none"){  
      reglist <- lapply(1:length(regvar),function(x)grep(paste(regvar[x],"$",sep=""),grep(regvar[x],udg,value=TRUE,fixed=TRUE),fixed=TRUE,value=TRUE))
      reglist <- lapply(reglist,function(x)strsplit(x,"$",fixed=TRUE))
      if(length(which(sapply(reglist,function(x)!length(x)>0)))!=0){
        reglist <- reglist[-(which(sapply(reglist,function(x)!length(x)>0)))]}
      if(length(reglist)!=0){ 
        reglistnames <- vector()
        for(i in 1:length(reglist)){
          regsublistnames<-vector()
          inregvar<-which(sapply(1:length(reglist[[i]]),function(x)!reglist[[i]][[x]][1]%in%regvar))
          if(length(inregvar)>0){
            reglist[[i]] <- reglist[[i]][-(inregvar)]
          }
          reglistnames[i] <- reglist[[i]][[1]][1]
          reglist[[i]]<-lapply(1:length(reglist[[i]]),function(x){
                reglist[[i]][[x]]<-reglist[[i]][[x]][-1]
                strsplit(reglist[[i]][[x]],": ")})
          regsublistnames <-lapply(1:length(reglist[[i]]),function(x){
                regsublistnames <- reglist[[i]][[x]][[1]][1]
                if("Leap Year"%in%regsublistnames){
                  regsublistnames<-gsub("Leap Year",replacement="leapyear",regsublistnames,fixed=TRUE)  
                }
                if("Trading Day"%in%regsublistnames){
                  regsublistnames<-gsub("Trading Day",replacement="td",regsublistnames,fixed=TRUE)  
                }
                regsublistnames<-regsublistnames
              })
          #
          reglist[[i]]<-lapply(1:length(reglist[[i]]),function(x){reglist[[i]][[x]][[1]]<-reglist[[i]][[x]][[1]][-1]
                reglist[[i]][[x]][[1]] <- gsub("^\\s+|\\s+$", "", reglist[[i]][[x]][[1]])})
          if("Leap Year"%in%reglistnames){
            
            reglistnames<-gsub("Leap Year",replacement="leapyear",reglistnames,fixed=TRUE) 
          }
          if("Trading Day"%in%reglistnames){
            
            reglistnames<-gsub("Trading Day",replacement="td",reglistnames,fixed=TRUE) 
          }
          names(reglist[[i]])<-regsublistnames
        }
        dg[(length(dg)+1):(length(dg)+length(reglist))]<-reglist
        names(dg)<-dglistnames<-c(dglistnames,reglistnames)
      }else{
        names(dg)<- dglistnames
        reglistnames <- NULL
      }}else{
      names(dg)<- dglistnames
      reglistnames <- NULL
    }
#End Extrawurst Teil 2 
#Checking for derived parameter estimates:
    if(any(grepl("nregderived:",udg,fixed=TRUE))){
      regderived<-grep("nregderived:",udg,value=TRUE,fixed=TRUE)	
      indnrd<-which(udg==regderived)
      nrd<-as.numeric(strsplit(regderived,": ",fixed=TRUE)[[1]][2])
      regder<-udg[(indnrd+1):(indnrd+nrd)]
      derived.coef<-sapply(strsplit(gsub("$",replacement="_",regder,fixed=TRUE),": ",fixed=TRUE),function(x)
            gsub("\\s+",replacement="",x[[1]]))
    }
  }
#End RegARIMA Option
  
# x11Regression Option:
  else{
#Informationen die aus udg File eingelesen werden sollen:
    dglist <- c("x11regress:",
        "samode:","finmode:","seasonalma:","trendma:","sfmsr:",
        "finalxreg","x11irrcrtval:",
        #"$",
        "$AO","User-defined$","Automatically Identified Outliers$",
        #"sfmsr:","finaltrendma:",
        "peaks.seas:","peaks.td:",  
        "f2.idseasonal:","d11.f:",
        "spcori",
        "spcsa",
        "spcirr",
        "f3.m01:","f3.m02:","f3.m03:","f3.m04:","f3.m05:","f3.m06:",
        "f3.m07:","f3.m08:","f3.m09:","f3.m10:","f3.m11:",
        "f3.q:","f3.qm2:","f3.fail:")
    
#Neue Variablennamen:
    dglistnames <-  c("x11regress",
        "samode","finalsamode","seasonalma","trendma","finalseasonalma",
        "regmdl","crit","outlier","userdefined","autooutlier",
        #Autokorrelationen
        #"sfmsr","finaltrendma",
        "peaks.seas","peaks.td",  
        "id.seas","id.rsdseas",
        "spcori",
        "spcsa",
        "spcirr",
        "m1","m2","m3","m4","m5","m6","m7","m8",
        "m9","m10","m11",
        "q","q2","nmfail") 
    numvariables <- c("crit",
        "spcori",
        "spcsa",
        "spcirr",
        "m1","m2","m3","m4","m5","m6","m7","m8",
        "m9","m10","m11",
        "q","q2","nmfail","seasonalma","trendma","finalseasonalma")
    
    if(slidingspans){
      dglist[(length(dglist)+1):(length(dglist)+8)]<-c("ssa:","ssfstab:","ssfmov:","ssm7:","ssident:",
          "ssran.","s2.","s3.")
      
      dglistnames[(length(dglistnames)+1):(length(dglistnames)+8)]<-c("ss.options","ss.stabseas","ss.movseas","ss.m7","ss.idseas",
          "ss.S1","ss.S2","ss.S3")
    }
    if(history){
      dglist[(length(dglist)+1):(length(dglist)+7)]<-c("historytarget",
          "r01.lag","r02.lag","r04.lag","r05.lag","r06","meanssfe")
      
      dglistnames[(length(dglistnames)+1):(length(dglistnames)+7)]<-c("h.target",
          "h.R1","h.R2","h.R4","h.R5","h.R6","h.meanssfe")
    }
    
    dg <- lapply(dglist,function(x)grep(x,udg,value=TRUE,fixed=TRUE))
#Extrawurst fuer Regression variables Teil 1:
    if(length(dg[[which(dglist=="finalxreg")]])>1){
      dg[[which(dglist=="finalxreg")]]<-dg[[which(dglist=="finalxreg")]][-grep("nfinalxreg",dg[[which(dglist=="finalxreg")]])]	
    }
    regvar<-unlist(strsplit(dg[[which(dglist=="finalxreg")]],": "))
    if(any(grepl("+",regvar,fixed=TRUE))){
      regvar <- strsplit(regvar,"+",fixed=TRUE)
      regvar<-unlist(lapply(regvar,function(x)gsub("^\\s+|\\s+$", "", x))[-(grep("finalxreg",regvar))])
      if(""%in%regvar){
        regvar<-regvar[-which(regvar=="")]}
    }
    othername<-which(regvar %in% c("User-defined","Automatically Identified Outliers"))
    if(length(othername)>0){
      regvar <- regvar[-othername]}
#End Extrawurst Teil 1
    
    empty <- which(lapply(dg,function(x)length(x))==0)
    dglist[which(dglist%in%grep("$",grep(":",dglist[empty],invert=TRUE,fixed=TRUE,value=TRUE),invert=TRUE,fixed=TRUE,value=TRUE))]<-paste(grep("$",grep(":",dglist[empty],invert=TRUE,fixed=TRUE,value=TRUE),invert=TRUE,fixed=TRUE,value=TRUE),":",sep="")  
    dg[empty] <- paste(gsub("$",replacement=":",dglist[empty],fixed=TRUE),"-")
    dg <- lapply(dg,function(x)strsplit(x,": "))
    names(dg)<- dglistnames
    
    if(length(which(dg[["outlier"]]%in%dg[["autooutlier"]]))>0){
      dg[["outlier"]]<-dg[["outlier"]][-which(dg[["outlier"]]%in%dg[["autooutlier"]])]}
    if(length(dg[["outlier"]])==0){
      dg["outlier"] <- list(strsplit(paste("outlier","-")," ",fixed=TRUE))
      empty <- c(empty,which(names(dg)=="outlier"))
    }
    
    grone <- which(lapply(1:length(dg),function(x){length(dg[[x]])})>1)
    grone.nodoll <- unlist(lapply(grep("$",dglist[grone],fixed=TRUE,value=TRUE,invert=TRUE),function(x){
              grep(x,dglist,fixed=TRUE)}))	
    for(i in c(1:length(dg))[-grone.nodoll]){
      names(dg[[i]])<-lapply(1:length(dg[[i]]),function(x){
            if(grepl("$",dg[[i]][[x]][1],fixed=TRUE)){			
              strsplit(dg[[i]][[x]][1],"$",fixed=TRUE)
              paste(dglistnames[i],"_",strsplit(dg[[i]][[x]][1],"$",fixed=TRUE)[[1]][2],sep="")
              #gsub(grep("$",dglist[i],fixed=TRUE,value=TRUE),replacement=paste(dglistnames[i],"_",sep=""),dg[[i]][[x]][1],fixed=TRUE)
            }else{
              gsub(dg[[i]][[x]][1],replacement=dglistnames[i],dg[[i]][[x]][1],fixed=TRUE)
            }})
      for(j in 1:length(dg[[i]])){
        dg[[i]][[j]] <- dg[[i]][[j]][-1]
        # returns string w/o leading or trailing whitespace:
        dg[[i]][[j]] <- gsub("^\\s+|\\s+$", "", dg[[i]][[j]])
      }}
    
#Extrawurst fuer Regression variables Teil 2:
    if(length(regvar)!=0 &&regvar!="none"){  
      reglist <- lapply(1:length(regvar),function(x)grep(paste(regvar[x],"$",sep=""),grep(regvar[x],udg,value=TRUE),fixed=TRUE,value=TRUE))
      reglist <- lapply(reglist,function(x)strsplit(x,"$",fixed=TRUE))
      if(length(which(sapply(reglist,function(x)!length(x)>0)))!=0){
        reglist <- reglist[-(which(sapply(reglist,function(x)!length(x)>0)))]}}
    
    if(length(reglist)!=0){
      reglistnames <- vector()
      for(i in 1:length(reglist)){
        regsublistnames<-vector()
        inregvar<-which(sapply(1:length(reglist[[i]]),function(x)!reglist[[i]][[x]][1]%in%regvar))
        if(length(inregvar)!=0){
          reglist[[i]] <- reglist[[i]][-(inregvar)]
        }
        reglistnames[i] <- reglist[[i]][[1]][1]
        reglist[[i]]<-lapply(1:length(reglist[[i]]),function(x){
              reglist[[i]][[x]]<-reglist[[i]][[x]][-1]
              strsplit(reglist[[i]][[x]],": ")})
        regsublistnames <-lapply(1:length(reglist[[i]]),function(x){
              regsublistnames <- reglist[[i]][[x]][[1]][1]
              if(any(grepl("AO",regsublistnames))){
                regsublistnames<-paste("outlier_",grep("AO",regsublistnames,value=TRUE),sep="")}
              if("Leap Year"%in%regsublistnames){
                regsublistnames<-gsub("Leap Year",replacement="leapyear",regsublistnames,fixed=TRUE)  
              }
              if("Trading Day"%in%regsublistnames){
                regsublistnames<-gsub("Trading Day",replacement="td",regsublistnames,fixed=TRUE)  
              }
              regsublistnames<-regsublistnames
            })
        
        reglist[[i]]<-lapply(1:length(reglist[[i]]),function(x){reglist[[i]][[x]][[1]]<-reglist[[i]][[x]][[1]][-1]
              reglist[[i]][[x]][[1]] <- gsub("^\\s+|\\s+$", "", reglist[[i]][[x]][[1]])})
        
        if("Leap Year"%in%reglistnames){ 
          reglistnames<-gsub("Leap Year",replacement="leapyear",reglistnames,fixed=TRUE) 
        }
        if("Trading Day"%in%reglistnames){
          reglistnames<-gsub("Trading Day",replacement="td",reglistnames,fixed=TRUE) 
        }
        if(any(grepl("AO",reglistnames))){
          reglistnames[i]<-"outlier"}	
        
        names(reglist[[i]])<-regsublistnames
      }
      dg[(length(dg)+1):(length(dg)+length(reglist))]<-reglist
      names(dg)<-dglistnames<-c(dglistnames,reglistnames)
    }else{
      names(dg)<- dglistnames
      reglistnames <- NULL
    }
#End Extrawurst Teil 2  
    
#Checking for derived parameter estimates:
    if(any(grepl("nxregderived:",udg,fixed=TRUE))){
      regderived<-grep("nxregderived:",udg,value=TRUE,fixed=TRUE)	
      indnrd<-which(udg==regderived)
      nrd<-as.numeric(strsplit(regderived,": ",fixed=TRUE)[[1]][2])
      regder<-udg[(indnrd+1):(indnrd+nrd)]
      derived.coef<-sapply(strsplit(gsub("$",replacement="_",regder,fixed=TRUE),": ",fixed=TRUE),function(x)
            gsub("\\s+",replacement="",x[[1]]))
    }
    
  }
#End x11Regression Option
  
#Extrawurst fuer Critical Value und die Spektren
  for(i in grone.nodoll){
    names(dg[[i]])<-lapply(1:length(dg[[i]]),function(x){names <- dg[[i]][[x]][1]})
    for(j in 1:length(dg[[i]])){
      dg[[i]][[j]] <- dg[[i]][[j]][-1]
      dg[[i]][[j]] <- gsub("^\\s+|\\s+$", "", dg[[i]][[j]])
    }
  }
  
#Numerische Daten auch als solche ausgeben lassen:
  suppressWarnings(for(i in which(names(dg) %in% numvariables)){
        for(j in 1:length(dg[[i]])){
          num.possible <- as.numeric(dg[[i]][[j]])
          if(length(which(is.na(num.possible)))==0){
            dg[[i]][[j]] <- num.possible	
          }else{	
            num.split <- strsplit(dg[[i]][[j]],"\\s+")
            if(class(num.split)=="list"){
              ex.num <- as.numeric(unlist(num.split))
              num.split <- as.list(num.split[[1]])
            }else{
              ex.num <- as.numeric(num.split)
            }
            if(length(which(!is.na(ex.num)))>0){
              dg[[i]][[j]]<- as.list(num.split)
              dg[[i]][[j]][[which(!is.na(ex.num))]]<- ex.num[which(!is.na(ex.num))]
            }	
          }
        }})	
  
  sel <- c("userdefined","autooutlier","outlier",reglistnames) #werden extra behandelt
  
  for(i in which(names(dg) %in% sel)[!which(names(dg) %in% sel)%in%empty]){
    for(j in 1:length(dg[[i]])){
#(Regular Expressions as used in R)
#+ :The preceding item will be matched one or more times.
#Symbols \d, \s, \D and \S denote the digit and space classes and their negations. 
      dg[[i]][[j]] <- strsplit(dg[[i]][[j]],"\\s+")
      dg[[i]][[j]][[1]]<-as.numeric(dg[[i]][[j]][[1]])
      names(dg[[i]][[j]][[1]])<- c("coef","stderr","tval")#"coef" oder "estimate"
    }}	
  
  if("almostoutlier"%in%names(dg)){
    if(!which(names(dg)=="almostoutlier")%in%empty){
      for(j in 1:length(dg[["almostoutlier"]])){
        #cat("j=",j, "\n")
        dg[["almostoutlier"]][[j]] <- strsplit(dg[["almostoutlier"]][[j]],"\\s+")
        dg[["almostoutlier"]][[j]][[1]]<-as.numeric(dg[["almostoutlier"]][[j]][[1]])
        names(dg[["almostoutlier"]][[j]][[1]])<- c("tval_AO","tval_LS","tval_TC")	
      }}
  }
  
  
  if(slidingspans){
    
    ### S 0.
    if(dg[[which(names(dg) %in% "ss.options")]][[1]]!="-"){
      ss.options <- as.list(as.numeric(unlist(strsplit(dg[[which(names(dg) %in% "ss.options")]][[1]],"\\s+"))))
      names(ss.options)<-c("nSpans","lSpans","period1span1","year1span1")
      dg[[which(names(dg) %in% "ss.options")]] <- ss.options
    }
    
    
    
#ss.out <- c("ss.stabseas","ss.movseas","ss.m7","ss.idseas")
#suppressWarnings(for(i in which(names(dg) %in% ss.out)){
#			 	out.vec <- unlist(strsplit(dg[[i]][[1]],"\\s+"))
#			 	if(length(which(is.na(as.numeric(out.vec))))==0){
#						 		dg[[i]][[1]] <- as.numeric(out.vec)	
#						 	}else{	
#						 	dg[[i]][[1]] <- out.vec
#						 	}
#			 })
    
    ss.out <- c("ss.stabseas","ss.movseas","ss.m7","ss.idseas")
    ss.S0 <-dg[which(names(dg) %in% ss.out)]
    rn.S0<-names(dg)[which(names(dg) %in% ss.out)]
    ss.S0 <- strsplit(unlist(ss.S0),"\\s+")			
    dims <- lapply(ss.S0,length)
    if(any(dims>1)){
      rn.S0 <- rn.S0[which(dims>1)]
      ss.S0 <- as.data.frame(t(data.frame(do.call(rbind,ss.S0[which(dims>1)]))))
      ss.S0[,-which(rn.S0=="ss.idseas")]<-apply(ss.S0[,-which(rn.S0=="ss.idseas")],2,as.numeric)	
      row.names(ss.S0)<- paste("span",1:(dim(ss.S0)[1]),sep="")
      colnames(ss.S0)<-gsub("ss.","",rn.S0)
      
      
      names(dg)[[which(names(dg)=="ss.stabseas")]]<-"ss.seasTests"
      dg <- dg[-which(names(dg) %in% c("ss.movseas","ss.m7","ss.idseas"))]
      dg[[which(names(dg)=="ss.seasTests")]]  <- ss.S0	
    }
    
    
    
    
    ### S 1.
    
    if(any(dg[[which(names(dg)=="ss.S1")]]!="-")){
      ss.S1split <- dg[[which(names(dg)=="ss.S1")]]
      ss.S1split <- strsplit(unlist(ss.S1split),"\\s+")
      dims <- lapply(ss.S1split,length)
      
      ss.S1 <- data.frame(do.call(rbind,ss.S1split[which(dims==unique(dims)[[1]])]))
      ss.S1[,1]<-as.character(ss.S1[,1])
      ss.S1[,-1]<-apply(ss.S1[,-1],2,as.numeric)
      row.names(ss.S1)<-NULL
      colnames(ss.S1)<-c("period", paste("span",1:(unique(dims)[[1]]-3),sep=""), "maxPercDiff", "allSpans")
      
      ss.summaryS1 <- data.frame(do.call(rbind,ss.S1split[which(dims==unique(dims)[[2]])]),stringsAsFactors=FALSE)
      ss.summaryS1<-apply(ss.summaryS1,2,as.numeric)
      row.names(ss.summaryS1) <- c(paste("span",1:(dim(ss.summaryS1)[1]-1),sep=""),"allSpans")
      colnames(ss.summaryS1)<-c("min","max","range")
      
      dg[[which(names(dg)=="ss.S1")]]  <- list(ss.S1=ss.S1,ss.summaryS1=ss.summaryS1)
    }
    
    ### S 2.
    if(any(dg[[which(names(dg)=="ss.S2")]]!="-")){
      ss.S2 <- dg[[which(names(dg)=="ss.S2")]]
      ss.S2split <- dg[[which(names(dg)=="ss.S2")]]
      ss.S2split <- strsplit(unlist(ss.S2split),"\\s+")
      dims <- lapply(ss.S2split,length)
      if(any(dims>1)){
        ss.S2 <- data.frame(apply(do.call(rbind,ss.S2split[which(dims>1)]),2,as.numeric))
        rn.S2 <- which(c("s2.a.per", "s2.b.per", "s2.c.per", "s2.d.per", "s2.e.per")%in%names(ss.S2split))
        row.names(ss.S2) <- c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")[rn.S2]
        colnames(ss.S2)<-c("nUnstable","nPeriods","percUnstable")
        dg[[which(names(dg)=="ss.S2")]]  <- ss.S2
      }
    }
    
    ### S 3.
#remove histograms
#remove hinges (no standard format)
    if(any(dg[[which(names(dg)=="ss.S3")]]!="-")){
      ss.S3 <- dg[[which(names(dg)=="ss.S3")]]
      ss.S3 <- ss.S3[-grep("thist",names(ss.S3))]
      ss.S3 <- ss.S3[-grep("hinge",names(ss.S3))]
      
      ss.S3split <- strsplit(unlist(ss.S3),"\\s+")
#dims <- lapply(ss.S3split,length)
      ss.S3 <- data.frame(do.call(rbind,ss.S3split),stringsAsFactors=FALSE)
      ss.S3[,-1] <- apply(ss.S3[,-1],2,as.numeric)
      colnames(ss.S3) <- c("periodYear","nBreakdowns","ampd")
      
      ss.S3.new <- list()
      ss.S3.names <- vector()
      for( i in 1:5){
        s3.names <- c("s3.a.brk", "s3.b.brk", "s3.c.brk", "s3.d.brk", "s3.e.brk")
        ss.S3.new.obj <-ss.S3[grep(s3.names[i],rownames(ss.S3)),]
        if(dim(ss.S3.new.obj)[1]!=0){
          rownames(ss.S3.new.obj)<-NULL
          ss.S3.new[[length(ss.S3.new)+1]] <-ss.S3.new.obj
          ss.S3.names[length(ss.S3.names)+1] <- c("a.seasFac","b.td","c.SA","d.period-period","e.year-year")[i]
        }
      }	
      ss.S3 <- ss.S3.new
      names(ss.S3) <- ss.S3.names
      dg[[which(names(dg)=="ss.S3")]]  <- ss.S3
    }
    
  }	
  
  if(history){
    ### R 1. - R 6.
    
    h.tablenames <- c("h.R1","h.R2","h.R4","h.R5","h.R6")
    for(i in 1:length(h.tablenames)){
      if(any(unlist(dg[h.tablenames[i]])!="-")){	
        lags <- unique(substr(names(dg[h.tablenames[i]][[1]]),1,9))	
        if(any(grepl("aarmo",lags)))
          lags <- lags[-grep("aarmo",lags)]
        colnames.lags <- gsub(paste("r0",substr(h.tablenames[i],4,4),".lag0",sep=""),"lag",lags)
        colnames.lags <- gsub(paste("r0",substr(h.tablenames[i],4,4),".lag",sep=""),"lag",colnames.lags)
        colnames.lags <- gsub(paste("r0",substr(h.tablenames[i],4,4),".proj.",sep=""),"proj",colnames.lags)
        
        dates <- date.list <- total.list <- hingeValues.list <- hinges <- list()
        for(j in 1:length(lags)){
          table <- strsplit(unlist(dg[h.tablenames[i]][[1]]),"\\s+")
          table <- table[grep(lags[j],names(table))]
          
          date <- table[-c(grep(".all",names(table)),grep(".hinge",names(table)))]
          date <- do.call(rbind,date)
          rownames(date)<-NULL
          date.list[[length(date.list)+1]]<-date[,2]
          dates[[length(dates)+1]] <- date[,1]
          
          total <- unlist(table[grep(".all",names(table))])
          rownames(total)<-NULL
          total.list[[length(total.list)+1]]<-total
          
          hingeValues <- table[grep(".hinge",names(table))]
          hingeValues <- do.call(rbind,hingeValues)
          hinges[[length(hinges)+1]] <- do.call(rbind,strsplit(rownames(hingeValues),"hinge."))[,2]
#rownames(hingeValues)<-NULL
          hingeValues.list[[length(hingeValues.list)+1]]<-hingeValues
        }
        
        if(length(lags)>1){
#date
          if(all(lapply(dates,length)==length(unique(unlist(dates))))){								
            date <- as.data.frame(cbind(dates[[1]],do.call(cbind,date.list)),stringsAsFactors=FALSE)
            colnames(date)<-c("date",colnames.lags)
            date[,-1]<-apply(date[,-1],2,as.numeric)
            rownames(date)<-NULL
          }else{
            all.dates.sorted <- all.dates <- unique(unlist(dates))
            months <-c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")	
            if(any(all.dates%in%months)){
              all.dates.sorted[1:length(suppressWarnings(which(is.na(as.numeric(all.dates)))))] <- all.dates[suppressWarnings(which(is.na(as.numeric(all.dates))))][match(months,all.dates[suppressWarnings(which(is.na(as.numeric(all.dates))))])]
            }else{
              quarters <- c("1st",  "2nd",  "3rd",  "4th")	
              all.dates.sorted[1:length(suppressWarnings(which(is.na(as.numeric(all.dates)))))] <- all.dates[suppressWarnings(which(is.na(as.numeric(all.dates))))][match(quarters,all.dates[suppressWarnings(which(is.na(as.numeric(all.dates))))])]		
            }	
            all.dates.sorted[(length(suppressWarnings(which(is.na(as.numeric(all.dates)))))+1):length(all.dates.sorted)] <- as.character(sort(suppressWarnings(as.numeric(all.dates[which(!is.na(as.numeric(all.dates)))]))))
            all.dates <- all.dates.sorted
            for(k in 1:length(dates)){
              date.list[[k]] <- c(date.list[[k]],rep(NA,length(all.dates[which(!all.dates%in%dates[[k]])])))
              dates[[k]] <- c(dates[[k]],all.dates[which(!all.dates%in%dates[[k]])])
              date.list[[k]] <- date.list[[k]][ match(all.dates,dates[[k]])]
            }	
            date <- as.data.frame(cbind(all.dates,do.call(cbind,date.list)),stringsAsFactors=FALSE)
            colnames(date)<-c("date",colnames.lags)
            date[,-1]<-apply(date[,-1],2,as.numeric)
            rownames(date)<-NULL
          }
#total
          total <- as.data.frame(do.call(cbind,total.list),stringsAsFactors=FALSE)	
          total[,1:(dim(total)[2])] <- as.numeric(total)
          rownames(total)<-NULL
          colnames(total)<-colnames.lags
#hingeValues			
          if(length(unique(unlist(hinges)))==length(unlist(hinges))/length(lags)){
            hingeValues <- as.data.frame(cbind(unlist(hinges[[1]]),do.call(cbind,hingeValues.list)),stringsAsFactors=FALSE)	
            rownames(hingeValues)<-NULL
            colnames(hingeValues)<-c("hinges",colnames.lags)
            hingeValues[,-1]<-apply(hingeValues[,-1],2,as.numeric)
          }else{
            hingeValues <- as.data.frame(do.call(cbind,lapply(1:length(lags),function(x)cbind(hinges[[x]],hingeValues.list[[x]]))),stringsAsFactors=FALSE)
            rownames(hingeValues)<-NULL
            colnames(hingeValues)[seq(from=1,by=2,length.out=length(lags))]<-rep("hinges",length(lags))
            colnames(hingeValues)[seq(from=2,by=2,length.out=length(lags))]<-colnames.lags	
          }
        }else{ 	
#date
          date <- as.data.frame(cbind(unlist(dates),unlist(date.list)),stringsAsFactors=FALSE)	
          colnames(date)<-c("date",colnames.lags)
          date[,1]<-as.character(date[,1])
          date[,2]<-as.numeric(date[,2])
          rownames(date)<-NULL
#total
          total <- as.data.frame(do.call(cbind,total.list),stringsAsFactors=FALSE)	
          total[1]<-as.numeric(total)
          colnames(total)<-colnames.lags
          rownames(total)<-NULL
#hingeValues
          hingeValues <- as.data.frame(cbind(unlist(hinges[[1]]),do.call(cbind,hingeValues.list)),stringsAsFactors=FALSE)	
          rownames(hingeValues)<-NULL
          colnames(hingeValues)<-c("hinges",colnames.lags)
          hingeValues[,-1]<-as.numeric(hingeValues[,-1])
        }
        dg[[which(names(dg) %in% h.tablenames[i])]] <- list(date=date,total=total,hingeValues=hingeValues)
      }}
    
    
    ### R 7.
    if(file.exists(paste(filename,".","lkh",sep=""))){
      h.R7 <- readLines(con=paste(filename,".","lkh",sep=""),n=-1)
      names.h.R7 <- c("date","loglikelihood","aicc")
      h.R7 <- data.frame(do.call(rbind,lapply(strsplit(h.R7[-(1:2)],"\t"),as.numeric)))
      colnames(h.R7)<-names.h.R7
      dg[[length(dg)+1]]<-h.R7
      names(dg)[[length(dg)]]<-"h.R7"
    }
    
    ### R 8.
    if(file.exists(paste(filename,".","fce",sep=""))){
      h.R8 <- readLines(con=paste(filename,".","fce",sep=""),n=-1)
      names.h.R8 <- gsub("Sum","sum",gsub("\\(|\\)","",unlist(strsplit(h.R8[1],"\t"))))	
      h.R8 <- data.frame(do.call(rbind,lapply(strsplit(h.R8[-(1:2)],"\t"),as.numeric)))
      colnames(h.R8)<-names.h.R8
      rownames(h.R8)<-NULL
      h.meanssfe <- as.data.frame(t(as.numeric(unlist(strsplit(unlist(dg[["h.meanssfe"]]),"\\s+")))))
      colnames(h.meanssfe)<-gsub("sumSqFcstError","lead",names.h.R8[-1])
      rownames(h.meanssfe)<-NULL	
      dg[[length(dg)+1]]<-list(h.R8=h.R8,meanSumSqFcstError=h.meanssfe)
      names(dg)[[length(dg)]]<-"h.R8"
      dg <- dg[-which(names(dg) %in% "h.meanssfe")]
    }
    
  }#end history
  
  ### identify: autocorrealtions of the residuals
  if(file.exists(paste(filename,".","iac",sep=""))){
    iac <- readLines(con=paste(filename,".","iac",sep=""),n=-1)
    i.tables <- grep("\\$diff",iac)
    #grep("\\$sdiff",iac,value=TRUE)
    name.table <- iac.list <- list()
    for(i in 1:length(i.tables)){
      name.table[[length(name.table)+1]] <- gsub("= ","",gsub("\\$","",paste(iac[i.tables[i]],iac[i.tables[i]+1],sep="_")))
      #names.iac <- gsub("\\(|\\)","",unlist(strsplit(iac[i.tables[i]+2],"\t")))	
      names.iac <- c("lag","sample.acf","stderr.acf","Ljung-Box.q","df.q","pval" )
      if(i.tables[i]!=i.tables[length(i.tables)])
        iac.table <- data.frame(do.call(rbind,lapply(strsplit(iac[(i.tables[i]+4):(i.tables[i+1]-1)],"\t"),as.numeric)))
      else
        iac.table <- data.frame(do.call(rbind,lapply(strsplit(iac[(i.tables[i]+4):length(iac)],"\t"),as.numeric)))
      colnames(iac.table)<-names.iac
      rownames(iac.table)<-NULL
      iac.list[[length(iac.list)+1]] <- iac.table
    }
    names(iac.list)<-name.table
    dg[[length(dg)+1]]<-iac.list
    names(dg)[[length(dg)]]<-"rsd.iac"	
  }
  
  ### identify: partial autocorrealtions of the residuals
  if(file.exists(paste(filename,".","ipc",sep=""))){
    ipc <- readLines(con=paste(filename,".","ipc",sep=""),n=-1)
    i.tables <- grep("\\$diff",ipc)
    #grep("\\$sdiff",ipc,value=TRUE)
    name.table <- ipc.list <- list()
    for(i in 1:length(i.tables)){
      name.table[[length(name.table)+1]] <- gsub("= ","",gsub("\\$","",paste(ipc[i.tables[i]],ipc[i.tables[i]+1],sep="_")))
      #names.ipc <- gsub("\\(|\\)","",unlist(strsplit(ipc[i.tables[i]+2],"\t")))	
      names.ipc <- c("lag","sample.pacf","stderr.pacf")
      if(i.tables[i]!=i.tables[length(i.tables)])
        ipc.table <- data.frame(do.call(rbind,lapply(strsplit(ipc[(i.tables[i]+4):(i.tables[i+1]-1)],"\t"),as.numeric)))
      else
        ipc.table <- data.frame(do.call(rbind,lapply(strsplit(ipc[(i.tables[i]+4):length(ipc)],"\t"),as.numeric)))
      colnames(ipc.table)<-names.ipc
      rownames(ipc.table)<-NULL
      ipc.list[[length(ipc.list)+1]] <- ipc.table
    }
    names(ipc.list)<-name.table
    dg[[length(dg)+1]]<-ipc.list
    names(dg)[[length(dg)]]<-"rsd.ipc"	
  }
  
  
  if(transform=="auto" && dg[["transform"]]!="Automatic selection" &&!x11regress){
    dg[["autotransform"]] <- dg[["transform"]] 
    names(dg[["autotransform"]]) <-"autotransform"
    dg[["transform"]][[1]]<-"Automatic selection"	
    names(dg[["transform"]])<-"transform"
  }
  
  if(any(grepl("derived.coef",ls()))){
    dg[[length(dg)+1]]<-derived.coef
    names(dg)[[length(dg)]]<-"derived.coef"
  }
  
  if(dg[["finalseasonalma"]]!="-"){
    dg[["seasonalma"]] <- c(dg[["seasonalma"]],dg[["finalseasonalma"]])
  }
  dg <- dg[-which(names(dg)=="finalseasonalma")]
  
  if(dg[["finalsamode"]]!="-"){
    dg[["samode"]] <- c(dg[["samode"]],dg[["finalsamode"]])
  }
  dg <- dg[-which(names(dg)=="finalsamode")]
  
  if(dg[["id.rsdseas"]]=="-")
    dg[["id.rsdseas"]]<-"Residual seasonality present"	
  
  if(is.null(outlier) &&!x11regress){
    dg[[length(dg)+1]]<-"No outlier detection performed"
    names(dg)[[length(dg)]]<-"ifout"
  }else if(!is.null(outlier) &&!x11regress){
    dg[[length(dg)+1]]<-"Outlier detection performed"
    names(dg)[[length(dg)]]<-"ifout"	
  }
  
#filename <- "M:/Meraner/Workspace/Saisonbereinigung_Test/gra/Rout"
#file<-"M:/Meraner/Workspace/Saisonbereinigung_Test/Rout"
#!if exists notwendig
  if(file.exists(paste(filename,".","acf",sep=""))){
#sample autocorrelations of residuals
    acf <- readLines(con=paste(filename,".","acf",sep=""),n=-1)
#names.acf <- unlist(strsplit(acf[1],"\t"))
    names.acf <- c("lag","sample.acf","stderr.acf","Ljung-Box.q","df.q","pval")
    acf <- data.frame(do.call(rbind,lapply(strsplit(acf[-(1:2)],"\t"),as.numeric)))
    colnames(acf)<-names.acf
    dg[[length(dg)+1]]<-acf
    names(dg)[[length(dg)]]<-"rsd.acf"
  }
  
  if(file.exists(paste(filename,".","pcf",sep=""))){
#sample partial autocorrelations of residuals
    pacf <- readLines(con=paste(filename,".","pcf",sep=""),n=-1)
#names.pacf <- unlist(strsplit(pacf[1],"\t"))
    names.pacf <- c("lag","sample.pacf","stderr.pacf")
    pacf <- data.frame(do.call(rbind,lapply(strsplit(pacf[-(1:2)],"\t"),as.numeric)))
    colnames(pacf)<-names.pacf
    dg[[length(dg)+1]]<-pacf
    names(dg)[[length(dg)]]<-"rsd.pacf"
  }
  
  if(file.exists(paste(filename,".","ac2",sep=""))){
#sample autocorrelations of squared residuals
    acf2 <-readLines(con=paste(filename,".","ac2",sep=""),n=-1)
    names.acf2 <- unlist(strsplit(acf2[1],"\t"))
    names.acf2 <- c("lag","sample.acf2","stderr.acf2","Ljung-Box.q","df.q","pval")
    acf2 <- data.frame(do.call(rbind,lapply(strsplit(acf2[-(1:2)],"\t"),as.numeric)))
    colnames(acf2)<-names.acf2[1:dim(acf2)[2]]#Box Ljung fehlt fuer acf2
#Box.test(,lag=1,type="Ljung-Box",)
    dg[[length(dg)+1]]<-acf2
    names(dg)[[length(dg)]]<-"rsd.acf2"
  }
  
  ##Forecasts als data.frame
#if(file.exists(paste(filename,".","fct",sep=""))){
  ##sample autocorrelations of residuals
#	fct <- readLines(con=paste(filename,".","fct",sep=""),n=-1)
#	names.fct <- unlist(strsplit(fct[1],"\t"))
  ##	names.fct <- c("..")
#	fct <- data.frame(do.call(rbind,lapply(strsplit(fct[-(1:2)],"\t"),as.numeric)))
#	colnames(fct)<-names.fct
  ##	dg[[length(dg)+1]]<-fct
  ##	names(dg)[[length(dg)]]<-"fct"
#out$forecast <- fct
#}
#
  ##Backcasts als data.frame
#if(file.exists(paste(filename,".","bct",sep=""))){
  ##sample autocorrelations of residuals
#	bct <- readLines(con=paste(filename,".","bct",sep=""),n=-1)
#	names.bct <- unlist(strsplit(bct[1],"\t"))
  ##	names.bct <- c("..")
#	bct <- data.frame(do.call(rbind,lapply(strsplit(bct[-(1:2)],"\t"),as.numeric)))
#	colnames(bct)<-names.bct
  ##	dg[[length(dg)+1]]<-bct
  ##	names(dg)[[length(dg)]]<-"bct"
#out$backcast <- bct	
#}
  dg[[length(dg)+1]]<- file
  names(dg)[[length(dg)]]<-"tsName"	
#out <- list()
  dg[[length(dg)+1]]<- freq_series
  names(dg)[[length(dg)]]<-"frequency"	
  
  if(any(grepl("span:",udg,fixed=TRUE))){
    span <- grep("span:",udg,fixed=TRUE,value=TRUE)
    if(length(span)>1){
      dg[[length(dg)+1]]<- span
      names(dg)[[length(dg)]]<-"span"
    }else{
      span <- str_trim(unlist(strsplit(span,":")))		
      dg[[length(dg)+1]]<- span[2]
      names(dg)[[length(dg)]]<-"span"
    }
  }	
  
  out[["dg"]] <- list()
  out[["dg"]] <- dg
  
#  out$seats <- seats
  out$file <- file
  out$tblnames <- tblnames
  out$Rtblnames <- Rtblnames
  class(out) <- "x12work"
  out  
  
  
}