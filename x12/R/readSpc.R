readSpc <- function(file,filename=TRUE){
  if(length(file)==1)
    return(readSpcS(file))
  else{
    tsl <- list()
    for(i in 1:length(file)){
      tsl[[i]] <- readSpcS(file[i])
    }
   return(new("x12Batch",tsList=tsl))
  }
  
}
readSpcS <- function(file,filename=TRUE){
  whichgrep <- function(x,txt){
    which(unlist(lapply(x,function(y)length(grep(txt,y,ignore.case=TRUE))>0)))
  }
  gsub1 <- function(x,word){
    x <- gsub(word,"",x,ignore.case=TRUE)
    x <- gsub("=","",x)
    x <- gsub("\\(","",x)
    x <- gsub("\\)"," ",x)
    x <- gsub("\\'","",x)
    x <- gsub("\\\"","",x)
    x <- gsub("\\{","",x)
    x <- gsub("\\}","",x)
    x <- unlist(strsplit(x," "))
    x[x!=""]
  }
  gsub2 <- function(x,word){
    x <- gsub(word,"",x,ignore.case=TRUE)
    x <- gsub("=","",x)
    x <- gsub("\\(","",x)
    x <- gsub("\\)","",x)
    x <- gsub("\\{","",x)
    x <- gsub("\\}","",x)
    x <- gsub("\\'","",x)
    x <- gsub("\\\"","",x)
    gsub(" ","",x)
  }
  yes <- function(x)
    length(gregexpr("yes",x,ignore.case=TRUE))>0
  getPart <- function(x,word){
    startreg <- whichgrep(x,word) 
    if(length(startreg>0)){
      endreg <- whichgrep(x,"}")
      endreg <- endreg[endreg>=startreg][1]
      return(c(startreg,endreg))
    }else
      return(vector())
  }
  para <- new("x12Parameter")
  Lines <- readLines(file)
###CLEANING COMMENTS
#comment line
  todel <- NULL
  for(i in 1:length(Lines)){
    Lines[i] <- str_trim(Lines[i])
    if(substr(Lines[i],0,1)=="#")
      todel <- c(todel,i)
  }
  if(!is.null(todel))
    Lines <- Lines[-todel]
#end of line comments
  commentX <- whichgrep(Lines,"#")
  for(l in commentX){
    Lines[l] <- substr(Lines[l],0,which(unlist(strsplit(Lines[l],""))=="#")-1)
  }
###SERIES
  ind <- getPart(Lines,"series")
  series <- Lines[ind[1]:ind[2]]
  Lines <- Lines[-c(ind[1]:ind[2])]  
#Series_name  
  if(!filename){
    namX <- whichgrep(series,"name")
    namX2 <- whichgrep(series,"title")
    if(length(namX)>0){
      name <- gsub2(series[namX],"name")
    }else if(length(namX2)>0){
      name <- gsub2(series[namX2],"title")
    }else{
      name <- "Series_1"
    }
  }else{
    name <- substr(file,0,nchar(file)-4)
  }
  perX <- whichgrep(series,"period")
  if(length(perX)>0){
    period <- as.numeric(gsub1(series[perX],"period"))
    if(!period%in%c(4,12))
      stop("Period argument wring?!")
  }else{
    period <- 12
  }
  startX <- whichgrep(series,"start")
  if(length(startX)>0){
    start <- gsub1(series[startX],"start")
    start <- unlist(strsplit(start,"\\."))
    if(!period%in%c(4,12))
      stop("Period argument wring?!")
  }else{
    start <- 1
  }
  start <- as.numeric(start)
  dataX <-  whichgrep(series,"data")
  fileX <-  whichgrep(series,"file")
  if(length(dataX)>0){
    dataXEnd <-  whichgrep(series,")")[1]
    dataXEnd <- dataXEnd[dataXEnd>=dataX][1]
    data <- c()
    for(i in dataX:dataXEnd){
      line <- gsub1(series[i],"data")
      if(length(line)>0){
        data <- c(data,unlist(strsplit(line," ")))
      }
    }
    data <- as.numeric(data[data!=""])
  }else if(length(fileX)>0){
    dataFile <- gsub2(series[fileX] ,"file")
    dataFile <- gsub("\\\\","/",dataFile)
    dataFile <- gsub("//","/",dataFile)
    data <- scan(dataFile)
  }else{
    warning("No data found, only the x12Parameter object will be created")
  }
  data <- ts(data,start=start,frequency=period)
  mspanX <-  whichgrep(series,"modelspan")
  if(length(mspanX)>0){
    modelspan <- gsub2(series[mspanX],"modelspan")
    modelspan <- unlist(lapply(unlist(strsplit(modelspan,",")),function(x)strsplit(x,"\\.")))
    para <- setP(para,list(series.modelspan=as.numeric(modelspan)))
  }
  spanX <-  whichgrep(series,"span")
  spanX <- spanX[!spanX%in%mspanX]
  if(length(spanX)>0){
    span <- gsub2(series[spanX],"span")
    span <- unlist(lapply(unlist(strsplit(span,",")),function(x)strsplit(x,"\\.")))
    para <- setP(para,list(series.span=as.numeric(span)))
  }
###TRANSFORM
  ind <- getPart(Lines,"transform")
  if(length(ind)>0){
    trans <- Lines[ind[1]:ind[2]]
    Lines <- Lines[-c(ind[1]:ind[2])]   
    func <- whichgrep(trans,"function")
    transF <- gsub2(trans[func],"function")
    para <- setP(para,list(transform.function=transF))
    adj <- whichgrep(trans,"adjust")
    trans <- gsub2(trans[adj],"adjust")
    if(length(trans>0))
      para <- setP(para,list(transform.adjust=trans))
  }
###OUTLIER
  ind <- getPart(Lines,"outlier")
  if(length(ind)>0){
    outlier <- Lines[ind[1]:ind[2]]
    Lines <- Lines[-c(ind[1]:ind[2])]  
    typX <- whichgrep(outlier,"types")
    if(length(typX)>0){
      types <- outlier[typX]
      types <- gsub("types","",types,ignore.case=TRUE)
      types <- gsub("=","",types)
      types <- gsub("\\(","",types)
      types <- gsub("\\)","",types)
      types <- unlist(strsplit(types," "))
      types <- types[types!=""]
      para <- setP(para,list(outlier.types=types))
    }
    spanX <-  whichgrep(outlier,"span")
    if(length(spanX)>0){
      out_span <- gsub2(outlier[spanX],"span")
      out_span <- unlist(lapply(unlist(strsplit(out_span,",")),function(x)strsplit(x,"\\.")))
      para <- setP(para,list(outlier.span=as.numeric(out_span)))
    }
    methodX <-  whichgrep(outlier,"method")
    if(length(methodX)>0){
      out_meth <- gsub2(outlier[methodX],"method")
      para <- setP(para,list(outlier.method=out_meth))
    }
    critX <-  whichgrep(outlier,"critical")
    if(length(critX)>0){
      crit <- gsub1(outlier[critX],"critical")
      if(length(crit)==1)
        crit <- as.numeric(crit)
      else{
        crit <- as.list(as.numeric(crit))
        names(crit) <- c("AO","LS","TC")[1:length(crit)]
      }
      para <- setP(para,list(outlier.critical=crit))
    }
  } 
#REGRESSION
  ind <- getPart(Lines,"regression")
  if(length(ind)>0){
    regression <- Lines[ind[1]:ind[2]]
    Lines <- Lines[-c(ind[1]:ind[2])]
    varX <-  whichgrep(regression,"variables")
    if(length(varX)>0){
      variables <- gsub1(regression[varX],"variables")
      para <- setP(para,list(regression.variables=variables))
    }
    centeruserX <-  whichgrep(regression,"centeruser")
    if(length(centeruserX)>0){
      centeruser <- gsub1(regression[centeruserX],"centeruser")
      para <- setP(para,list(regression.centeruser=centeruser))
      regression <- regression[-centeruserX]
    }
    usertypeX <-  whichgrep(regression,"usertype")
    if(length(usertypeX)>0){
      usertype <- gsub1(regression[centeruserX],"usertype")
      para <- setP(para,list(regression.usertype=usertype))
      regression <- regression[-usertypeX]
    }
    userX <-  whichgrep(regression,"user")
    if(length(userX)>0){
      reguser <- gsub1(regression[userX],"user")
      para <- setP(para,list(regression.user=reguser))
    }
    fileX <-  whichgrep(regression,"file")
    if(length(fileX)>0){
      regfile <- gsub1(regression[fileX],"file")
      para <- setP(para,list(regression.file=regfile))
    }
  }
#ARIMA
  ind <- getPart(Lines,"arima")
  if(length(ind)>0){
    arima <- Lines[ind[1]:ind[2]]
    Lines <- Lines[-c(ind[1]:ind[2])]
    modelX <-  whichgrep(arima,"model")
    
    if(length(modelX)>0){
      #aa <<- arima[modelX]
      arr <- arima[modelX]
      arr <- gsub(","," ",arr)
      model <- as.numeric(gsub1(arr,"model"))
      if(length(model)==6)
        para <- setP(para,list(arima.model=model[1:3],arima.smodel=model[4:6]))
      else if(length(model)==3)
        para <- setP(para,list(arima.model=model[1:3],arima.smodel=c(0,0,0)))
      else
        stop("Problem in reading ARIMA specification")
    }
  }
#AUTOMDL
  ind <- getPart(Lines,"automdl")
  if(length(ind)>0){
    para <- setP(para,list=(automdl=TRUE))
    automdl <- Lines[ind[1]:ind[2]]
    Lines <- Lines[-c(ind[1]:ind[2])]
    orderX <-  whichgrep(automdl,"automdl.maxorder")
    if(length(orderX)>0){
      maxorder <- as.numeric(unlist(strsplit(gsub1(automdl[orderX],"maxorder"),",")))
      para <- setP(para,list(maxorder=maxorder))
    }
    diffX <-  whichgrep(automdl,"automdl.maxdiff")
    if(length(diffX)>0){
      maxdiff <- as.numeric(unlist(strsplit(gsub1(automdl[diffX],"maxdiff"),",")))
      para <- setP(para,list(maxdiff=maxdiff))
    }
    acceptX <-  whichgrep(automdl,"acceptdefault")
    if(length(acceptX)>0){
      accept <- yes(gsub1(automdl[acceptX],"acceptdefault"))
      para <- setP(para,list(automdl.acceptdefault=accept))
    }
    balancedX <-  whichgrep(automdl,"balanced")
    if(length(balancedX)>0){
      balanced <- yes(gsub1(automdl[balancedX],"balanced"))
      para <- setP(para,list(automdl.balanced=balanced))
    }
  }else{
    para <- setP(para,list(automdl=FALSE))
  }
#ESTIMATE
  ind <- getPart(Lines,"estimate")
  if(length(ind)>0){
    para <- setP(para,list=(estimate=TRUE))
    estimate <- Lines[ind[1]:ind[2]]
    Lines <- Lines[-c(ind[1]:ind[2])]
    oosX <-  whichgrep(estimate,"outofsample")
    if(length(oosX)>0){
      oos <- yes(gsub1(estimate[oosX],"outofsample"))
      if(oos)
        para <- setP(para,list(estimate=TRUE))
      para <- setP(para,list(estimate.outofsample=oos))
    }
  }
#SLIDINGSSPAN
  ind <- getPart(Lines,"slidingspans")
  if(length(ind)>0){
    para <- setP(para,list=(slidingspans=TRUE))
    Lines <- Lines[-c(ind[1]:ind[2])]
  }
#Forecast
  ind <- getPart(Lines,"forecast")
  if(length(ind)>0){
    forecast <- Lines[c(ind[1]:ind[2])]
    Lines <- Lines[-c(ind[1]:ind[2])]
    maxleadX <-  whichgrep(forecast,"maxlead")
    if(length(maxleadX)>0){
      fc_y <- as.numeric(gsub1(forecast[maxleadX],"maxlead"))/period
      para <- setP(para,list(forecast_years=fc_y))
    }
    maxbackX <-  whichgrep(forecast,"maxback")
    if(length(maxbackX)>0){
      bc_y <- as.numeric(gsub1(forecast[maxbackX],"maxback"))/period
      para <- setP(para,list(backcast_years=bc_y))
    }
    probX <-  whichgrep(forecast,"PROBABILITY")
    if(length(probX)>0){
      prob_conf <- as.numeric(gsub1(forecast[probX],"PROBABILITY"))
      para <- setP(para,list(forecast_conf=prob_conf))
    }
    
  }
#X11
  ind <- getPart(Lines,"x11")
  if(length(ind)>0){
    x11 <- Lines[c(ind[1]:ind[2])]
    Lines <- Lines[-c(ind[1]:ind[2])]
    sigmalimX <-  whichgrep(x11,"sigmalim")
    if(length(sigmalimX)>0){
      sigmalim <- as.numeric(unlist(strsplit(gsub2(x11[sigmalimX],"sigmalim"),",")))
      para <- setP(para,list(x11.sigmalim=sigmalim))
    }
    calsigX <-  whichgrep(x11,"calendarsigma")
    if(length(calsigX)>0){
      calendarsigma <- gsub1(x11[calsigX],"calendarsigma")
      para <- setP(para,list(x11.calendarsigma=calendarsigma))
    }
    modeX <-  whichgrep(x11,"mode")
    if(length(modeX)>0){
      samode <- gsub1(x11[modeX],"mode")
      para <- setP(para,list(x11.samode=samode))
    }
    seasonalmaX <-  whichgrep(x11,"seasonalma")
    if(length(seasonalmaX)>0){
      seasonalma <- gsub1(x11[seasonalmaX],"seasonalma")
      para <- setP(para,list(x11.seasonalma=seasonalma))
    }
    trendmaX <-  whichgrep(x11,"trendma")
    if(length(trendmaX)>0){
      trendma <- gsub1(x11[trendmaX],"trendma")
      para <- setP(para,list(x11.trendma=trendma))
    }
    appendfcstX <-  whichgrep(x11,"appendfcst")
    if(length(appendfcstX)>0){
      appendfcst <- yes(gsub1(x11[appendfcstX],"appendfcst"))
      para <- setP(para,list(x11.appendfcst=appendfcst))
    }
    appendbcstX <-  whichgrep(x11,"appendbcst")
    if(length(appendbcstX)>0){
      appendbcst <- yes(gsub1(x11[appendbcstX],"appendbcst"))
      para <- setP(para,list(x11.appendbcst=appendbcst))
    }
    excludefcstX <-  whichgrep(x11,"excludefcst")
    if(length(excludefcstX)>0){
      excludefcst <- yes(gsub1(x11[excludefcstX],"excludefcst"))
      para <- setP(para,list(x11.excludefcst=excludefcst))
    }
    finalX <-  whichgrep(x11,"final")
    if(length(finalX)>0){
      final <- gsub1(x11[finalX],"final")
      para <- setP(para,list(x11.final=final))
    }
  } 
  cat("File:",file,"processed\n")
  return(new("x12Single",ts=data,x12Parameter=para,tsName=name))
}