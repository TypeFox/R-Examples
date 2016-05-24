ms_sev <- function(data, type="global_armss", range=2, matrix=F, omsss=FALSE ){
  
  
  names <- colnames(data)
  
  if(type=="local_armss"){
    
    if(!is.na(match("edss",names)) && !is.na(match("ageatedss", names)))
      results <- local_armss(data, range)
    else stop("Variable name edss or ageatedss missing in data")
  }else if(type=="global_armss"){
    if(!is.na(match("edss",names)) && !is.na(match("ageatedss", names)))
      results <- global_armss(data, matrix) 
    else stop("Variable name edss or ageatedss missing in data")
  }else if(type=="local_msss"){
    if(!is.na(match("edss",names)) && !is.na(match("dd", names)))
      results <- local_msss(data, range)
    else stop("Variable name edss or dd missing in data")
    
  }else if(type=="global_msss"){
    if(!is.na(match("edss",names)) && !is.na(match("dd", names)))
      results <- global_msss(data, matrix, omsss)
    else stop("Variable name edss or dd missing in data")
    
    
  }
  
  
  return(results)
}


global_msss <- function(data, matrix=FALSE, omsss=FALSE){
  
  oldMsss <- NULL
  globalMsss <- NULL
  
  if(omsss==TRUE){
    data(oldMsss, package="ms.sev", envir=environment())

  }
  
  data(globalMsss, package="ms.sev", envir=environment())

  
  data$uGMSSS <- NA
  if(omsss) data$oGMSSS <- NA
  ## truncate the disease duration to 30 years
  
  data$ddT <- data$dd
  data$ddT[data$ddT>29]<-30
  
  data$tmp <- data$edss/0.5+1
  
  data$tmp[data$tmp==1] <- 2
  
  for(i in 1:nrow(data)){
    
    if(floor(data$ddT[i])!=0 && data$edss[i]<10 & !is.na(data$ddT[i]) & !is.na(data$edss[i])){
      data$uGMSSS[i] <- globalMsss[floor(data$ddT[i]), data$tmp[i]]
      if(omsss) {
        data$oGMSSS[i] <- oldMsss[floor(data$ddT[i]), data$tmp[i]]
      }
    }
    
  }
  data$ddT <- NULL
  data$tmp <- NULL
  if(matrix==TRUE){
    if(!omsss){
      out = list("data"=data, "globalMsss"=globalMsss)
    }else{
      out = list("data"=data, "globalMsss"=globalMsss, "oldGlobalMsss"=oldMsss)
    }
  }else{
    out = list("data"=data)
  }
  
  return(out)
}

## returns a matrix for the local msss
## keep id
## assumes variables are names dd and edss
local_msss <- function(data, range=2){
  
  data$lid <- rep(1:nrow(data))
  data$lMSSS <- NA
  ## truncate the data to age 30 as max
  data$dd[data$dd>30] <- 30
  ## sort the data based on diseas duration and edss in a decreasing order
  ## initial check of the data should make sure the parameter name ddl and edss is always correct
  data.sorted <- data[order(data$dd, data$edss, decreasing=FALSE, na.last=NA),]
  ## we are only interested in the first 30 years. 
  for(i in 1:30){
    
    ## noone is interesting before they've had the disease for at least one year
    if(i-range <= 0) lower <- 1
    else lower <- i-range
    
    ## and 29 is the maximum disease duration
    if(i+range > 30 ) upper <- 30
    else upper <- i+range
    
    ## create a subset of all that have had the disease for i +/- range years.
    # include the lower limit as well
    
    data.range <- subset(data.sorted, (data.sorted$dd >= lower) & (data.sorted$dd <= upper))
    ## assign an averaged rank based on the edss score. 
    data.range$rank <- rank(data.range$edss, ties.method="average")
    ## keep only individuals during the i:th year of disease
    year <- subset(data.range, floor(data.range$dd)==i)
    
    ## do only this for the years where we have any patients
    if(dim(year)[1]>0){
      # get number of patients the rank is based on
      n <- dim(data.range)[1] #*dim(data.range)[2]
      
      year$msss <- (year$rank/(1+n))*10
      
      ## paste back the msss some how... 
      ## not pretty, but will work for now
      for(i in 1:nrow(year)){
        #if(!is.na(data$dd[i])){ ## need to add edss??
          data[data$lid[year$lid[i]],]$lMSSS <- round(year$msss[i],2)
       # }
        
        
      }
      
    }
    
  }
  data$lid <- NULL
  return(data)
}
local_armss <- function(data, range=2){
  
  data$lid <- rep(1:nrow(data))
  data$lARMSS <- rep(NA, nrow(data))
  
  
  ## sort the data based on diseas duration and edss in a decreasing order
  ## initial check of the data should make sure the parameter name ddl and edss is always correct
  data.sorted <- data[order(data$ageatedss, data$edss, decreasing=FALSE, na.last=NA),]
  
  ## getting the max age 
  top <- ceiling(max(data$ageatedss, na.rm=T))
  for(i in 1:top){
    
    ## noone is interesting before they've had the disease for at least one year
    if(i-range <= 0) lower <- 1
    else lower <- i-range
    
    ## truncating to the max age
    if(i+range > top ) upper <- top
    else upper <- i+range
    
    ## create a subset of all that are within the age-range
    data.range <- subset(data.sorted, (data.sorted$ageatedss >= lower) & (data.sorted$ageatedss <= upper))
    
    ## assign an averaged rank based on the edss score. 
    data.range$rank <- rank(data.range$edss, ties.method="average")
    
    ## keep only individuals during the i:th year of disease
    year <- subset(data.range, floor(data.range$ageatedss)==i)
    
    ## do only this for the years where we have any patients
    if(dim(year)[1]>0){
      # get number of patients the rank is based on
      n <- dim(data.range)[1] #*dim(data.range)[2]
      
      year$armss <- (year$rank/(1+n))*10
      
      ## paste back the msss some how... 
      ## not pretty, but will work for now
      for(i in 1:nrow(year)){
        
        if(!is.na(data$edss[i]) & !is.na(data$ageatedss[i])){
          data[data$lid[year$lid[i]],]$lARMSS <- round(year$armss[i],2)
        }
        
        
      }
      
    }
    
  }
  data$lid <- NULL
  
  return(data)
  
}
global_armss <- function(data, matrix=FALSE){
  globalArmss <- NULL

  data(globalArmss, package="ms.sev", envir=environment())
  #colnames(globalArmss) <- c(0,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5)
  
  data$gARMSS <- NA
 
  data$tmp <- data$edss/0.5+1
  data$tmp[data$tmp==1] <- 2
  data$ageT <- data$ageatedss-17
  for(i in 1:nrow(data)){
    ## check first that the age is in the range
    
    if(data$ageT[i] <= 58 & data$ageT[i] >= 1 & data$edss[i]<10 & !is.na(data$ageT[i]) & !is.na(data$edss[i])) {
      
      data$gARMSS[i] <- globalArmss[floor(data$ageT[i]), data$tmp[i]]
    }
  }
  data$tmp <- NULL
  data$ageT <- NULL
  
  if(matrix==TRUE){
    out = list("data"=data, "globalArmss"=globalArmss)
  }else{
    out = list("data"=data)
  }
  
  return(out)
}
