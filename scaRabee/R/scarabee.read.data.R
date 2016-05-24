
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

scarabee.read.data <- function(files=NULL, method=NULL){
  
  ## Checks inputs
  if (is.null(files)){
    stop('files argument is NULL.')
  }
  
  if (is.null(files$data)){
    stop('files argument does not have any data level or files$data is NULL.')
  }
  
  if (!file.exists(files$data)){
    stop('data file does not exist.')
  }
  
  if (is.null(method)){
    stop('method argument is NULL.')
  }
  
  if (length(method)!=1 | !method%in%c('population','subject')){
    stop('method argument should be \'population\' or \'subject\'.')
  }
  
  ## Read-in data file
  cat('Processing data file:\n')
  data <- read.csv(file=files$data,
                   header=TRUE,
                   as.is=TRUE,
                   na.string=c('.',NA))
  
  # Check data names
  if (paste(names(data),collapse='')!=paste(toupper(names(data)),collapse=''))
    stop('headers in data file should be in upper case.')
  
  names(data) <- toupper(names(data))
  
  ## Check dimension
  if (size(data,1)<1)
    stop('the data file does not contain any data.')
  
  if (length(names(data))!=length(unique(names(data))))
    stop('the data file should not contain duplicate variables names.')
  
  ## Check variable headers
  if (!'OMIT'%in%names(data))
    data <- cbind(data.frame(omit=0),data)
  
  required.vars <- c('OMIT','TRT','ID','TIME','AMT','RATE',
                     'CMT','EVID','DV','DVID','MDV')
  
  if(any(is.na(match(required.vars,names(data))))){
    which.missing <- required.vars[is.na(match(required.vars,names(data)))]
    stop(paste('the dataset does not contain the required variable(s):\n  ',
               paste(which.missing,collapse=' '),sep=''))
  }
  
  ## Check for reserved variable names
  reserved.names <- c('firstTIME','newTIME')
  if (any(names(data)%in%reserved.names)){
    stop(paste('the dataset should not contain variables using',
               'reserved names:\n  ',
               paste(reserved.names,collapse=' '),sep=''))
  }
  
  ## Check numerical variables
  num.vars <- required.vars[-c(1,2)]
  numvar.class <- sapply(num.vars,function(x,...) class(data[,x]),data)
  if (!all(numvar.class=="numeric" | numvar.class=="integer")){
    nonnumvars.index <- which(!(numvar.class=="numeric" |
                                numvar.class=="integer"))
    stop(paste('\  non numerical values are contained in the ',
               'following variable(s):\n  ',
               paste(num.vars[nonnumvars.index],collapse=' '),sep=''))
  }
  
  ## Coerced reserved variables to integers
  data$OMIT <- as.integer(data$OMIT)
  data$ID <- as.integer(data$ID)
  data$CMT <- as.integer(data$CMT)
  data$EVID <- as.integer(data$EVID)
  data$MDV <- as.integer(data$MDV)
  data$DVID <- as.integer(data$DVID)
  
  ## Check for illegal missing
  if (sum(is.na(data$ID))>0)
    stop('ID variable cannot be NA or \'.\'.')
  if (sum(is.na(data$CMT))>0)
    stop('CMT variable cannot be NA or \'.\'.')
  if (sum(is.na(data$EVID))>0)
    stop('EVID variable cannot be NA or \'.\'.')
  if (sum(is.na(data$MDV))>0)
    stop('MDV variable cannot be NA or \'.\'.')
  
  ## Check ID column contain sequence of continuous integers
  uniqueIDs <- unique(data$ID)
  minID <- min(uniqueIDs)
  maxID <- max(uniqueIDs)
  
  if (minID!=1 | length(uniqueIDs)!=maxID)
    stop(paste('the subject IDs do not consist of a ',
               'continuous sequence of integers starting from 1.',sep=''))
  
  ## Check EVID
  uniqueEVIDs <- sort(unique(data$EVID))
  
  if (!all(uniqueEVIDs%in%c(0,1)))
    stop('EVID variable should be set to 0 or 1.')
  
  ## Check DVID content
  if (length(which(is.na(data$DVID)))> 0)
    data$DVID[which(is.na(data$DVID))] <- 0
  dosingDVID <- unique(data$DVID[which(data$EVID==1)])
  obsDVID <- unique(data$DVID[which(data$EVID==0)])
  
  if (length(which(data$EVID==1)) > 0 && 
      (length(dosingDVID)!=1 | dosingDVID[1]!=0))
    stop('DVID variable must be set to 0 for dosing events.')
  
  if (any(obsDVID[1]==0))
    stop('DVID variable must not be set to 0 for observation events.')
  
  ## Check MDV
  uniqueMDVs <- sort(unique(data$MDV))
  
  if (!all(uniqueMDVs%in%c(0,1)))
    stop('MDV variable should be set to 0 or 1.')
  
  dosingMDV <- unique(data$MDV[which(data$EVID==1)])
  
  if (length(which(data$EVID==1)) > 0 && 
      (length(dosingMDV)!=1 | dosingMDV[1]!=1))
    stop('MDV variable must be to 1 for dosing events.')
  
  ## Exclude any data record with OMIT not equal to 0 or NA
  if (length(which(is.na(data$OMIT))) > 0)
    data$OMIT[which(is.na(data$OMIT))] <- 0
  
  if (length(which(data$OMIT==0)) > 0) {
    data <- data[which(data$OMIT==0),]
  } else {
    stop('all records set for omission.')
  }
  data <- data[,-which(names(data)=="OMIT")]
  required.vars <- required.vars[-1]
  
  ## Check that rate is not set to -1 and -2 in the same dataset
  if (all(c(-1,-2)%in%unique(data$RATE))){
    stop('the RATE variable cannot be set to -1 and -2 in the dataset.')
  }
  
  # Sort data by subject, treatment, time, and DVID
  data <- data[order(data$ID,
                     data$TRT,
                     data$TIME,
                     data$DVID),]
  
  row.names(data) <- 1:size(data,1)
  
  # Assess time variable (assume the time of the first record with same ID and
  # TRT should be set to 0)
  uniques <- data[!duplicated(data[,c('ID','TRT')]),]
  names(uniques)[which(names(uniques)=='TIME')] <- 'firstTIME'
  
  data$unique <- paste(data$ID,data$TRT)
  uniques$unique <- paste(uniques$ID,uniques$TRT)
  
  data <- merge(data,uniques[,c('firstTIME','unique')],by='unique')[,-1]
  data$newTIME <- data$TIME - data$firstTIME
  
  # Create a new data file if TIME and newTIME are different
  new.data.file <- FALSE
  
  if (!all(data$TIME==data$newTIME)){
    cat(paste('  Note: The time of the first event in at least one ',
              'treatment period for at\n        least one individual was not 0.',
              ' The modified dataset created for this\n        analysis ',
              'was stored in the working directory:\n        ',files$data,
              '.new\n', sep=''))
    
    data$TIME <- data$newTIME
    data <- data[,-which(names(data)%in%c('firstTIME','newTIME'))]
    
    write.csv(data,
              file=paste(files$data,'.new',sep=''),
              quote=FALSE,
              row.names=FALSE)
    
    new.data.file <- TRUE
    
  }
  if (length(which(names(data)%in%c('firstTIME','newTIME'))>0))
    data <- data[,-which(names(data)%in%c('firstTIME','newTIME'))]
  
  # Find variable indices
  itt.index <- match(c('ID','TRT','TIME'),names(data))
  data.index <- match(c('DVID','DV'),names(data))
  dose.index <- match(c('CMT','AMT','RATE'),names(data))
  cmt.index <- match('CMT',names(data))
  cov.index <- which(!names(data)%in%required.vars)
  
  ## Process data according to method
  ###### Start of subject-based data processing
  if (method=='subject') {
    
    # Sort data by subject, treatment, DVID, and time
    data <- data[order(data$ID,
                       data$TRT,
                       data$DVID,
                       data$TIME),]
    
    ## Start point for real data processing
    uniqueIDs <- unique(data$ID)
    
    # Create backbone of function output list (with one slot per id)
    out.data <- sapply(1:length(uniqueIDs),
                       function(x,tmp) c(tmp,list()),
                       tmp=list())
    names(out.data) <- uniqueIDs
    out.data <- c(out.data, list(ids=uniqueIDs))
    
    for (id in uniqueIDs){
      id.data <- data[which(data$ID==id),]
      idTRTs <- unique(id.data$TRT)
      out.data[[as.character(id)]] <- sapply(1:length(idTRTs),
                                             function(x,tmp) c(tmp,list()),
                                             tmp=list())
      names(out.data[[as.character(id)]]) <- idTRTs
      out.data[[as.character(id)]] <- c(out.data[[as.character(id)]],
                                        list(trts=idTRTs))
      out.data[[as.character(id)]]$id <- id
    }
    
    ## Extract dosing information
    dose.data <- data[which(data$EVID==1),]
    data <- data[which(data$EVID!=1),]
    if (size(data,1)>1) {
      row.names(data) <- 1:size(data,1)
    } else {
      stop(paste('the dataset does not contain any data ',
                 'for the dependent variable(s).',sep=''))
    }
    dose.data <- dose.data[,c(itt.index,dose.index,cov.index)]
    
    # Exclude sample records with MDV not set to 0
    if (length(which(data$MDV==0)) > 0) {
      data <- data[which(data$MDV==0),]
    } else {
      stop('no sample record with MDV=0 was found.')
    }
    
    # Split data between analysis data and covariates data by ID and TRT
    ana.data <- data[,c(itt.index,data.index)]
    cov.data <- data[,c(itt.index,cmt.index,cov.index)]
    
    for (id in uniqueIDs){
      id.ana.data <- ana.data[which(ana.data$ID==id),]
      id.cov.data <- cov.data[which(cov.data$ID==id),]
      uniqueTRTs <- unique(id.ana.data$TRT)
      
      for (treat in uniqueTRTs){
        browser()
        out.data[[as.character(id)]][[as.character(treat)]]$ana <-
          id.ana.data[which(id.ana.data$TRT==treat), names(id.ana.data)[-c(1,2)]]
        tmp <- id.cov.data[which(id.cov.data$TRT==treat),names(id.cov.data)]
        names(tmp) <- toupper(names(tmp))
        out.data[[as.character(id)]][[as.character(treat)]]$cov <- tmp
        out.data[[as.character(id)]][[as.character(treat)]]$trt <- treat
        
      }
    }
    
    # Convert dosing information to scarabee-readable format
    for (id in uniqueIDs) {
      id.data <- dose.data[which(dose.data$ID==id),]
      id.data <- id.data[,-1]
      
      if (size(id.data,1) > 0) {  # There is dosing information for id
        
        for (treat in uniqueTRTs) {
          trt.data <- id.data[which(id.data$TRT==treat),]
          
          if (size(trt.data,1) > 0) {  # There is dosing information for treat
            bol.data <- trt.data[which(trt.data$RATE==0),
                                 c('TIME','CMT','AMT','RATE')]
            inf.data <- trt.data[which(trt.data$RATE!=0),-1]
            
            bol.data <- bol.data[order(bol.data$TIME,bol.data$CMT),]
            inf.data <- inf.data[order(inf.data$TIME,inf.data$CMT),]
            
            if (size(bol.data,1) > 0){
              out.data[[as.character(id)]][[as.character(treat)]]$bolus <-
                bol.data
            } else {
              out.data[[as.character(id)]][[as.character(treat)]]$bolus <-
                data.frame(TIME=numeric(0),
                           CMT=numeric(0),
                           AMT=numeric(0),
                           RATE=numeric(0))
            }
            if (size(inf.data,1) > 0){
              out.data[[as.character(id)]][[as.character(treat)]]$infusion <-
                inf.data
            } else {
              out.data[[as.character(id)]][[as.character(treat)]]$infusion <-
                data.frame(TIME=numeric(0),
                           CMT=numeric(0),
                           AMT=numeric(0),
                           RATE=numeric(0))
            }
            
          } else {                     # There is no dosing information for treat
            out.data[[as.character(id)]][[as.character(treat)]]$bolus <-
              data.frame(TIME=numeric(0),
                         CMT=numeric(0),
                         AMT=numeric(0),
                         RATE=numeric(0))
            out.data[[as.character(id)]][[as.character(treat)]]$infusion <-
              data.frame(TIME=numeric(0),
                         CMT=numeric(0),
                         AMT=numeric(0),
                         RATE=numeric(0))
          }
        }
        
      } else {                           # There is no dosing information for id
        for (treat in uniqueTRTs) {
          out.data[[as.character(id)]][[as.character(treat)]]$bolus <-
            data.frame(TIME=numeric(0),
                       CMT=numeric(0),
                       AMT=numeric(0),
                       RATE=numeric(0))
          out.data[[as.character(id)]][[as.character(treat)]]$infusion <-
            data.frame(TIME=numeric(0),
                       CMT=numeric(0),
                       AMT=numeric(0),
                       RATE=numeric(0))
        }
      }
    }
  
  ###### Start of population-based data processing
  } else {
    # Sort data by treatment, DVID, and time
    data <- data[order(data$TRT,
                       data$DVID,
                       data$TIME),]
    
    # Create backbone of function output list (with one slot per id)
    out.data <- list(list())
    names(out.data) <- '1'
    out.data <- c(out.data, list(ids=1))
    
    idTRTs <- unique(data$TRT)
    out.data[[1]] <- sapply(1:length(idTRTs),
                            function(x,tmp) c(tmp,list()),
                            tmp=list())
    names(out.data[[1]]) <- idTRTs
    out.data[[1]] <- c(out.data[[1]],list(trts=idTRTs))
    out.data[[1]]$id <- 1
    
    ## Extract dosing information
    dose.data <- data[which(data$EVID==1),]
    data <- data[which(data$EVID!=1),]
    if (size(data,1)>1) {
      row.names(data) <- 1:size(data,1)
    } else {
      stop(paste('the dataset does not contain any data ',
                 'for the dependent variable(s).',sep=''))
    }
    dose.data <- dose.data[,c(itt.index,dose.index,cov.index)]
    
    # Exclude sample records with MDV not set to 0
    if (length(which(data$MDV==0)) > 0) {
      data <- data[which(data$MDV==0),]
    } else {
      stop('no sample record with MDV=0 was found.')
    }
    
    # Split data between analysis data and covariates data by ID and TRT
    ana.data <- data[,c(itt.index,data.index)]
    cov.data <- data[,c(itt.index,cmt.index,cov.index)]
    
    uniqueTRTs <- unique(ana.data$TRT)
    
    for (treat in uniqueTRTs){
      if (length(which(ana.data$TRT==treat)) < 2 ) {
        stop(paste('there was less than two samples for treatment ',
                   treat, '.',sep=''))
      }
      out.data[[1]][[as.character(treat)]]$ana <-
        ana.data[which(ana.data$TRT==treat), names(ana.data)[-c(1,2)]]
      tmp <- cov.data[which(cov.data$TRT==treat), names(cov.data)]
      names(tmp) <- toupper(names(tmp))
      out.data[[1]][[as.character(treat)]]$cov <- tmp
      out.data[[1]][[as.character(treat)]]$trt <- treat
    }
    
    ## Convert dosing information to scarabee-readable format
    dose.data <- dose.data[,c('ID','TRT','TIME','CMT','AMT','RATE')]
    
    id.data <- dose.data[which(dose.data$ID==1),]
    id.data <- id.data[,-1]
    
    if (size(id.data,1) > 0) {  # There is dosing information
      
      for (treat in uniqueTRTs) {
        trt.data <- id.data[which(id.data$TRT==treat),]
        
        if (size(trt.data,1) > 0) {  # There is dosing information for treat
          bol.data <- trt.data[which(trt.data$RATE==0),
                               c('TIME','CMT','AMT','RATE')]
          inf.data <- trt.data[which(trt.data$RATE!=0),-1]
          
          bol.data <- bol.data[order(bol.data$TIME,bol.data$CMT),]
          inf.data <- inf.data[order(inf.data$TIME,inf.data$CMT),]
          
          if (size(bol.data,1) > 0){
            out.data[[1]][[as.character(treat)]]$bolus <- bol.data
          } else {
            out.data[[1]][[as.character(treat)]]$bolus <-
              data.frame(TIME=numeric(0),
                         CMT=numeric(0),
                         AMT=numeric(0),
                         RATE=numeric(0))
          }
          if (size(inf.data,1) > 0){
            out.data[[1]][[as.character(treat)]]$infusion <- inf.data
          } else {
            out.data[[1]][[as.character(treat)]]$infusion <-
              data.frame(TIME=numeric(0),
                         CMT=numeric(0),
                         AMT=numeric(0),
                         RATE=numeric(0))
          }
          
        } else {                     # There is no dosing information for treat
          out.data[[1]][[as.character(treat)]]$bolus <-
            data.frame(TIME=numeric(0),
                       CMT=numeric(0),
                       AMT=numeric(0),
                       RATE=numeric(0))
          out.data[[1]][[as.character(treat)]]$infusion <-
            data.frame(TIME=numeric(0),
                       CMT=numeric(0),
                       AMT=numeric(0),
                       RATE=numeric(0))
        }
      }
      
    } else {                           # There is no dosing information
      for (treat in uniqueTRTs) {
        out.data[[1]][[as.character(treat)]]$bolus <-
          data.frame(TIME=numeric(0),
                     CMT=numeric(0),
                     AMT=numeric(0),
                     RATE=numeric(0))
        out.data[[1]][[as.character(treat)]]$infusion <-
          data.frame(TIME=numeric(0),
                     CMT=numeric(0),
                     AMT=numeric(0),
                     RATE=numeric(0))
      }
    }
  }
  
  cat('  Done\n\n')
  
  return(list(data=out.data,new=new.data.file))
  
}
