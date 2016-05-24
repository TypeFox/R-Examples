
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

make.dosing <- function(allparms=NULL,
                        bolus=NULL,
                        infusion=NULL,
                        check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(bolus) | is.null(infusion))
      stop('bolus or infusion argument is NULL.')
    
    # Check infusion and bolus
    if (any(infusion$RATE==-1)){
      # Check for expected parameters
      cmts <- unique(infusion[which(infusion$RATE==-1),'CMT'])
      expected <- paste('R',cmts,sep='')
      missings <- expected[which(!expected%in%names(allparms))]
      if (length(missings)!=0){
        stop(paste('parameter(s) expected from dosing history not defined:\n  ',
                   paste(missings,collapse=', '),sep=''))
      }
      for (cmt in cmts){
        rcmt <- paste('allparms$R',cmt,sep='')
        # Check dimension of expected parameter
        if (size(eval(parse(text=rcmt)),1)>1){
          stop(paste(rcmt,' variable has invalid dimensions.',sep=''))
        }
        if (length(eval(parse(text=rcmt)))%in%c(1,size(infusion,2))){
          stop(paste(rcmt,' variable has invalid dimensions.',sep=''))
        }
      }
    } else if (any(infusion$RATE==-2)) {
      # Check for expected parameters
      cmts <- unique(infusion[which(infusion$RATE==-2),'CMT'])
      expected <- paste('D',cmts,sep='')
      missings <- expected[which(!expected%in%names(allparms))]
      if (length(missings)!=0){
        stop(paste('parameter(s) expected from dosing history not defined:\n  ',
                   paste(missings,collapse=', '),sep=''))
      }
      
      for (cmt in cmts){
        dcmt <- paste('allparms$D',cmt,sep='')
        # Check dimension of expected parameter
        if (size(eval(parse(text=dcmt)),1)>1){
          stop(paste(dcmt,' variable has invalid dimensions.',sep=''))
        }
        if (!(length(eval(parse(text=dcmt)))%in%c(1,size(infusion,2)))){
          stop(paste(dcmt,' variable has invalid dimensions.',sep=''))
        }
      }
    }
  }
  
  # Get estimated rate
  if (any(infusion$RATE==-1)){
    cmts <- unique(infusion[which(infusion$RATE==-1),'CMT'])
    
    for (cmt in cmts){
      rcmt <- paste('allparms$R',cmt,sep='')
      index <- which(infusion$CMT==cmt)
      if (length(eval(parse(text=rcmt)))==1){
        infusion$RATE[index] <- eval(parse(text=rcmt))
      } else {
        infusion$RATE[index] <- eval(parse(text=rcmt))[index]
      }
    }
  } else if (any(infusion$RATE==-2)) {
    cmts <- unique(infusion[which(infusion$RATE==-2),'CMT'])
    
    for (cmt in cmts){
      dcmt <- paste('allparms$D',cmt,sep='')
      index <- which(infusion$CMT==cmt)
      if (length(eval(parse(text=dcmt)))==1){
        infusion$RATE[index] <- infusion$AMT[index]/eval(parse(text=dcmt))
      } else {
        infusion$RATE[index] <- 
          infusion$AMT[index]/eval(parse(text=dcmt))[index]
      }
    }
  }
  
  # Subset infusion to first 4 columns (i.e, remove covariates)
  infusion <- infusion[,1:4]
  
  # Create dosing
  dosing <- data.frame(TIME=numeric(0),
                       CMT=numeric(0),
                       AMT=numeric(0),
                       RATE=numeric(0),
                       TYPE=numeric(0))
  
  cmts <- unique(c(unique(bolus$CMT),unique(infusion$CMT)))
  
  for (cmt in cmts){
    cmt.bolus <- bolus[bolus$CMT==cmt,]
    cmt.infusion <- infusion[infusion$CMT==cmt,]
    if (size(cmt.infusion,1) > 0){  # There is infusion data
      tmp.data <- convert.infusion(cmt.infusion)
      tmp.data$TYPE <- 0
      if(size(cmt.bolus,1) > 0){     # There is bolus data (need to infer rate)
        cmt.bolus$RATE <- approx(x=tmp.data$TIME,
                                 y=tmp.data$RATE,
                                 xout=cmt.bolus$TIME,
                                 method='constant',
                                 yleft=0,
                                 rule=2,
                                 f=0,
                                 ties='ordered')$y
        cmt.bolus$TYPE <- 1
        tmp.data <- rbind(tmp.data, cmt.bolus)
        tmp.data <- tmp.data[order(tmp.data$TIME,
                                   tmp.data$CMT,
                                   tmp.data$TYPE),]
      }
      tmp.data <- tmp.data[,c('TIME','CMT','AMT','RATE','TYPE')]
    } else {                    # There is no infusion data
      tmp.data <- cmt.bolus[,c('TIME','CMT','AMT','RATE')]
      tmp.data$TYPE <- 1
    }
    dosing <- rbind(dosing,tmp.data)
  }
  
  # Add lag-times on dosing events if necessary
  abslags <- allparms[grep('^ALAG[[:digit:]]*$',names(allparms))]
  if (length(abslags)>0){
    cmtlags <- gsub('ALAG','',names(abslags))
    for (cmtlag in cmtlags){
      index <- which(dosing$CMT==cmtlag)
      dosing$TIME[index] <- dosing$TIME[index] + 
        eval(parse(text=paste('allparms$ALAG',cmtlag,sep='')))
    }
  }
  
  # Re-order dosing
  dosing <- dosing[order(dosing$TIME,
                         dosing$CMT,
                         dosing$TYPE),]
  
  return(as.matrix(dosing))
}
