
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

convert.infusion <- function(infusion.data=NULL){
  
  # Check inputs
  if (!is.data.frame(infusion.data)){
    stop('infusion.data argument is not a data frame.')
  }
  
  if (!all(names(infusion.data)==c('TIME','CMT','AMT','RATE'))){
    stop(paste('variables in infusion.data are missing or not',
               'in the proper order (TIME, CMT,\n  AMT, RATE).',sep=' '))
  }
  
  # Check uniqueness of CMT
  if (length(unique(infusion.data$CMT))>1) {
    stop('the CMT variable must contain identical values.')
  }
  
  # Get number of infusion events
  inf.nb <- size(infusion.data,1)
  
  if (inf.nb==0){
    return(infusion.data)
  }
  
  # Calculate time of end of infusion
  infusion.data$TEND <- infusion.data$TIME+infusion.data$AMT/infusion.data$RATE
  
  # Create processing matrix
  inf.mat <- matrix(0,nrow=2*inf.nb,ncol=1+inf.nb)
  
  # Set start and end times of infusions 
  inf.mat[2*(1:inf.nb)-1] <-infusion.data$TIME
  inf.mat[2*(1:inf.nb)] <-infusion.data$TEND
  inf.mat <- inf.mat[order(inf.mat[,1]),]
  
  # Set rate from start of infusion time to end of infusion
  for (inf in 1:inf.nb) {
    inf.mat[which((inf.mat[,1]>=infusion.data$TIME[inf]) &
                  (inf.mat[,1]<infusion.data$TEND[inf])),inf+1] <- 
                     infusion.data$RATE[inf]
  }
  
  # Create copy of the processing matrix
  copy.mat <- matrix(NA,nrow=2*inf.nb,ncol=1+inf.nb)
  
  # Set start and end times of infusions 
  copy.mat[2*(1:inf.nb)-1] <-infusion.data$TIME
  copy.mat[2*(1:inf.nb)] <-infusion.data$TEND
  
  
  # Set rate to RATE at start of infusion time and to 0 at end of infusion
  for (inf in 1:inf.nb) {
    copy.mat[2*inf-1,inf+1] <- infusion.data$RATE[inf]
    copy.mat[2*inf,inf+1] <- 0
  }
  copy.mat <- copy.mat[order(copy.mat[,1]),]
  
  # Set rate to 0 before start of infusions
  for (inf in 1:inf.nb) {
    tmp <- which(!is.na(copy.mat[,inf+1]))[1]
    if (tmp > 1) copy.mat[1:(tmp-1),inf+1] <- 0
  }
  
  # Consolidate inf.mat and copy.mat
  if (length(which(!is.na(copy.mat)))!=0) {
    inf.mat[which(!is.na(copy.mat))] <- copy.mat[which(!is.na(copy.mat))]
  }
  
  # Calculate the overall infusion rate by summing all rates
  if (inf.nb==1) {
    processed.infusion <- data.frame(TIME=inf.mat[,1],
                                     RATE=inf.mat[,2])
  } else {
    processed.infusion <- data.frame(TIME=inf.mat[,1],
                                     RATE=apply(inf.mat[,2:(inf.nb+1)],1,sum))
  }
  
  # Add variables
  processed.infusion$TRT <- infusion.data$TRT[1]
  processed.infusion$CMT <- infusion.data$CMT[1]
  processed.infusion$AMT <- 0
  
  # Return the processed infusion data frame
  return(processed.infusion[,c('TIME','CMT','AMT','RATE')])
  
}
