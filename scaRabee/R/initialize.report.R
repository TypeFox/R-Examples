
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

initialize.report <- function(problem=NULL,
                              param=NULL,
                              files=NULL,
                              isgrid=0){
  
  # Creates and populates the report file
  if (!isgrid){
    tmp <- 'Maximum likelihood estimation report'
    tmp <- c(tmp,paste('\nEstimation started at: ',Sys.time(),'\n',sep=''))
  } else {
    tmp <- 'Direct grid search report'
    tmp <- c(tmp,paste('\nGrid search started at: ',Sys.time(),'\n',sep=''))
  }
  
  tmp <- c(tmp,paste('Model defined in:',
                     paste(rep(' ',24),collapse=''),
                     files$model,sep=''))
  tmp <- c(tmp,paste('Data imported from:',
                     paste(rep(' ',22),collapse=''),
                     files$data,
                     sep=''))
  tmp <- c(tmp,paste('Parameter characteristics imported from: ',
                     files$param,
                     sep=''))
  tmp <- c(tmp,paste('\nLevel of analysis:', problem$method))
  
  write(tmp,file=files$report,sep='\n')
  
  if (!isgrid){
    # Determines which parameters will be estimated
    x <- param[which(param$isfix==0),]
    fixparam <- param[which(param$isfix==1),]
    
    # Creates and populates the iteration log
    mystr <- paste('ID','Iteration,Objective function,Number of evaluations,Procedure',
                   paste(x$names,collapse=','),
                   'Time',
                   sep=',')
    
    write(mystr,file=files$iter)
    
    # Creates and populates the prediction file
    write('ID,TRT,NUM,DVID,TIME,DV,IPRED,IRES,W,WRES',
          file=files$pred,append=FALSE,sep='\n')
    
  }
}

