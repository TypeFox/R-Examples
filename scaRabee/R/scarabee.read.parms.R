
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

scarabee.read.parms <- function(files=NULL){
  
  # Check inputs
  if (is.null(files)){
    stop('files argument is NULL.')
  }
  
  if (is.null(files$param)){
    stop('files argument does not have any param level or files$param is NULL.')
  }
  
  if (!file.exists(files$param)){
    stop('parameter file does not exist.')
  }
  
  # Read parameter file
  cat('Processing parameter file:\n')
  param <- read.csv(files$param,
                    header=TRUE,
                    as.is=TRUE,
                    col.names=c('names','type','value','isfix','lb','ub'))
  
  # Check parameter
  if (size(param,1)==0)
    stop(sprintf(paste('the matrix of parameters created from %s',
                       'is empty. Please, check the content of your file.'),
                 files$param))
  
  if (class(param$names)!='character')
    param$names <- as.character(param$names)
  
  if (class(param$type)!='character')
    param$type <- as.character(param$type)
  
  if (class(param$value)!='numeric')
    param$value <- as.numeric(param$value)
  
  if (class(param$isfix)!='integer')
    param$isfix <- as.numeric(param$isfix)
  
  if (class(param$lb)!='numeric')
    param$lb <- as.numeric(param$lb)
  
  if (class(param$ub)!='numeric')
    param$ub <- as.numeric(param$ub)
  
  if (any(is.na(param$names)))
    stop('parameter names cannot be NA.')
  
  if (any(is.na(param$type)))
    stop('parameter type cannot be NA.')
  
  if (any(is.na(param$value)))
    stop('parameter value cannot be NA.')
  
  if (any(is.na(param$isfix)))
    stop('parameter isfix cannot be NA.')
  
  if (any(is.na(param$lb)))
    stop('parameter lb cannot be NA.')
  
  if (any(is.na(param$ub)))
    stop('parameter ub cannot be NA.')
  
  if (!all(param$type%in%c('P','IC','L','V'))){
    wrongtype <- which(!param$type%in%c('P','IC','L','V'))
    stop(paste('parameters with incorrect type; check index:\n  ',
               paste(wrongtype,collapse=', '),sep=''))
  }
  
  cat('  Done\n\n')
  
  return(param)
  
}
