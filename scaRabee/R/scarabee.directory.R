
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

scarabee.directory <- function (curwd=getwd(),
                                files=NULL,
                                runtype=NULL,
                                analysis=NULL){
  
  # Check inputs
  if (!is.null(curwd)){
    if (!file.exists(curwd)){
      stop('curwd argument is expected to be an existing path.')
    }
  } else {
    stop('curwd is expected to be an existing path, but is NULL.')
  }
  
  # Set length of substring for directory name 
  if (runtype=='gridsearch'){
    lentype <- 4
  } else {
    lentype <- 3
  }
  
  # Interactive call
  if (interactive()){
    repeat{
      wd <- readline(sprintf(paste('Change the current working directory: %s?\n',
                  '(Enter a new path if you want to use a new ',
                  'directory)\n',
                  sep=''),
              curwd))
      if (wd!=''){
        if (file.exists(wd)){
          setwd(wd)
          break
        }
      } else {
        break
      }
    }
  }
  
  # Create new directory and set working directory
  dirCreated <- 0 ; i <- 1
  
  while (dirCreated==0){
    if (i<10){
      dirName <- sprintf('%s.%s.0%d',
          analysis,
          substring(runtype,1,lentype),
          i)
    } else {
      dirName <- sprintf('%s.%s.%d',
          analysis,
          substring(runtype,1,lentype),
          i)
    }
    if (file.exists(dirName)){
      i <- i+1
    } else {
      dir.create(dirName)
      dirCreated <- 1
    }
  }
  
  newwd <- paste(getwd(),'/',dirName,'/',sep='')
  setwd(newwd)
  
  # Create sub-directories
  dir.create(paste(newwd,'run.config.files/',sep=''))
   
  # Create files
  file.create(paste(newwd,files$data,sep=''))
  if (!is.null(files$newdata)) file.create(paste(newwd,files$newdata,sep=''))
  file.create(paste(newwd,files$param,sep=''))
  file.create(paste(newwd,'run.config.files/',files$model,sep=''))
  file.create(paste(newwd,'run.config.files/',files$data,sep=''))
  file.create(paste(newwd,'run.config.files/',files$param,sep=''))
  
  # Copy configuration files to working directory
  file.copy(from=paste('../',files$data,sep=''),
      to=paste(newwd,files$data,sep=''),overwrite=TRUE)
  if (!is.null(files$newdata)) {
    file.copy(from=paste('../',files$newdata,sep=''),
        to=paste(newwd,files$newdata,sep=''),overwrite=TRUE)
  }
  file.copy(from=paste('../',files$param,sep=''),
      to=paste(newwd,files$param,sep=''),overwrite=TRUE)
  file.copy(from=paste('..',files$model,sep=''),
      to=paste(newwd,files$model,sep=''),overwrite=TRUE)
  
  # Copy configuration files to backup directory
  file.copy(from=paste('../',files$data,sep=''),
      to=paste(newwd,'run.config.files/',files$data,sep=''),overwrite=TRUE)
  if (!is.null(files$newdata)) {
    file.copy(from=paste('../',files$newdata,sep=''),
      to=paste(newwd,'run.config.files/',files$newdata,sep=''),overwrite=TRUE)
  }
  file.copy(from=paste('../',files$param,sep=''),
      to=paste(newwd,'run.config.files/',files$param,sep=''),overwrite=TRUE)
  file.copy(from=paste('../',files$model,sep=''),
      to=paste(newwd,'run.config.files/',files$model,sep=''),overwrite=TRUE)
  
  # Copy original R script if possible
  script <- paste('../',analysis,'.R',sep='')
  
  if (file.exists(script)){
    file.create(paste(newwd,analysis,'.R',sep=''))
    file.create(paste(newwd,'run.config.files/',analysis,'.R',sep=''))
    file.copy(from=script,to=paste(newwd,analysis,'.R',sep=''),overwrite=TRUE)
    file.copy(from=script,
        to=paste(newwd,'run.config.files/',analysis,'.R',sep=''),
        overwrite=TRUE)
  }
  
}
