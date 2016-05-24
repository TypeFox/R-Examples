
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

scarabee.new <- function(name = 'myanalysis',
                         path = NULL,
                         type = 'simulation',
                         method = 'population',
                         template = 'ode'){
  
  # Check inputs
  if (!is.null(name)){
    if (!is.vector(name,'character')){
      warning(paste('name argument was not a character vector and was ',
                    'coerced to \'myanalysis\'',sep=''))
      name <- 'myanalysis'
    }
    
    if (length(name)!=1){
      warning('name argument coerced to its first element')
      name <- name[1]
    }
    
    name <- as.character(name)
    
  }
  
  if (is.null(path)) path <- paste(getwd(),name,sep='/')
  
  if (!(type=='estimation' | type=='simulation' | type=='gridsearch')){
    warning(paste('type argument was neither \'estimation\', ',
                  '\'simulation\', nor \'gridsearch\'. It was\n  coerced to ',
                  '\'simulation\'',sep=''))
    type <- 'simulation'
  }
  
  if (!(method=='subject' | method=='population')){
    warning(paste('method argument was neither \'population\' nor ',
                  '\'subject\'. It was coerced to \'population\'',sep=''))
    method <- 'population'
  }
  
  if (!is.null(template)){
    if (!is.vector(template,'character')){
      warning(paste('template argument was not a character vector and was ',
                    'coerced to \'ode\'',sep=''))
      template <- 'ode'
    }
    
    if (length(template)!=1){
      warning('template argument coerced to its first element')
      template <- template[1]
    }
    
    if (!is.element(template,c('explicit','ode','dde'))){
      warning(paste('template argument was not a valid selection and was ',
                    'coerced to \'ode\'.',sep=''))
      template <- 'ode'
    }
  }
  
  # Files names
  if (substring(path,nchar(path),nchar(path))!='/'){
    path <- paste(path,'/',sep='')
  }
  
  ana.file <- paste(path,name,'.R',sep='')
  data.file <- paste(path,'data.csv',sep='')
  param.file <- paste(path,'initials.csv',sep='')
  model.file <- paste(path,'model.txt',sep='')
  
  # Create directory
  if (!file.exists(path)){
    
    tmp <- try(dir.create(path,showWarnings=FALSE,recursive=TRUE),silent=TRUE)
    
    if (any(!tmp))
      stop(paste('new scaRabee working could not be created. Please, check path',
                 'argument or permission\n  to target directory.'))
  } else {
    stop('path argument must be different because this directory already exists.')
  }
  
  # Create main analysis script
  tmp <- sprintf(
    paste(
      '#Copyright (c) 2009-2014 Sebastien Bihorel',
      '#All rights reserved.',
      '#',
      '#This file is part of scaRabee.',
      '#',
      '#    scaRabee is free software: you can redistribute it and/or modify',
      '#    it under the terms of the GNU General Public License as published by',
      '#    the Free Software Foundation, either version 3 of the License, or',
      '#    (at your option) any later version.',
      '#',
      '#    scaRabee is distributed in the hope that it will be useful,',
      '#    but WITHOUT ANY WARRANTY; without even the implied warranty of',
      '#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the',
      '#    GNU General Public License for more details.',
      '#',
      '#    You should have received a copy of the GNU General Public License',
      '#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.',
      '#',
      '################################################################################',
      '#',
      '# User settings:',
      '#',
      '# (Optional) Please provide a working directory if working in a interactive R',
      '# session',
      '# setwd()',
      '#',
      '# Please provide the names of the following files',
      '#   - data:  the dataset; must be a .csv file',
      '#   - param: the initial guess for model and residual variability',
      '#            parameters; must be a .csv file with 6 columns',
      '#   - model: the model; a .R file located in the model.definition',
      '#            subdirectory',
      '',
      'files <- list(data  = \'data.csv\',',
      '              param = \'initials.csv\',',
      '              model = \'model.txt\')',
      '',
      '# Please, enter the type of run: \'simulation\', \'estimation\', or \'gridsearch\'',
      '  runtype <- \'%s\'',
      '',
      '# Please, define the analysis method: \'subject\' or \'population\'',
      '  method <- \'%s\'',
      '',
      '# Please, indicate if you want to run the analysis in debug mode: TRUE or FALSE',
      '  debugmode <- TRUE',
      '',
      '# Please, enter the maximum number of iterations or function evaluations',
      '# for the estimation',
      '  estim <- list(maxiter=500,',
      '                maxfunc=5000)',
      '',
      '# Please, enter the number of points per grid dimension and the dispersion',
      '# factor',
      '  npts <- 3',
      '  alpha <- NULL',
      '',
      '################################################################################',
      'require(scaRabee)',
      '',
      'scarabee.analysis(files=files,',
      '                  runtype=runtype,',
      '                  method=method,',
      '                  debugmode=debugmode,',
      '                  estim.options=estim,',
      '                  npts=npts,',
      '                  alpha=alpha)',
      '',
      'warnings()',
      'traceback()',sep='\n'),
    type, method)
    
  write.table(tmp,
              file=ana.file,
              sep='\n',
              quote=FALSE,
              row.names=FALSE,
              col.names=FALSE)
  
   # Create data file
  write('OMIT,ID,TRT,TIME,AMT,RATE,CMT,EVID,DV,DVID,MDV',
        file=data.file,
        sep='\n')
  
  # Create initials file
  write('Parameter,Type,Value,Fixed,Lower bound,Upper bound',
        file=param.file,
        sep='\n')
  
  # Create model file
  if (template == 'explicit'){
    tmp <- sprintf(paste(
        '$ANALYSIS %s',
        '',
        '$OUTPUT',
        '  y <- rbind()',
        '',
        '$VARIANCE',
        '  v <- rbind()',
        '',
        '$SECONDARY',
        sep='\n'),
      name)
  } else if (template == 'ode'){
    tmp <- sprintf(paste(
        '$ANALYSIS %s',
        '',
        '$DERIVED',
        '',
        '$IC',
        '  init <- c()',
        '',
        '$SCALE',
        '  scale <- c()',
        '',
        '$ODE',
        '  dadt <- rbind()',
        '',
        '$OUTPUT',
        '  y <- rbind()',
        '',
        '$VARIANCE',
        '  v <- rbind()',
        '',
        '$SECONDARY',
        sep='\n'),
      name)
  } else if (template == 'dde'){
    tmp <- sprintf(paste(
        '$ANALYSIS %s',
        '',
        '$DERIVED',
        '',
        '$IC',
        '  init <- c()',
        '',
        '$SCALE',
        '  scale <- c()',
        '',
        '$LAGS',
        '',
        '$DDE',
        '  dadt <- rbind()',
        '',
        '$OUTPUT',
        '  y <- rbind()',
        '',
        '$VARIANCE',
        '  v <- rbind()',
        '',
        '$SECONDARY',
        sep='\n'),
      name)
  }
  
  write.table(tmp,
              file=model.file,
              sep='\n',
              quote=FALSE,
              row.names=FALSE,
              col.names=FALSE)
  
  setwd(path)
            
  # Display message
  cat(sprintf('\nA new scaRabee directory has been created at:\n%s\n',
              path))
  
  cat(sprintf('\nWorking directory set to:\n%s\n\n',
              path))
  
}