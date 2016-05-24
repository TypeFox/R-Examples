require(scaRabee)

if (.Platform$OS.type=="windows"){
  end.of.folder <- ''
} else {
  end.of.folder <- '/'
}

# User-prompt: define target directory
if (interactive()){
  cat('\nExample 8 - Direct grid search for a model defined with delay differential\n')
  cat('equations\n\n')
  
  repeat{
    wd <- readline('Enter a path to store the demo files:\n>')
    if (wd!=''){
      if (substring(wd,nchar(wd),nchar(wd))!='/'){
        wd <- paste(wd,end.of.folder,sep='')
      }
      if (!file.exists(wd)){
        action <- readline(sprintf(paste('\nDirectory \'%s\' does not exist:\n',
              '  [c] Continue with current working directory: %s\n',
              '  [r] Retry\n',
              '  [a] Abort\n>',sep=''),wd,getwd()))
        if (action == 'a') {
          stop('Demo aborted',call.=FALSE)
        } else if (action == 'c') {
          wd <- getwd()
          if (substring(wd,nchar(wd),nchar(wd))!='/'){
            wd <- paste(wd,end.of.folder,sep='')
          }
          options(warn=-1)
          nd <- try(file.create(paste(wd,'test.R',sep='')))
          options(warn=0)
          
          if (nd) { # User has permission on directory
            file.remove(paste(wd,'test.R',sep=''))
            break
          }
          cat('\nYou don\'t have permissions on this directory.\n')
        }
        
      } else {
        if (substring(wd,nchar(wd),nchar(wd))!='/'){
          wd <- paste(wd,end.of.folder,sep='')
        }
        options(warn=-1)
        nd <- try(file.create(paste(wd,'test.R',sep='')))
        options(warn=0)
        
        if (nd) { # User has permission on directory
          file.remove(paste(wd,'test.R',sep=''))
          break
        }
        cat('\nYou don\'t have permissions on this directory.\n')
      }
      
    } else {
      
      wd <- getwd()
      if (substring(wd,nchar(wd),nchar(wd))!='/'){
        wd <- paste(wd,end.of.folder,sep='')
      }
      options(warn=-1)
      if (file.exists(wd))
        nd <- try(file.create(paste(wd,'test.R',sep='')))
      options(warn=0)
      
      if (nd) { # User has permission on directory
        file.remove(paste(wd,'test.R',sep=''))
        break
      }
      
      cat('\nYou don\'t have permissions on this directory.\n')
      
    }
  }
} else {
  return(NULL)
}

# Set files
old.wd <- getwd()
wd <- paste(wd,'/example8/',sep='')
ana.file <- paste(wd,'example8.R',sep='')
data.file <- paste(wd,'data.csv',sep='')
param.file <- paste(wd,'initials.csv',sep='')
model.file <- paste(wd,'model.txt',sep='')

# Copy files
data(example8.data,
  example8.initials)

# Create main and model.definition directory
if (file.exists(wd)) {
  stop(sprintf(paste('\nDirectory \'%s\' already exists.\nDemo aborted. ',
        'Please retry using a different target directory.\n',sep=''),wd),
    call.=FALSE)
}
scarabee.new(name='example8',
  path=wd,
  type='gridsearch',
  method='population',
  template='dde')
setwd(wd)

# Update the data file
write.table(example8.data,
  file=data.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  append=FALSE)

# Update the parameter file
write.table(example8.initials,
  file=param.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE,
  append=TRUE)

# Update the model file
tmp <- c('$ANALYSIS example8',
  '',
  '$DERIVED',
  '  kin <- 100/Tr',
  '$IC',
  '  init <- c(0,100)',
  '',
  '$SCALE',
  '  scale <- c(1,1)',
  '',
  '$LAGS',
  '  lag1 <- Tm',
  '  lag2 <- Tp+Tm',
  '  lag3 <- Tm+Tr',
  '  lag4 <- Tp+Tm+Tr',
  '$DDE',
  '  # PK explicit',
  '  flag <- ifelse((t-3)>0,1,0)',
  '',
  '  Cp <- ratio*(L1/lbd1)*(1-flag*(1-exp(-lbd1*(t-3)))-exp(-lbd1*t)) +',
  '        ratio*(L2/lbd2)*(1-flag*(1-exp(-lbd2*(t-3)))-exp(-lbd2*t)) +',
  '        ratio*(L3/lbd3)*(1-flag*(1-exp(-lbd3*(t-3)))-exp(-lbd3*t))',
  '',
  '  # PD',
  '  dadt <- rbind(kmax*Cp/(kc50+Cp),',
  '                kin*exp(-(alag.lag1[1]-alag.lag2[1]))-',
  '                  kin*exp(-(alag.lag3[1]-alag.lag4[1])))',
  '',
  '$OUTPUT',
  '  # PK explicit',
  '  flag <- ifelse((times-3)>0,1,0)',
  '',
  '  Cp <- ratio*(L1/lbd1)*(1-flag*(1-exp(-lbd1*(times-3)))-exp(-lbd1*times)) +',
  '        ratio*(L2/lbd2)*(1-flag*(1-exp(-lbd2*(times-3)))-exp(-lbd2*times)) +',
  '        ratio*(L3/lbd3)*(1-flag*(1-exp(-lbd3*(times-3)))-exp(-lbd3*times))',
  '',
  '  y <- rbind(f[2,])',
  '',
  '$VARIANCE',
  '  v <- rbind((SD1+CV1*y[1,])^2)'
)

write(tmp,
  file=model.file,
  sep='\n',
  append=FALSE)

# Update the master script
tmp <- sprintf(
  paste(
    '#Copyright (c) 2009-2011 Sebastien Bihorel',
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
    '  alpha <- 2',
    '',
    '################################################################################',
    'require(scaRabee)',
    '',
    paste('scarabee.analysis(files=files,',
      '                  runtype=runtype,',
      '                  method=method,',
      '                  debugmode=debugmode,',
      '                  estim.options=estim,',
      '                  npts=npts,',
      '                  alpha=alpha,',
      '                  dde.options=dde.options)',
      sep='\n'),
    '',
    'warnings()',
    'traceback()',sep='\n'),
  'gridsearch', 'population')

write(tmp,
  file=ana.file,
  sep='\n',
  append=FALSE)

# Run example 8
source(ana.file)
setwd(old.wd)
