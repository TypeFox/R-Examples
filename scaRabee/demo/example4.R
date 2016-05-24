require(scaRabee)

if (.Platform$OS.type=="windows"){
  end.of.folder <- ''
} else {
  end.of.folder <- '/'
}

# User-prompt: define target directory
if (interactive()){
  cat('\nExample 4 - Simulation of a model defined with delay differential equations \n')
  cat('the population level\n\n')
  
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
wd <- paste(wd,'/example4/',sep='')
ana.file <- paste(wd,'example4.R',sep='')
data.file <- paste(wd,'data.csv',sep='')
param.file <- paste(wd,'initials.csv',sep='')
model.file <- paste(wd,'model.txt',sep='')

# Copy files
data(example4.data,
  example4.initials)

# Create main and model.definition directory
if (file.exists(wd)) {
  stop(sprintf(paste('\nDirectory \'%s\' already exists.\nDemo aborted. ',
        'Please retry using a different target directory.\n',sep=''),wd),
    call.=FALSE)
}
scarabee.new(name='example4',
  path=wd,
  type='simulation',
  method='population',
  template='dde')
setwd(wd)

# Update the data file
write.table(example4.data,
  file=data.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  append=FALSE)

# Update the parameter file
write.table(example4.initials,
  file=param.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE,
  append=TRUE)

# Update the model file
tmp <- c('$ANALYSIS example4',
  '',
  '$DERIVED',
  '  ke <- CL/V1',
  '  k12 <- Q/V1',
  '  k21 <- Q/V2',
  '',
  '$IC',
  '  init <- c(IC1,IC2)',
  '',
  '$SCALE',
  '  scale <- c(1,1)',
  '',
  '$LAGS',
  '',
  '$DDE',
  '  dadt <- rbind(-(ke+k12)*a1+k21*a2,',
  '                k12*alag.XYZ[1]-k21*a2)',
  '',
  '$OUTPUT',
  '  y <- rbind(f[1,]/V1)',
  '',
  '$VARIANCE',
  '  v <- rbind((CV*f[1,])^2))',
  '',
  '$SECONDARY')

write(tmp,
  file=model.file,
  sep='\n',
  append=FALSE)

# Run example 4
source(ana.file)
setwd(old.wd)
