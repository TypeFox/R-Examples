require(scaRabee)

if (.Platform$OS.type=="windows"){
  end.of.folder <- ''
} else {
  end.of.folder <- '/'
}

# User-prompt: define target directory
if (interactive()){
  cat('\nExample 2 - Simulation of inputs into a model defined with differential\n')
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
wd <- paste(wd,'/example2/',sep='')
ana.file <- 'example2.R' #paste(wd,'example2.R',sep='')
data.file <- 'data.csv' # paste(wd,'data.csv',sep='')
param.file <- 'initials.csv' #paste(wd,'initials.csv',sep='')
model.file <- 'model.txt' #paste(wd,'model.txt',sep='')

# Copy files
data(example2.data,
  example2.initials)

# Create main and model.definition directory
if (file.exists(wd)) {
  stop(sprintf(paste('\nDirectory \'%s\' already exists.\nDemo aborted. ',
        'Please retry using a different target directory.\n',sep=''),wd),
    call.=FALSE)
}
scarabee.new(name='example2',
  path=wd,
  type='simulation',
  method='population',
  template='ode')


# Update the data file
write.table(example2.data,
  file=data.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  append=FALSE)

# Update the parameter file
write.table(example2.initials,
  file=param.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE,
  append=TRUE)

# Update the model file
tmp <- c('$ANALYSIS example2',
  '',
  '$DERIVED',
  '',
  '$IC',
  '  init <- c(IC1,IC2)',
  '',
  '$SCALE',
  '  scale <- c(1,2)',
  '',
  '$ODE',
  '  dadt <- rbind(0,0)',
  '',
  '$OUTPUT',
  '  y <- f')

write(tmp,
  file=model.file,
  sep='\n',
  append=FALSE)

# Run example 2
source(ana.file)
setwd(old.wd)
