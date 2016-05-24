require(scaRabee)

if (.Platform$OS.type=="windows"){
  end.of.folder <- ''
} else {
  end.of.folder <- '/'
}

# User-prompt: define target directory
if (interactive()){
  cat('\nExample 1 - Simulation of a model defined with algebraic equations at the \n')
  cat('population level\n\n')
  
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
wd <- paste(wd,'/example1/',sep='')
ana.file <- 'example1.R'
data.file <- 'data.csv'
param.file <- 'initials.csv'
model.file <- 'model.txt'

# Copy files
data(example1.data,
  example1.initials)

# Create main and model.definition directory
if (file.exists(wd)) {
  stop(sprintf(paste('\nDirectory \'%s\' already exists.\nDemo aborted. ',
        'Please retry using a different target directory.\n',sep=''),wd),
    call.=FALSE)
}
scarabee.new(name='example1',
  path=wd,
  type='simulation',
  method='population',
  template='explicit')

# Update the data file
write.table(example1.data,
  file=data.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  append=FALSE)

# Update the parameter file
write.table(example1.initials,
  file=param.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE,
  append=TRUE)

# Update the model file
tmp <- c('$ANALYSIS example1',
  '',
  '$OUTPUT',
  '  ke <- cl/vc',
  '  flag <- ifelse((times-2)>0,1,0)',
  '',
  '  Cp <- (DOSE[1]*500/vc)*((1/ke)*(1-flag*(1-exp(-ke*(times-2)))-exp(-ke*times)))',
  '  y <- rbind(log(Cp),',
  '             base*(1-imax*Cp/(ic50+Cp)))')
write(tmp,
  file=model.file,
  sep='\n',
  append=FALSE)

# Run example 1
source(ana.file)
setwd(old.wd)
