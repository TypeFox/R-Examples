require(scaRabee)

if (.Platform$OS.type=="windows"){
  end.of.folder <- ''
} else {
  end.of.folder <- '/'
}

# User-prompt: define target directory
if (interactive()){
  cat('\nExample 6 - Simulation of a model defined with ordinary differential equations at the \n')
  cat('subject level\n\n')
  
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
wd <- paste(wd,'/example6/',sep='')
ana.file <- paste(wd,'example6.R',sep='')
data.file <- paste(wd,'data.csv',sep='')
param.file <- paste(wd,'initials.csv',sep='')
model.file <- paste(wd,'model.txt',sep='')

# Copy files
data(example6.data,
  example6.initials)

# Create main and model.definition directory
if (file.exists(wd)) {
  stop(sprintf(paste('\nDirectory \'%s\' already exists.\nDemo aborted. ',
        'Please retry using a different target directory.\n',sep=''),wd),
    call.=FALSE)
}
scarabee.new(name='example6',
  path=wd,
  type='simulation',
  method='subject',
  template='explicit')
setwd(wd)

# Update the data file
write.table(example6.data,
  file=data.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  append=FALSE)

# Update the parameter file
write.table(example6.initials,
  file=param.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE,
  append=TRUE)

# Update the model file
tmp <- c('$ANALYSIS example6',
  '',
  '$DERIVED',
  '  ke <- CL/V',
  '  kin <- kout*R0',
  '$IC',
  '  init <- c(0,R0,R0)',
  '',
  '$SCALE',
  '  scale <- c(1,1,1)',
  '',
  '$ODE',
  '  inh <- 1 - Imax*(a1/V)/(IC50+(a1/V))',
  '  dadt <- rbind(-ke*a1,                 # CP',
  '                kin-kt*inh*a2,          # P',
  '                kt*inh*a2 - kout*a3)    # R',
  '',
  '$OUTPUT',
  '  y <- rbind(f[1,]/V,',
  '             f[2,],',
  '             f[3,])'
)

write(tmp,
  file=model.file,
  sep='\n',
  append=FALSE)

# Run example 6
source(ana.file)
setwd(old.wd)
