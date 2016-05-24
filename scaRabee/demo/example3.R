require(scaRabee)

if (.Platform$OS.type=="windows"){
  end.of.folder <- ''
} else {
  end.of.folder <- '/'
}

# User-prompt: define target directory
if (interactive()){
  cat('\nExample 3 - Simulation of a model defined with ordinary differential equations\n')
  cat('at the population level\n\n')
  
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
wd <- paste(wd,'/example3/',sep='')
ana.file <- 'example3.R' #paste(wd,'example3.R',sep='')
data.file <- 'data.csv' # paste(wd,'data.csv',sep='')
param.file <- 'initials.csv' #paste(wd,'initials.csv',sep='')
model.file <- 'model.txt' #paste(wd,'model.txt',sep='')

# Copy files
data(example3.data,
  example3.initials)

# Create main and model.definition directory
if (file.exists(wd)) {
  stop(sprintf(paste('\nDirectory \'%s\' already exists.\nDemo aborted. ',
        'Please retry using a different target directory.\n',sep=''),wd),
    call.=FALSE)
}
scarabee.new(name='example3',
  path=wd,
  type='simulation',
  method='population',
  template='ode')
setwd(wd)

# Update the data file
write.table(example3.data,
  file=data.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  append=FALSE)

# Update the parameter file
write.table(example3.initials,
  file=param.file,
  sep=',',
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE,
  append=TRUE)

# Update the model file
tmp <- c('$ANALYSIS example3',
  '',  
  '$DERIVED',
  '  if (DOSE[1]==1){',
  '    F <- F1',
  '  } else if (DOSE[1]==3){',
  '    F <- F2',
  '  } else {',
  '    F <- F3',
  '  }',
  '  if (SC[1]==1) {',
  '    F <- F/(20000*380)',
  '  } else {',
  '    F <- 1/(20000*380)',
  '  }',
  '',
  '$IC',
  '  init <- c(0,0,0,0,0)',
  '',
  '$SCALE',
  '  scale <- c(1/F,1,1/F,1,1)',
  '',
  '$ODE',
  '  Rf <- Rmax - a5',
  '  dadt <- rbind(-ka*a1,                                              # SC',
  '                ka*a1-ka2*a2,                                        # AL',
  '                ka2*a2+ktp*a4+koff*a5-(kon/Vc)*a3*Rf-(kpt+kloss)*a3, # AP',
  '                kpt*a3-ktp*a4,                                       # AT',
  '                (kon/Vc)*a3*Rf-(koff+kint)*a5)                       # DR',
  '',
  '$OUTPUT',
  '  y <- rbind(f[3,]*(20000*380)/Vc)')

write(tmp,
  file=model.file,
  sep='\n',
  append=FALSE)

# Run example 3
source(ana.file)
setwd(old.wd)
