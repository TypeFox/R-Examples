package.dir <- function(package='base', lib.loc=NULL,
      exclude=c('chtml', 'data', 'help', 'html', 'latex', 'libs',
                'man', 'Meta', 'po', 'R', 'R-ex', 'src'),
      include=NULL, pattern=NULL, recursive=FALSE){
##  
##  1.  fullPath <- system.file(package=package, lib.loc=lib.loc)
##
  fullPath <- system.file(package=package, lib.loc=lib.loc)
##  
##  2.  Dir <- dir(fullPath)
##
#   2.1.  Directory?    
  Dir <- dir(fullPath)
  Dir. <- dir(fullPath, full.names=TRUE)
#   2.2.  Identify subdirectories
  DirInfo <- sapply(Dir., function(x)file.info(x)$isdir)
  subDirs <- Dir[DirInfo]   
  fullDir <- Dir.[DirInfo]   
##  
##  3.  Restrict Dir only to 'include' if provided and to all but
##  'exclude' otherwise.
##
  {
    if(is.null(include)){
      excl <- subDirs %in% exclude
      pacDir <- subDirs[!excl]
      pacFull <- fullDir[!excl] 
    }
    else {
      incl <- subDirs %in% include
      pacDir <- subDirs[incl]
      pacFull <- fullDir[incl] 
    }
  }
##  
##  4.  If recursive:  
##
#  4.1.  recursive = TRUE
  {
    if(recursive){
      nSub <- length(pacDir) 
      pacList <- vector('list', nSub)
      names(pacList) <- pacDir
      if(nSub>0){
        for(i in seq(1, length=nSub))
          pacList[[i]] <- dir(pacFull[i], pattern=pattern)
      }
      return(pacList)
    }
    else
      return(pacDir)
  }
}
