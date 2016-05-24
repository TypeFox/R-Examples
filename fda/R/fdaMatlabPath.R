fdaMatlabPath <- function(R.matlab) {
##
## 1.  path2fdaM = path to ~R/library/fda/Matlab/fdaM
##
  path2fdaM <- system.file('Matlab/fdaM', package='fda')
##
## 2.  dirs2add = dirs(path2fdaM, ...)
##
  d2a. <- dirs(path2fdaM, exclude='^@|^private$',
             full.names=TRUE, recursive=TRUE)
  dirs2add <- c(path2fdaM, d2a.) 
##
## 3.  requires(R.matlab)?
##
  missRmat <- missing(R.matlab)
  if(missRmat)R.matlab <- TRUE 
  if(R.matlab){
    if(require(R.matlab))
      dirs2add <- c(system.file('externals', package='R.matlab'),
                    dirs2add)
    else
      if(!missRmat)
        warning('Package R.matlab is not installed and can not be',
                ' included in the Matlab path.')
  }
##
## 4.  Create Matlab 'addpath' commands.
##
  d2a <- paste("addpath('", dirs2add, "');", sep='')
##
## 5.  write file
##
  writeLines(d2a, 'fdaMatlabPath.m')
##
## 6.  Done
##
  invisible(d2a)  
}
