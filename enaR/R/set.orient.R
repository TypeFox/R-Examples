#' set.orient --- globally reorients matrices
#' INPUT = matrix orientation (rc or cr)
#' OUTPUT = sets the expected orientation of matrices
#' 
#' M. Lau | Feb 2013
#' ---------------------------------------------------

set.orient <- local({
  orientation <- 'rc'
  warn <- ''
  f <- function(x=c('rc','school')){
    if (any(x%in%c('rc','school','internal'))){
      orientation <<- x[1]
      if (x[1] == 'school'){
        warning('NOTE: output of functions from a particular analytical school will be returned in the standard orientation of that school.')
      }else if (x[1]=='internal'){
        orientation <<- 'school'
      }
    }else{
      warning('Unknown orientation.')
    }
  }
})
