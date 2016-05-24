pmatch2 <- function(x, table){
##
## 1.  Check x  
##
  if(is.list(x)){
    stop('x is a list; not allowed.')
  }
##
## 2.  Check table 
##
  if(is.list(table)){
    stop('table is a list; not allowed.')
  }    
##
## 3.  create out list 
##
  nx <- length(x)
  out <- vector('list', nx) 
  names(out) <- x
##
## 4.  process element by element 
##
  for(ix in seq(length=nx)){
#    xi <- which(x[ix] %in% table)
    xi <- which(table %in% x[ix])
    if(length(xi)<1){
      oops <- c(grep('(', x[ix], fixed=TRUE),  
                grep('[', x[ix], fixed=TRUE), 
                grep('{', x[ix], fixed=TRUE) )
      if(length(oops)>0){
        warning('strange character string = ', 
                x[ix]) 
        xi <- grep(x[ix], table, fixed=TRUE)
      } else {
        xi1 <- paste0('^', x[ix])
        xi <- grep(xi1, table)
        if(length(xi)<1)xi <- grep(x[ix], table)
      }
    }
    out[[ix]] <- xi 
  }
##
## 5.  done 
##
  out
}