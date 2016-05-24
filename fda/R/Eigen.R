Eigen <- function(x, symmetric, only.values = FALSE, EISPACK = FALSE,
      valuenames ){
##
## 1.  symmetric?  
##
  N <- nrow(x)
  if (missing(symmetric)) 
    symmetric <- isSymmetric.matrix(x)  
##
## 2.  eigen
##
  ev <- eigen(x, symmetric, only.values = FALSE, EISPACK = FALSE)
##
## 3.  rNames
##
  rNames <- rownames(x)
  if(is.null(rNames))
    rNames <- {
      if(symmetric) paste('x', 1:N, sep='')
      else paste('xrow', 1:N, sep='')
    }
##
## 4.  parse valuenames
##   
  {
    if(missing(valuenames)){
      cNames <- colnames(x)
      if(is.null(cNames))
        cNames <- {
          if(symmetric) paste('x', 1:N, sep='')
          else paste('xcol', 1:N, sep='')
        }
      if(symmetric){
        valuenames <- {
          if(all(rNames==cNames))paste('ev', 1:N, sep='')
          else cNames
        }
      }
      else
        valuenames <- cNames
    }
    else{
      if(length(valuenames)<N)
        valuenames <- paste(valuenames, 1:N, sep='')
      else {
        if(length(valuenames)>N)
          warning('length(valuenames) = ', length(valuenames),
                  ' > nrow(x) = ', N,
                  '; using only the first ', N)
        valuenames <- valuenames[1:N]
      }
    }
  }  
##
## 5.  rNames
##
  names(ev$values) <- valuenames
  if(!only.values)
    dimnames(ev$vectors) <- list(rNames, valuenames)
##
## 6.  Done
##
  ev
}
