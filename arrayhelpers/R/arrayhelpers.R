##' Little helper functions to work with arrays
##' @name arrayhelpers-package
##' @docType package
##' 
{

  if (require ("svUnit", quietly = TRUE)){
    `.test<-` <- svUnit::`test<-`
  } else {
    `.test<-` <- function (f, value) {
      class (value) <-  c ("svTest", "function")
      attr (f, "test") <- value
      f
    }
  }

    
### test data
    ## vector
    v <- 1 : 3
    names (v) <- letters [1 : 3]

    ## matrix
    m <- matrix (1:6, 2,
                 dimnames = list (rows = letters [1:2],
                   columns = LETTERS [1:3])) 

    ## array
    a <- array (1 : 24, 4 : 2,
                dimnames = list (rows = letters [1:4],
                  columns = LETTERS [1:3],
                  d3 = 1:2)
                )
}



