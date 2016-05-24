summary.bpca <- function(object,
                         presentation=FALSE,...)
{
  if (!inherits(object, 'bpca'))
    stop("Use this function only with 'bpca' class!")

  if(!presentation){

    d <- length(object$number)

    x <- list('Eigenvalue(s)' = object$eigenvalues,
              'Considered on reduction' = object$eigenvalues[object$number[1]:object$number[d]],
              'Variance retained by each' = object$eigenvalues[object$number[1]:object$number[d]]^2 / 
              sum(object$eigenvalues^2),
              'Cumulative variance retained' = cumsum(object$eigenvalues[object$number[1]:object$number[d]]^2) /
              sum(object$eigenvalues^2),
              'Prop. of total variance retained' = object$importance[1]) 

    if(object$importance[1] != object$importance[2]) 
      x$'Prop. of partial variance retained' = object$importance[2]

    class(x) <- c('summary.bpca', 'listof')

    x

  } else {

    d <- length(object$number)

    cat(' Eigenvalue(s):\t\t\t\t',
        object$eigenvalues)

    cat('\n  - Considered on reduction:\t\t',
        object$eigenvalues[object$number[1]:object$number[d]])

    cat('\n  - Variance retained by each:\t\t',
        object$eigenvalues[object$number[1]:object$number[d]]^2 /
        sum(object$eigenvalues^2))

    cat('\n  - Cumulative variance retained:\t',
        cumsum(object$eigenvalues[object$number[1]:object$number[d]]^2) /
        sum(object$eigenvalues^2))

    cat('\n  - Prop. of total variance retained:\t',
        object$importance[1])

    if(object$importance[1] != object$importance[2])
      cat('\n  - Prop. of partial variance retained:\t',
          object$importance[2])

    cat('\n')
  }
}
