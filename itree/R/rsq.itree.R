#ALG: this code is based on the function 'rsq.rpart' but
# does an additional check to stop if a penalized tree is passed.

## This function plots the approximate r-square for the different
## splits (assumes using anova method).

## SCCS @(#)rsq.rpart.s	1.6 08/28/97

rsq.itree <- function(x) {

  if(!inherits(x,'itree')) stop("Not a legitimate itree object")

  p.rpart <- printcp(x)
  xstd <- p.rpart[,5L]
  xerror <- p.rpart[,4L]
  rel.error <- p.rpart[,3L]
  nsplit <- p.rpart[,2L]
  method <- x$method
  if(x$penalty!='none'){stop("Cannot use this function with penalization.\n")}
  if(!method=='anova') cat("Method is not anova... rsq.tree may not be applicable!\n")

  plot(nsplit, 1-rel.error, xlab='Number of Splits', ylab='R-square',
       ylim=c(0,1), type='o')
  par(new=TRUE)
  plot(nsplit, 1-xerror, type='o', ylim=c(0,1),lty=2, xlab=' ', ylab=' ')
  legend(0,1, c('Apparent','X Relative'), lty=1:2)


  ylim <- c(min(xerror-xstd) -.1, max(xerror + xstd) + .1)
  plot(nsplit, xerror, xlab='Number of Splits', ylab='X Relative Error',
       ylim=ylim, type='o')
  segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
  invisible()

  }

