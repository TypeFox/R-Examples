# $Id: glh.test.R 625 2005-06-09 14:20:30Z nj7w $

glh.test <- function( reg, cm, d=rep(0, nrow(cm)) )
{

  if( !is.matrix(cm) && !is.data.frame(cm) )
    cm <- matrix(cm, nrow=1)

  if ( !( "lm" %in% class(reg) ) )
    stop("Only defined for lm,glm objects")

  bhat <- summary.lm(reg)$coefficients[,1,drop=FALSE]
  XpX <- summary.lm(reg)$cov.unscaled
  df <- reg$df.residual
  msr <- summary.lm(reg)$sigma  # == SSE / (n-p)
  r <- nrow(cm)


  if ( ncol(cm) != length(bhat) ) stop(
                   paste( "\n Dimension of ",
                         deparse( substitute( cm ) ), ": ",
                         paste( dim(cm), collapse="x" ),
                         ", not compatible with no of parameters in ",
                         deparse( substitute( reg ) ), ": ",
                         length(bhat), sep="" ) )


  #                        -1
  #     (Cb - d)' ( C (X'X)   C' ) (Cb - d) / r
  # F = ---------------------------------------
  #                 SSE / (n-p)
  #

  Fstat <- t(cm %*% bhat - d) %*% solve((cm %*% XpX %*% t(cm))) %*% (cm %*% bhat - d) / r / msr^2

  p <- 1-pf(Fstat,r,df)

  retval <- list()
  retval$call <- match.call()
  retval$statistic <- c(F=Fstat)
  retval$parameter <- c(df1=r,df2=df)
  retval$p.value <- p
  retval$conf.int <- NULL
  retval$estimate <- cm%*%bhat
  retval$null.value <- d
  retval$method <- "Test of General Linear Hypothesis"
  retval$data.name <- deparse(substitute(reg))
  retval$matrix <- cm
  colnames(retval$matrix) <- names(reg$coef)

  class(retval) <- c("glh.test","htest")

  retval
}

print.glh.test <- function(x, digits = 4, ... )
{
    cat("\n")
    cat("\t",x$method, prefix = "\t")
    cat("\n")
    cat("Call:\n")
    print(x$call)

    if (!is.null(x$statistic))
        cat(names(x$statistic), " = ", format(round(x$statistic,
            4)), ", ", sep = "")
    if (!is.null(x$parameter))
        cat(paste(names(x$parameter), " = ", format(round(x$parameter,
            3)), ",", sep = ""), "")
    cat("p-value =",
        format.pval(x$p.value, digits = digits),
        "\n")
    cat("\n")
  }



summary.glh.test <- function(object, digits = 4, ... )
{
    cat("\n")
    cat("\t",object$method, prefiobject = "\t")
    cat("\n")
    cat("Regression: ", object$data.name, "\n")
    cat("\n")
    cat("Null Hypothesis: C %*% Beta-hat = d \n")
    cat("\n")
    cat("C matrix: \n")
    print(object$matrix, digits=digits)
    cat("\n")
    cat("d vector: \n")
    print(object$null.value, digits=digits)
    cat("\n")
    cat("C %*% Beta-hat: \n")
    print(c(object$estimate))
    cat("\n")

    if (!is.null(object$statistic))
        cat(names(object$statistic), " = ", format(round(object$statistic,
            4)), ", ", sep = "")
    if (!is.null(object$parameter))
        cat(paste(names(object$parameter), " = ", format(round(object$parameter,
            3)), ",", sep = ""), "")
    cat("p-value =",
        format.pval(object$p.value, digits = digits),
        "\n")
    cat("\n")
  }



