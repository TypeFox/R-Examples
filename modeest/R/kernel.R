# Author: P. Poncet

.kernel <-
function(x,
         ...)
{
  UseMethod(".kernel")
}


.kernel.default <-
function(x,
         ...)
{
  x <- tolower(x)
  x <- match.arg(x, .kernelList)  # '.kernelList' is defined in 'K.R'
  K <- eval(parse(text = paste(".kernel.", x, "(...)", sep = "")))
  if (inherits(K, "try-error")) {
    stop(paste("argument 'x' is incorrect. Function .kernel.", x, " does not exist", sep = ""))
  }

  ## Output
  return(invisible(structure(list(x = K$x,
                                  k = K$k,
                                  name = x,
                                  properties = K$properties,
                                  call = match.call()),
                             class = ".kernel")))
}


.kernel.character <- .kernel.default


#plot..kernel <-
#function(x,
#         xlab,
#         ylab,
#         ...)
#{
#  if (!inherits(x, ".kernel")) {
#    stop("argument 'x' must be of class '.kernel'")  
#  }
#  
#  if (missing(xlab)) xlab <- "x"
#  if (missing(ylab)) ylab <- paste(x$name, "'s kernel", sep = "")
#  plot.default(x$x, x$k, xlab = xlab, ylab = ylab, ...)
#}


#print..kernel <-
#function(x,
#         digits = NULL,
#         ...)
#{
#  if (!inherits(x, ".kernel")) stop("argument 'x' must be of class '.kernel'")
#  k <- x$k
#  prop <- x$properties
#  cat("class: .kernel\n")
#  cat("name of the kernel used: ", x$name, "\n")
#  cat("Int K(x) dx = ", format(prop$integral.K, digits = digits), "\n")
#  cat("Int K^2(x) dx = ", format(prop$integral.K2, digits = digits), "\n")
#  dd <- prop$differentiable
#  cc <- prop$continuous
#  if (dd == cc-1) {
#    if (cc == 1) {
#      cat("this is a continuous kernel\n")
#    } else cat("this is a ", dd, " times continuously differentiable kernel\n")
#  } else {
#    if (cc <= dd) cat("this is a ", cc, " times differentiable kernel\n")
#  }
#  cat("------------\n")
#  cat("summary statistics of the instance:\n")
#  print(c(Length = length(k), summary(k)))
#  cat("\n")
#}

