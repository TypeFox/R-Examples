findSubsample <- 
function(x, size, fun = sd, minimize = TRUE, niter = 10000)
{
   dataname <- deparse(substitute(x))
   stopifnot(is.numeric(x))
   x <- as.vector(x)
   n <- length(x)
   if (size > n)
     stop("'size' must be less than the length of 'x'!")
   stopifnot(is.function(fun))
   if (niter > choose(n, size))
      stop("'niter' misdefined!")
   if (is.null(names(x))) names(x) <- 1:n

   l <- list()
   stat <- numeric()
   for(i in 1:niter) {
      l[[i]] <- sample(x, size)
      stat[i] <- fun(l[[i]])      
   }

   if(minimize) {
      index <- which.min(stat)
      minmax <- 'minimized'
   } else {
      index <- which.max(stat)
      minmax <- 'maximized'
   }

   out <- list(dataname = dataname,
      niter = niter, fun = deparse(substitute(fun)), 
      stat = stat[index], criterion = minmax, 
      subsample = l[[index]], labels = names(l[[index]]))
   class(out) <- "findSubsample"
   return(out)
}

print.findSubsample <- function(x, digits = 4L, quote = TRUE, ...)
{
   cat("\n")
   cat(strwrap("Finding an Optimized Subsample", prefix = "\t"), sep = "\n")
   cat("\nData:", x$dataname, 
      "\nIterations:", x$niter,
      "\nObjective function:", x$fun, 
      "=", x$stat, "(", x$criterion, ")",
      "\nSubsample: \n")
   print(x$subsample)
   invisible(x)
}
