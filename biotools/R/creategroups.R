creategroups <-
function(x, ngroups, sizes, fun = mean,
tol = 1e-2, maxit = 200)
{
   stopifnot(is.numeric(x))
   if (is.null(names(x)))
      stop("'x' must be a named vector!")
   stopifnot(sum(sizes) == length(x))
   stopifnot(ngroups == length(sizes))
   iter <- 0
   repeat {
      g <- sample(rep(1:ngroups, sizes))
      s <- split(x, g)
      label <- split(names(x), g)
      dif <- diff(sapply(s, fun))
      dif <- mean(abs(dif))
      iter <- iter + 1
      if (iter > maxit)
         stop("'maxit' reached!")
      if (dif <= tol)
         break()
   }

   out <- list(covar = deparse(substitute(x)),
      func = deparse(substitute(fun)),
      val.func = sapply(s, fun), niter = iter,
      labels = label, groups = s)
   class(out) <- "creategroups"
   return(out)
}

# -------------------------------------------
# print method
print.creategroups <- 
function (x, digits = 4L, quote = TRUE, ...) 
{
   cat("\nCreating homogeneous groups \n")
   cat("------------------------------------------------------------------\n")
   cat("Covariate:", x$covar, 
      "\n\nGroups: \n")
   print(x$labels)
   cat("Objective function (equality of):", x$func, "\n")
   print(x$val.func)
   cat("\nNumber of iterations to convergence:", x$niter, "\n")
   invisible(x)
}