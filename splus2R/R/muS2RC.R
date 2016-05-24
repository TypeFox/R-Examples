# Some functions that were previously in package muS2RC
#   "Splus to R Compatibility for package muStat"
# Authors of muS2RC are Knut M. Wittkowski and Tingting Song
# Modifications here by Tim Hesterberg.

if(is.R()){
  `MC` <-
    function(f, env = NULL){
      # MC stands for "make closure". It allows the same function to
      # be used in S-PLUS and R, avoiding scoping problems in the
      # former.
      f
    }
} else {
  `MC` <-
    function(f, env = NULL) {
      # MC stands for "make closure". It allows the same function to
      # be used in S-PLUS and R, avoiding scoping problems in the
      # former.
      if (mode(f) != "function")
        stop("Not a function: ", deparse(substitute(f)))
      if (length(env) > 0 && any(names(env <- as.list(env)) == ""))
        stop("Contains unnamed arguments: ", deparse(substitute(env)))
      fargs <- ifelse1(length(f) > 1, f[1:(length(f) - 1)], NULL)
      if (any(duplicated(names(fargs <- c(fargs,env)))))
        stop(paste("duplicated arguments:", paste(names(fargs)),
              collapse = ", "))
      cf <- as.function(c(fargs, f[length(f)]))
      return(cf)
    }
}

