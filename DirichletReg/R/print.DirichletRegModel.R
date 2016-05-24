print.DirichletRegModel <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  .wd <- getOption("width")

  names(x$coefficients) <- x$coefnames

  if(x$optimization$convergence == 3) cat("\n",strwrap("CAUTION! Possible convergence problems!",.wd),"\n",sep="")
  if(x$optimization$convergence > 3) stop("\n",paste(strwrap(paste("\nOptimization did not converge in",x$optimization$bfgs.it,"+",x$optimization$iterations,"iterations and exited with code",x$optimization$convergence),.wd),sep="\n",collapse="\n"))

  if(interactive()) writeLines("")

  writeLines("Call:")
  writeLines(strwrap(deparse(x$call, width.cutoff=500), .wd))
  writeLines(paste0("using the ", x$parametrization, " parametrization"))

  writeLines("")

  writeLines(paste0("Log-likelihood: ",format(x$logLik,digits=digits)," on ",x$npar," df (", x$optimization$bfgs.it," BFGS + ",x$optimization$iterations," NR Iterations)"))

  writeLines("")

  coef.ind <- cumsum(x$n.vars)

  if(x$parametrization == "common"){

    for(i in 1:length(x$varnames)){
      writeLines(paste0(rep("-", min(41L, .wd)),collapse=""))
      writeLines(paste0("Coefficients for variable no. ",i,": ",x$varnames[i]))
      print.default(format(x$coefficients[ifelse(i==1,1,coef.ind[i-1]+1):coef.ind[i]], digits = digits), print.gap = 2L, quote = FALSE)
    }
    writeLines(paste0(rep("-", min(41L, .wd)), collapse=""))

  } else {

    printed.var <- 1
    set.size    <- ncol(x$X[[1]])

    writeLines("MEAN MODELS:")

    for(i in 1:length(x$varnames)){
      if(i == x$base){
        writeLines(paste0(rep("-", min(41L, .wd)),collapse=""))
        writeLines(paste0("Coefficients for variable no. ",i , ": ", x$varnames[i]))
        writeLines("- variable omitted (reference category) -")
      } else {
        writeLines(paste0(rep("-", min(41L, .wd)),collapse=""))
        writeLines(paste0("Coefficients for variable no. ",i , ": ", x$varnames[i]))
        print.default(format(x$coefficients[printed.var:(printed.var+set.size-1)], digits=digits), print.gap=2, quote=F)

        printed.var <- printed.var + set.size
      }
    }

    writeLines(paste0(rep("-", min(41L, .wd)), collapse=""))

    writeLines("")

    writeLines("PRECISION MODEL:")

    writeLines(paste0(rep("-", min(41L, .wd)), collapse=""))
    print.default(format(x$coefficients[printed.var:length(x$coefficients)], digits=digits), print.gap=2, quote=F)
    writeLines(paste0(rep("-", min(41L, .wd)), collapse=""))
  }

  if(interactive()) writeLines("")

}
