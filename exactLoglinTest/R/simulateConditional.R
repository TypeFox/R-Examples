simulateConditional <- function(formula,
                                 data,
                                 dens = hyper,
                                 nosim = 10 ^ 3,
                                 method = "bab",
                                 tdf = 3,
                                 maxiter = nosim,
                                 p = NULL,
                                 y.start = NULL){
  args <- build.mcx.obj(formula = formula,
                        data = data,
                        dens = dens,
                        nosim = nosim,
                        method = method,
                        tdf = tdf,
                        maxiter = maxiter,
                        p = p,
                        batchsize = 10#have to have this for error checking
                        )
  
  if (method == "bab")
    simtable.bab(args)
  else if (method == "cab"){
    if (is.null(y.start))
      return(simtable.cab(args))
    else
      return(simtable.cab(args, y.start = y.start))
  }
}
