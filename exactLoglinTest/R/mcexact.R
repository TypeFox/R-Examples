mcexact <- function(formula,
                    data,
                    stat = gof,
                    dens = hyper,
                    nosim = 10 ^ 3,
                    method = "bab",
                    savechain = FALSE,
                    tdf = 3,
                    maxiter = nosim,
                    p = NULL,
                    batchsize = NULL){  
  args <- build.mcx.obj(formula,
                        data,
                        stat,
                        dens,
                        nosim,
                        method,
                        savechain,
                        tdf,
                        maxiter,
                        p,
                        batchsize)
  update(args, savechain = savechain)
}

##This is no longer necessary
#.First.lib <- function(lib, pkg){
#  library.dynam("exactLoglinTest", pkg, lib)
#}
