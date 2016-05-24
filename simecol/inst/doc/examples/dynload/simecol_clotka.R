################################################################################
##
## Example on using compiled code efficiently *with* simecol
##
##
## thomas.petzoldt@tu-dresden.de
##
################################################################################

library("simecol")

# compile C++ code within R
# (requires installed compiler)
# on Windows: http://www.murdoch-sutherland.com/Rtools/
system("R CMD SHLIB clotka.c")


modeldll <- dyn.load("clotka.dll")

clotka <- new("odeModel",

  ## note that this main does not contain the equations directly
  ## but returns information where these can be found
  main = function(time, init, parms) {
     # list with dllname, func, nout, [jacfunc]
     list(lib     = "clotka",
          func    = "dlotka",
          jacfunc = NULL,
          nout    = 2)
  },

  ## parms, times, init are provided as usual, enabling
  ## scenario control like for "ordinary" simecol models
  parms  = c(k1=0.2, k2=0.2, k3=0.2),
  times  = c(from=0, to=100, by=0.5),
  init   = c(prey=0.5, predator=1),

  ## a special solver interface that evaluates funclist
  ## and passes its contents directly to the lsoda
  ## in the "compiled function" mode
  solver = function(init, times, funclist, parms, ...) {
    f <- funclist()
    as.data.frame(lsoda(init, times, func=f$func,
      parms = parms, dllname = f$lib, jacfunc=f$jacfunc, nout = f$nout, ...)
    )
  }
)

clotka <- sim(clotka)

## the two graphivs on top are the states
## the other are additional variables returned by the C code
## (for demonstration purposes here)
plot(clotka)

## Another simulation with more time steps
times(clotka)["to"] <- 1000
plot(sim(clotka))

## another simulation with intentionally reduced accuracy
## for testing
plot(sim(clotka, atol=1))

dyn.unload(as.character(modeldll[2]))
