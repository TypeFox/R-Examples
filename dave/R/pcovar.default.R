pcovar.default <-
function(veg,y,...) {
    o.pcovar<- pcoatest(veg,y)
    class(o.pcovar)<- "pcovar"
    o.pcovar$call<- match.call()
    cat("Call:\n")
    print(o.pcovar$call)
    o.pcovar
  }
