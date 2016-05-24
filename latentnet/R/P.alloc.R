## Lets the user call P.free.all from R. Only used for debugging.
P.free.all<-function() invisible(try(.C("P_free_all", PACKAGE="latentnet")))
