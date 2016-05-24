"MCbound.tsprt" <-
function(parms,conf.level=.99){
    bparms<-tSPRT.to.Bvalue(parms)
    #print(bparms)
    out<-MCbound.Bvalue(bparms["Nmax"],bparms["alpha"],bparms["e0"],bparms["e1"],conf.level)
    ### output tsprt parameterization instead of Bvalue parameterization
    out$type<-"tsprt"
    out$parms<-parms
    out
}

