
library(CIDnetworks)

data(dolphins)
model.er <- CID.Gibbs(dolphins,thin=1)
#model.er <- CID.Gibbs(dolphins,verbose=1)
#model.er <- CID.Gibbs(dolphins,verbose=0)
#model.lsm <- CID.Gibbs(dolphins,burnin=200,thin=10,
#                       components=list(LSM(2)),
#                       verbose=0)

#data(Lazega)
#net <- Lazega$Friendship
#model.er <- CID.Gibbs(net)
#model.er <- CID.Gibbs(net,verbose=0)
#model.lsm <- CID.Gibbs(net,burnin=200,thin=10,
#                        components=list(LSM(2)),
#                        verbose=0)

#model.er <- CID.Gibbs(net,auto.converge=TRUE)
#model.lsm <- CID.Gibbs(net,draws=200,thin=10,
#                       components=list(LSM(2)),
#                       verbose=2,auto.converge=TRUE)
#plot(model.lsm)
#test <- model.lsm$CID.mean$gibbs.value(model.lsm$results)


#model.sbm <- CID.Gibbs(net,draws=300,thin=10,
#                       components=list(SBM(3)),
#                       verbose=2,auto.converge=TRUE)
#plot(model.sbm)
#model.lsm$DIC
#model.lsm$CID.mean$log.likelihood
#model.lsm$CID.mean$components[[1]]$value()
#model.lsm$CID.mean$components[[1]]$mult.factor
#model.lsm$CID.mean$intercept
#model.lsm$CID.object$components[[1]]$mult.factor
