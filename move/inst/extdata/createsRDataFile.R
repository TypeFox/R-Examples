##creates the 
require(move)
leroy <- move(system.file("extdata","leroy.csv.gz", package="move"))
ricky <- move(system.file("extdata","ricky.csv.gz", package="move"))
stack <- moveStack(list(leroy, ricky))
leroydbbmm <- brownian.bridge.dyn(spTransform(leroy, center=T), dimSize=45, location.error=12, time.step=5)
#rickydbbmm <- brownian.bridge.dyn(spTransform(ricky, center=T), dimSize=45, location.error=12, time.step=9)
dbbmmstack <- brownian.bridge.dyn(spTransform(stack, center=T), dimSize=45, location.error=12, ext=.3, time.step=600)
save(leroy, ricky,  leroydbbmm, dbbmmstack, file="move.RData", compress='xz')
fishers<-stack
save(fishers,file="../../data/fishers.RData", compress='xz')
save(leroy,file="../../data/leroy.RData", compress='xz')
save(ricky,file="../../data/ricky.RData", compress='xz')
save(leroydbbmm,file="../../data/leroydbbmm.RData", compress='xz')
save(dbbmmstack,file="../../data/dbbmmstack.RData", compress='xz')

