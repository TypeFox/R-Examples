require(DoE.wrapper)
## try out all available designs with and without factor names
set.seed(1234)
bbd.design(3)
bbd.design(4)
bbd.design(5)
bbd.design(6)
bbd.design(7)
bbd.design(4,block.name="blocks")
bbd.design(5,block.name="blocks")
bbd.design(3,factor.names=Letters[23:25])
bbd.design(3,factor.names=list(X=c(0,10),Y=c(-10,+10),Z=c(223,277)))
bbd.design(3,block.name="blocks",factor.names=Letters[23:25])
plan <- bbd.design(3,block.name="blocks",factor.names=list(X=c(0,10),Y=c(-10,+10),Z=c(223,277)))
design.info(plan)
run.order(plan)
desnum(plan)
bbd.design(4,factor.names=Letters[22:25])
bbd.design(4,factor.names=list(W="",X=c(0,10),Y=c(-10,+10),Z=c(223,277)))
bbd.design(4,block.name="blocks",factor.names=Letters[22:25])
bbd.design(4,block.name="blocks",factor.names=list(W="",X=c(0,10),Y=c(-10,+10),Z=c(223,277)))

## randomize=FALSE
run.order(bbd.design(7,randomize=FALSE))
run.order(bbd.design(4,block.name="blocks",randomize=FALSE))

## randomize with seed
plan1 <- run.order(bbd.design(7,seed=28672))
plan2 <- run.order(bbd.design(7,seed=28672))
identical(plan1, plan2)

## default levels
bbd.design(3, default.levels=c(0,100))
bbd.design(4,factor.names=list(W="",X=c(0,10),Y=c(-10,+10),Z=c(223,277)),
       default.levels=c(0,100))
