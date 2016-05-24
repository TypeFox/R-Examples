require(DoE.wrapper)

ccd.augment(FrF2(8,4,randomize=FALSE),randomize=FALSE)
ccd.augment(FrF2(8,3,randomize=FALSE),randomize=FALSE)

## resolution V design as basis
plan <- FrF2(16,5,ncenter=4,randomize=FALSE)
set.seed(2424)
y <- rexp(20)
r.plan <- add.response(plan,y)
ccd.augment(plan,4,randomize=FALSE)
## augmenting design with response
ccd.augment(r.plan,4,randomize=FALSE)

## estimable design
## basic order
basic <- FrF2(8,4, estimable=c("CD"),res3=TRUE,randomize=FALSE)
ccd.augment(basic,4,randomize=FALSE)
## reshuffled
## only one swap
reshuffled <- FrF2(8,4, estimable=c("AD"),res3=TRUE,randomize=FALSE)
ccd.augment(reshuffled,6,randomize=FALSE)
reshuffled <- FrF2(8,4, estimable=c("AD"),res3=TRUE, default.levels=c(0,100),
    factor.names=Letters[22:25],randomize=FALSE)
reshuffled <- FrF2(8,4, estimable=c("AD"),res3=TRUE, default.levels=c(0,100),
    factor.names=list(T=c(30,50),U=c(24,26),V=c(100,400),W=c(30,75)),randomize=FALSE)
ccd.augment(reshuffled,6,randomize=FALSE)
desnum(ccd.augment(reshuffled,6,randomize=FALSE))
## more reshuffling
reshuffled.big <- FrF2(32,7, estimable=c("AC","BC","AB"),randomize=FALSE,
   factor.names=Letters[19:25],default.levels=c(10,30))
ccd.reshuffled.big <- ccd.augment(reshuffled.big,6,randomize=FALSE)
reshuffled.big <- FrF2(32,7, estimable=c("AC","BC","AB"),randomize=FALSE,
   factor.names=list(T=c(30,50),U=c(24,26),V=c(100,400),W=c(30,75), X=c(0.1,0.7),Y="",Z=""),
   default.levels=c(10,30),repl=2)
ccd.reshuffled.big <- ccd.augment(reshuffled.big,6,randomize=FALSE)
ccd.reshuffled.big
desnum(ccd.reshuffled.big)

## a properly designed plan
plan <- ccd.augment(FrF2(16,5,randomize=FALSE),6,randomize=FALSE)
set.seed(23232)
y <- rexp(38)
r.plan <- add.response(plan, y)
rsm(y~SO(A,B,C,D,E), r.plan)

## replicated designs
plan <- ccd.augment(FrF2(16,5,repl=2,randomize=FALSE),6,randomize=FALSE)

## blocked design
plan <- ccd.augment(FrF2(32,5,blocks=4,randomize=FALSE),6,randomize=FALSE)
plan <- ccd.augment(FrF2(16,5,blocks=2,randomize=FALSE),6, bbreps=c(1,2),randomize=FALSE)

## two different versions of big designs because of different content of base.design
planblockpickbig <- FrF2(64,gen=c(7,11,14),blocks=16,alias.block.2fis=TRUE,randomize=FALSE)
set.seed(2323)
y <- rnorm(64)
planblockpickbig <- add.response(planblockpickbig,y)
plan <- ccd.augment(planblockpickbig,n0=1,randomize=FALSE)
planblockpickbig <- FrF2(64,nfactors=9,blocks=16,alias.block.2fis=TRUE,randomize=FALSE)
plan <- ccd.augment(planblockpickbig,randomize=FALSE)