require(FrF2)

summary(FrF2(16,7),brief=TRUE)
summary(FrF2(32,7,blocks=4),brief=TRUE)
summary(FrF2(32,9,blocks=4,alias.block.2fis=TRUE),brief=TRUE)
summary(FrF2(32,9,blocks=c("AB","AC"),alias.block.2fis=TRUE),brief=TRUE)
summary(FrF2(64,7,blocks=4),brief=TRUE)
summary(plan1 <- FrF2(32,7,estimable=c("AB","CD")),brief=TRUE)
summary(plan2 <- FrF2(nruns= 64 ,nfactors= 10 , estimable= c( "AB","AC","AD","AE","AF",
       "AG","AH","AJ","AK","BC","BD","BE","BF","BG","BH","BJ","BK" ) ), brief=TRUE)

C1 <- compromise(10, 1:2, msg=FALSE)
summary(plan3 <- FrF2(nruns= 64 ,nfactors= 10 , 
    estimable= C1$requirement, perms=C1$perms.full, clear=FALSE ), brief=TRUE)

### blockpick.big, a design that uses map
plan <- FrF2(64,15,blocks=16,alias.block.2fi=TRUE)
summary(plan,brief=TRUE)

### various split-plot setups
planfull <-FrF2(32,5,WPs=4,nfac.WP=2,factor.names=Letters[21:25])
summary(planfull, brief=TRUE)
design.info(planfull)[-which(names(design.info(planfull))=="FrF2.version")]
try(planfull <-FrF2(32,5,WPs=4,nfac.WP=3,factor.names=Letters[21:25]))

plan0 <-FrF2(32,7,gen=c(6,15),WPs=4,WPfacs=c("B","C","F"),factor.names=Letters[19:25])
summary(plan0,brief=TRUE)
design.info(plan0)[-which(names(design.info(plan0))=="FrF2.version")]

plan <-FrF2(32,7,WPs=4,nfac.WP=3,factor.names=Letters[19:25])
summary(plan, brief=TRUE)
generators(plan)
design.info(plan)[-which(names(design.info(plan))=="FrF2.version")]

plan2 <- FrF2(32,7,gen=c(6,15),WPs=4,nfac.WP=3,WPfacs=c(2,3,6),factor.names=Letters[19:25])
summary(plan2, brief=TRUE)
generators(plan2)
design.info(plan2)[-which(names(design.info(plan2))=="FrF2.version")]
