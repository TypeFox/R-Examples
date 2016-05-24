require("DoE.base")
### test programs for crossing designs
### and parameter designs
oa12 <- oa.design(nlevels=c(2,6,2), factor.names=c("first","second","third"),randomize=FALSE)
oa4 <- oa.design(nlevels=c(2,2,2), factor.names=Letters[4:6],randomize=FALSE)
oa4rep <- oa.design(nlevels=c(2,2,2), factor.names=Letters[7:9], repl=3,randomize=FALSE)
oa4reprepeat.only <-oa.design(nlevels=c(2,2,2), repl=3, factor.names=Letters[10:12], repeat.only=TRUE,randomize=FALSE)
cross1 <- cross.design(oa12,oa4,oa4rep,randomize=FALSE)
cross1
summary(cross1)
cross2 <- cross.design(oa12,oa4rep,oa4,randomize=FALSE)
cross2
summary(cross2)
altern <- c(35,55,80)
#alter <- factor(alter,levels=alter)
cross3 <- cross.design(oa12,oa4,altern,randomize=FALSE)
cross3
summary(cross3)
alterc <- c("jung","mittel","alt")
#alter <- factor(alter,levels=alter)
cross4 <- cross.design(oa12,oa4,alterc,randomize=FALSE)
cross4
summary(cross4)

cross5 <- cross.design(oa12,oa4,oa4reprepeat.only,randomize=FALSE)
cross5
summary(cross5)
cross6 <- cross.design(oa4rep,oa4reprepeat.only,randomize=FALSE)
cross6
summary(cross6)
design.info(cross2)
factor.names(cross2) <- Letters[1:design.info(cross2)$nfactors]
summary.data.frame(cross2)

param1 <- param.design(oa12, oa4, direction="wide")
param1
summary(param1)