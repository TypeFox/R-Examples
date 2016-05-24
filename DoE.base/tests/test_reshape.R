## test reshape functions in DoE.base
require(DoE.base)

a <- fac.design(nlevels=c(2,4,3),repl=2,repeat.only=TRUE,randomize=FALSE)
b <- oa.design(nlevels=c(2,6,2),repl=3,repeat.only=TRUE,randomize=FALSE)
c <- a
response.names(c) <- c("Y1", "Y2")
d <- b
d$Y <- rexp(12)
d$Z <- runif(12)
response.names(d) <- c("Y","Z")

aw <- reptowide(a)
bw <- reptowide(b)
cw <- reptowide(c)
dw <- reptowide(d)

al <- reptolong(aw)
bl <- reptolong(bw)
cl <- reptolong(cw)
dl <- reptolong(dw)

aw
al
bw
bl

design.info(aw)$responselist
design.info(bw)$responselist
design.info(cw)$responselist
design.info(dw)$responselist

design.info(aw)$response.names
design.info(bw)$response.names
design.info(cw)$response.names
design.info(dw)$response.names

design.info(al)$response.names
design.info(bl)$response.names
design.info(cl)$response.names
design.info(dl)$response.names
