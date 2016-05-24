### R code from vignette source 'cll.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: genimage
###################################################
library(potts)
nrow <- 32
ncol <- 32
ncolor <- 4
theta.true <- rep(0, ncolor+1)
theta.true[ncolor+1] <- log(1 + sqrt(ncolor))
x <- matrix(sample(ncolor, nrow*ncol, replace=TRUE), nrow=nrow,
            ncol=ncol)
out <- potts(packPotts(x, ncolor), theta.true, nbatch=1000, blen=1)
x <- unpackPotts(out$final)


###################################################
### code chunk number 2: potts-fig
###################################################
image(x)


###################################################
### code chunk number 3: calct
###################################################
t_stat <- calc_t(x, ncolor)
t_stat


###################################################
### code chunk number 4: calctcache
###################################################
t_cache_mple <- generate_t_cache(x, ncolor, t_stat, nrow*ncol, 1, 
                                  singleton)
t_cache_mple[[1]]
t_cache_two <- generate_t_cache(x, ncolor, t_stat, nrow*ncol/2, 2, 
                                twopixel.nonoverlap)
t_cache_two[[1]]


###################################################
### code chunk number 5: calccll
###################################################
composite.ll(theta.true[-1], t_stat, t_cache_mple)
gr.composite.ll(theta.true[-1], t_stat, t_cache_mple)

composite.ll(theta.true[-1], t_stat, t_cache_two)
gr.composite.ll(theta.true[-1], t_stat, t_cache_two)


###################################################
### code chunk number 6: optimcll
###################################################
theta.initial <- 1:ncolor
optim.mple <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                    t_stat, t_cache_mple, method="BFGS", 
                    control=list(fnscale=-1))
optim.mple$par

optim.two <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                   t_stat, t_cache_two, method="BFGS", 
                   control=list(fnscale=-1))
optim.two$par

theta.true


