# sequential bootstrapping

source('nuclearBootstrapInit.R')

R <- 20000
stime <- system.time(boot(nuke.data, nuke.fun, R=R, m=1, fit.pred=new.fit, x.pred=new.data))
cat('Elapsed time:', stime[3], '\n')
