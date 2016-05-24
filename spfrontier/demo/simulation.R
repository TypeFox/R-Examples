params000 <- list(n=c(50,100,200,300),
                  beta0=5,
                  beta1=2,
                  beta2=3,
                  sigmaV=0.1, 
                  sigmaU=0.5)

params000T <- c(params000, list(mu=1))
params100 <- c(params000, list(rhoY=0.2))
params100T <- c(params000T, list(rhoY=0.2))

params110 <- c(params100, list(rhoV=0.4))
params101 <- c(params100, list(rhoU=0.4))
params111 <- c(params110, list(rhoU=0.4))

params010 <- params110
params010$rhoY <- NULL

params011 <- params111
params011$rhoY <- NULL

params001 <- params011
params001$rhoV <- NULL

ctrl <- list(true.initial=F, seed=999, cores=detectCores())
res000 <- ezsimspfrontier(100, params = params000,  inefficiency = "half-normal",logging = "info", control=ctrl)
save(res000, file="res000.rData")

res000T <- ezsimspfrontier(100, params = params000T, inefficiency = "truncated",logging = "info", control=ctrl)
save(res000T, file="res000T.rData")

res100 <- ezsimspfrontier(100, params = params100,  inefficiency = "half-normal",logging = "info", control=ctrl)
save(res100, file="res100.rData")

res100A <- ezsimspfrontier(100, params = params100,  inefficiency = "half-normal",logging = "info", 
                               control=c(ctrl,list(ignoreWy=T)))
save(res100A, file="res100A.rData")

res100T <- ezsimspfrontier(100, params = params100T, inefficiency = "truncated",logging = "info", control=ctrl)
save(res100T, file="res100T.rData")

params001A <- c(params001, list(rhoY=0))
res001A <- ezsimspfrontier(100, params = params001A, inefficiency = "half-normal",logging = "info", 
                           control=c(ctrl,list(ignoreWu=T,replaceWyWu=T)))
save(res001A, file="res001A.rData")

params010A <- c(params010, list(rhoY=0))
res010A <- ezsimspfrontier(1, params = params010A, inefficiency = "half-normal",logging = "info", 
                           control=c(ctrl,list(ignoreWv=T,replaceWyWv=T)))
save(res010A, file="res010A.rData")


ctrl <- list(true.initial=TRUE, seed=0, cores=detectCores()-1)

#Takes ~20 hrs on 8 cores
params001$n <- c(50,100)
res001_50_100 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res001_50_100, file="res001_50_100.rData")
params001$n <- c(200)
res001_200 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res001_200, file="res001_200.rData")
params001$n <- c(300)
res001_300 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res001_300, file="res001_300.rData")

params010$n <- c(50,100)
res010_50_100 <- ezsimspfrontier(100, params = params010, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res010_50_100, file="res010_50_100.rData")
params010$n <- c(200)
res010_200 <- ezsimspfrontier(100, params = params010, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res010_200, file="res010_200.rData")
params010$n <- c(300)
res010_300 <- ezsimspfrontier(100, params = params010, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res010_300, file="res010_300.rData")

#All tests above work as appropriate
res010 <- ezsimspfrontier(100, params = params010, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res010, file="res010.rData")

params010$n <- c(300)
res010 <- ezsimspfrontier(100, params = params010, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res010, file="SimE05_res010.rData")

res110 <- ezsimspfrontier(100, params = params110, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res110, file="SimE07_res110.rData")

res101 <- ezsimspfrontier(100, params = params101, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res101, file="SimE08_res101.rData")

res011 <- ezsimspfrontier(100, params = params011, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res011, file="SimE09_res011.rData")

res111 <- ezsimspfrontier(100, params = params111, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res111, file="SimE10_res111.rData")
