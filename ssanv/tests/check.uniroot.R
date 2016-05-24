library(ssanv)
## test function when root fall on an integer
rootfunc<-function(x){ 150 - x }

## test: f(lower)=0
uniroot.integer(rootfunc,c(150,200),print.steps=TRUE)
## test: f(upper)=0
uniroot.integer(rootfunc,c(15,150),print.steps=TRUE)
## test: step into zero on way up f(lower+64)=0
uniroot.integer(rootfunc,c(150-64,10^3),print.steps=TRUE)
## test: step into zero on the way back
uniroot.integer(rootfunc,c(150-64+32,10^3),print.steps=TRUE)
## test: step into zero after switching more than once
uniroot.integer(rootfunc,c(150-64+32-16,10^3),print.steps=TRUE)
## test: other checks
uniroot.integer(rootfunc,c(150-64+32-15,10^3),print.steps=TRUE)

### rerun those checks for root not exactly on integer

rootfunc<-function(x){ 150.5 - x }

uniroot.integer(rootfunc,c(150,200),print.steps=FALSE)$root
uniroot.integer(rootfunc,c(15,151),print.steps=FALSE)$root
uniroot.integer(rootfunc,c(150-64,10^3),print.steps=FALSE)$root
uniroot.integer(rootfunc,c(150-64+32,10^3),print.steps=FALSE)$root
uniroot.integer(rootfunc,c(150-64+32-16,10^3),print.steps=FALSE)$root
uniroot.integer(rootfunc,c(150-64+32-15,10^3),print.steps=FALSE)$root

### rerun those checks for pos.side


uniroot.integer(rootfunc,c(150,200),pos.side=TRUE)$root
uniroot.integer(rootfunc,c(15,151),pos.side=TRUE)$root
uniroot.integer(rootfunc,c(150-64,10^3),pos.side=TRUE)$root
uniroot.integer(rootfunc,c(150-64+32,10^3),pos.side=TRUE)$root
uniroot.integer(rootfunc,c(150-64+32-16,10^3),pos.side=TRUE)$root
uniroot.integer(rootfunc,c(150-64+32-15,10^3),pos.side=TRUE)$root




rootfunc<-function(x){ 166-x }
## test: problem
uniroot.integer(rootfunc,c(5,10^5),print.steps=TRUE)

# test: make sure maxiter is working
#uniroot.integer(rootfunc,c(5,10^5),maxiter=3,print.steps=TRUE)

