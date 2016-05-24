#This is a longer version of the example to demostrate more #functions relating to MCMC parameter tracing, output #custimization and interaction with R envionment during model #fitting.

rm(list = ls())

require(NPBayesImpute)
data('NYexample')

#create the model in one step
#CreateModel(X,MCZ,K,Nmiss_Max, alphaa,aplhab)
model <- CreateModel(X,MCZ,50,200000,0.25,0.25)


#retrieve the disjointed MCZ
model$MCZ


#list all traceable parameters
model$traceable

#set parameters to be traced (optional)
#trace the quoted 4 parameters and set maximum number of traces to 1000
model$SetTrace(c("index","nu","alpha","psi"),1000)

#run 100 burnins, 200 mcmc iterations and thinning every 1 iterations
model$Run(100,200,1)

#run another 100 iterations, but thinning every 10 iterstions
model$Run(0,100,10)

#one can stop the run using R stop button and then resume using
model$Resume()

#retireve traces
trace <- model$GetTrace() #(optional)

#list the traced parameters
model$traced

#show the current iteration
model$CurrentIteration

#enable/disable Tracer
model$EnableTracer <- TRUE
#odel$EnableTracer <- FALSE

#check tracer enable/disable satus
model$EnableTracer

#retrieve parameters from the final iteration
result <- model$snapshot

#convert ImputedX matrix to dataframe, using proper factors/names etc.
IX <- GetDataFrame(result$ImputedX,X)
View(IX)

#names(result)

#retrieve the selected parameters
r <- model$Parameters(c("nu","alpha"))

#use coda for a simple traceplot
require(coda)
mcmc_obj = mcmc(trace$alpha)
traceplot(mcmc_obj)


