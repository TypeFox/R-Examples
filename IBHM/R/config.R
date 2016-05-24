
# IBHM configuration ------------------------------------------------------

ConfigureIBHM  <- function(stop.criterion = IterationSC(3),                           
                           weighting.function = function(y, w.par){ 0.01+dnorm(y,sd=abs(w.par))},
                           scal.optim = 'multi.CMAES',
                           scal.optim.params = list(retries=3, inner=list(maxit=50, stopfitness=-1)),
                           scal.candidates = c('dot.pr','radial','root.radial'),                                                      
                           activ.optim = 'multi.CMAES',
                           activ.optim.params = list(retries=3, inner=list(maxit=100, stopfitness=-1)),
                           activ.candidates = c('tanh','logsig','lin'),
                           jit=TRUE,                           
                           verbose=FALSE,
                           final.estimation = 'all',
                           final.estimation.x = NULL,
                           final.estimation.y = NULL,
                           final.estimation.maxit = 100
){    
  list( sc = stop.criterion,
        wf = weighting.function,
        final.estimation = switch(final.estimation,
                                  weights = OptimizeAllWeights,
                                  all = OptimizeAllParams,
                                  stop('Unknow final.estimation type: ',final.estimation)),
        final.estimation.x = final.estimation.x,
        final.estimation.y = final.estimation.y,
        final.estimation.maxit = final.estimation.maxit,
        jit=jit,
        verbose=verbose,
        scal = list(                  
          optim = switch(scal.optim, 
                         CMAES = OptimizeScalCMAES,
                         multi.CMAES = OptimizeScalMultiCMAES,
                         DE = OptimizeScalDE,
                         multi.DE = OptimizeScalMultiDE,
                         multi.NM = OptimizeScalMultiNM,
                         stop('Unknown scal optimization method: ',activ.optim)),
          optim.params = scal.optim.params,
          candidates = ScalFunctions(scal.candidates)
        ),
        activ = list(          
          optim = switch(activ.optim, 
                         CMAES = OptimizeActivCMAES,
                         multi.CMAES = OptimizeActivMultiCMAES,
                         DE = OptimizeActivDE,
                         multi.DE = OptimizeActivMultiDE,
                         multi.NM = OptimizeActivMultiNM,
                         stop('Unknown activ optimization method:',activ.optim)),
          optim.params = activ.optim.params,
          candidates = ActivationCandidates(activ.candidates)
        )
  )
}