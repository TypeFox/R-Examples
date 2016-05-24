##########################################################################
## sample from a K-dimensional two-parameter item response model with
## probit link. This is just a wrapper function that calls
## MCMCordfactanal.
##
## Andrew D. Martin
## Washington University
##
## Kevin M. Quinn
## Harvard University
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## June 8, 2003
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

"MCMCirtKd" <-
  function(datamatrix, dimensions, item.constraints=list(),
           burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA,
           alphabeta.start = NA, b0=0, B0=0,
           store.item=FALSE, store.ability=TRUE,
           drop.constant.items=TRUE, ... ) {

    datamatrix <- as.matrix(datamatrix)   
    
    post <- MCMCordfactanal(x=datamatrix, factors=dimensions,
                            lambda.constraints=item.constraints,
                            burnin=burnin, mcmc=mcmc, thin=thin,
                            tune=NA, verbose=verbose, seed=seed,
                            lambda.start=alphabeta.start,
                            l0=b0, L0=B0, store.lambda=store.item,
                            store.scores=store.ability,
                            drop.constantvars=drop.constant.items,
                            model="MCMCirtKd")
    return(post)
  }

