# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
gcv.sreg<- function(out, lambda.grid = NA, cost = 1, 
    nstep.cv = 80, rmse = NA, offset = 0, trmin = NA, trmax = NA, 
    verbose = FALSE, tol = 1e-05,
    give.warnings=TRUE) {
    shat.pure.error <- out$shat.pure.error
    pure.ss <- out$pure.ss
    nt <- 2
    np <- out$np
    N <- out$N
    out$cost <- cost
    out$offset <- offset
    # find good end points for lambda coarse grid.
    if (is.na(trmin)) 
        trmin <- 2.05
    if (is.na(trmax)) 
        trmax <- out$np * 0.95
    if (is.na(lambda.grid[1])) {
        l2 <- sreg.df.to.lambda(trmax, out$xM, out$weightsM)
        l1 <- sreg.df.to.lambda(trmin, out$xM, out$weightsM)
        lambda.grid <- exp(seq(log(l2), log(l1), , nstep.cv))
    }
    if (verbose) {
        cat("endpoints of coarse lamdba grid", fill = TRUE)
        cat(l1, l2, fill = TRUE)
    }
    # build up table of coarse grid serach results for lambda
    # in the matrix gcv.grid
    nl <- length(lambda.grid)
    V <- V.model <- V.one <- trA <- MSE <- RSS.model <- rep(NA, 
        nl)
    #   loop through lambda's and compute various quantities related to
    #   lambda and the fitted spline.
    for (k in 1:nl) {
        temp <- sreg.fit(lambda.grid[k], out, verbose = verbose)
        RSS.model[k] <- temp$rss
        trA[k] <- temp$trace
        V[k] <- temp$gcv
        V.one[k] <- temp$gcv.one
        V.model[k] <- temp$gcv.model
    }
    # adjustments to columns of gcv.grid
    RSS <- RSS.model + pure.ss
    shat <- sqrt(RSS/(N - trA))
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat)
    dimnames(gcv.grid) <- list(NULL, c("lambda", "trA", "GCV", 
        "GCV.one", "GCV.model", "shat"))
        gcv.grid<- as.data.frame( gcv.grid)
    if (verbose) {
        cat("Results of coarse grid search", fill = TRUE)
        print(gcv.grid)
    }
    lambda.est <- matrix(NA, ncol = 5, nrow = 5,
           dimnames = list(
           c("GCV","GCV.model", "GCV.one", "RMSE", "pure error"),
           c("lambda","trA", "GCV", "shat", "converge")))
    # now do various refinements for different flavors of finding
    # a good value for lambda the smoothing parameter
    ##### traditional leave-one-out
    IMIN<- rep( NA, 5)
    IMIN[1]<- which.min(    gcv.grid$GCV ) 
    IMIN[2]<- ifelse( is.na(shat.pure.error), NA,
                 which.min(gcv.grid$GCV.model) )
    IMIN[3]<- which.min(    gcv.grid$GCV.one)
    if( is.na( rmse)){
    	IMIN[4] <- NA
    }
    else{
       rangeShat<-  range( gcv.grid$shat) 
       IUpcross<- max( (1:nl)[gcv.grid$shat< rmse] )
      IMIN[4]<- ifelse( (rangeShat[1]<= rmse)&(rangeShat[2] >=rmse),
                                        IUpcross, NA)
    }
    IMIN[5]<- ifelse( is.na(shat.pure.error), NA,
                       which.min(abs(gcv.grid$shat-shat.pure.error)) ) 
    # NOTE IMIN indexes from smallest lambda to largest lambda in grid.        
    warningTable<- data.frame(
                    IMIN, IMIN == nl, IMIN==1,
                    gcv.grid$lambda[IMIN],
                    gcv.grid$trA[IMIN],
     row.names = c("GCV","GCV.model", "GCV.one", "RMSE", "pure error") ) 
    warning<- (warningTable[,2]|warningTable[,3])&
                      (!is.na(warningTable[,1]))
    indRefine<- (!warningTable[,2]) & (!warningTable[,3]) & 
                        (!is.na(warningTable[,1]))   
    warningTable<- cbind( warning, indRefine, warningTable ) 
    names( warningTable)<- c("Warning","Refine","indexMIN", "leftEndpoint", "rightEndpoint",
                             "lambda","effdf")
     if( verbose){
     	print(warningTable)
     }
   # fill in grid search estimates
      for( k in 1:5){
      	if( !is.na(IMIN[k])){
      		lambda.est[k,1]<-  gcv.grid$lambda[IMIN[k]]
      	}
      }                              
    # now optimze the search producing refined optima
    #
    # now step through the many different ways to find lambda
    # This is the key to these choices:
    #  1- the usual GCV proposed by Craven/Wahba
    #  2- GCV where data fitting is collapsed to the mean for
    #     each location and each location is omitted 
    #  3- True leave-one-out even with replicated observations
    #  4- Match estimate of sigma to external value supplied (RMSE)
    #  5- Match estimate of sigma from the estimate based the 
    #     pure error sum of squares obtained by the observations
    #     replicated at the same locations
    #test<- sreg.fit(.1, out)
    #print( test)
     if(indRefine[1]){    
        starts <- lambda.grid[IMIN[1] + c(-1,0,1)]
        outGs <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=sreg.fgcv, f.extra = out, tol = tol)               
        lambda.est[1,1]<-  outGs$x
        lambda.est[1,5]<-  outGs$iter
        }
    if( indRefine[2]) {
        starts <- lambda.grid[IMIN[2] + c(-1,0,1)]
        outGs <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=sreg.fgcv.model, f.extra = out, tol = tol)               
        lambda.est[2,1]<-  outGs$x 
        lambda.est[2,5]<-  outGs$iter     
    }
    if( indRefine[3]) {
        starts <- lambda.grid[IMIN[3] + c(-1,0,1)]
        outGs <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=sreg.fgcv.one, f.extra = out, tol = tol) 
         lambda.est[3, 1] <-outGs$x 
          lambda.est[3,5]<-  outGs$iter                 
        }                               
     if (  indRefine[4] ) {
            guess<- gcv.grid$lambda[IMIN[4]]
            lambda.rmse <- find.upcross(sreg.fs2hat, out,
                          upcross.level = rmse^2, 
                          guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
        } 
         if (  indRefine[5] ) { 	    
            guess <- gcv.grid$lambda[IMIN[5]]     
            lambda.pure.error <- find.upcross(sreg.fs2hat, out, 
                    upcross.level = shat.pure.error^2, guess = guess, 
                    tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
    }
   if (verbose) {
        cat("All forms of estimated lambdas so far", fill = TRUE)
        print(lambda.est)
    }
    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            temp <- sreg.fit(lam, out)
            lambda.est[k, 2] <- temp$trace
            if ((k == 1) | (k > 3)) {
                lambda.est[k, 3] <- temp$gcv
            }
            if (k == 2) {
                lambda.est[k, 3] <- temp$gcv.model
            }
            if (k == 3) {
                lambda.est[k, 3] <- temp$gcv.one
            }
            lambda.est[k, 4] <- temp$shat
        }
    }
    if( give.warnings & any(warningTable$Warning)){
    	cat("Methods at endpoints of grid search:", fill=TRUE)
    	print(warningTable[warningTable$Warning,])
    }
    list(gcv.grid = gcv.grid, lambda.est = lambda.est,
     warningTable=warningTable)
}
