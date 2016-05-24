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
"gcv.Krig" <- function(out, lambda.grid = NA, cost = 1, 
    nstep.cv = 200, rmse = NA, verbose = FALSE, tol = 1e-05, 
    offset = 0, y = NULL, give.warnings = TRUE) {
    nt <- out$nt
    np <- out$np
    N <- out$N
    D <- out$matrices$D
    # Yet another monster function called by Krig
    #    but there just many simple steps ...
    #
    # if a y data vector is not supplied then
    # use the one in the Krig object
    if (is.null(y)) {
        u <- out$matrices$u
        shat.pure.error <- out$shat.pure.error
        pure.ss <- out$pure.ss
    }
    else {
        #with new data need to update some statistics.
        out2 <- Krig.make.u(out, y = y)
        u <- out2$u
        shat.pure.error <- out2$shat.pure.error
        pure.ss <- out2$pure.ss
    }
    if (verbose) {
        cat("u used:", fill = TRUE)
        print(u)
    }
    #
    # generate a reasonable grid of lambda based on equally spaced
    # effective degrees of freedom
    if (is.na(lambda.grid[1])) {
        temp.df <- seq(nt, (np - offset) * 0.95, , nstep.cv)
        temp.df[1] <- temp.df[1] + 0.001
        for (k in 1:nstep.cv) {
            lambda.grid[k] <- Krig.df.to.lambda(temp.df[k], D)
        }
    }
    # make sure that the grid is in sorted order
    lambda.grid <- sort(lambda.grid)
    nl <- length(lambda.grid)
    nd <- length(D)
    V <- V.model <- V.one <- lplike <- trA <- shat <- rep(NA, 
        nl)
    Dl <- rep(NA, nd)
    #
    # this is small little list used to pass information to the
    # objective functions
    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, shat.pure.error = shat.pure.error, 
        offset = offset)
    #
    # loop over lambda values for the grid search
    for (k in 1:nl) {
        #
        # all the wonderful things calculated for each lambda
        # note the use of the info list.
        V[k]       <-   Krig.fgcv(lambda.grid[k], info)
        V.one[k]   <-   Krig.fgcv.one(lambda.grid[k], info)
        V.model[k] <-   Krig.fgcv.model(lambda.grid[k], info)
        lplike[k] <-    Krig.flplike(lambda.grid[k], info)
        shat[k] <- sqrt(Krig.fs2hat(lambda.grid[k], info))
        trA[k] <-       Krig.ftrace(lambda.grid[k], D)
    }
    #
    # reformat  as a matrix with all these values.
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat, 
        lplike)
    gcv.grid <- as.data.frame(gcv.grid)
    names(gcv.grid) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model", 
        "shat", "-lnLike Prof")
    # find minima over grid ifelse used to avoid 0 length vector from which.min
    IMIN<- rep( NA, 6)
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
    IMIN[6]<- which.min( gcv.grid[["-lnLike Prof"]])  
    # NOTE IMIN indexes from smallest lambda to largest lambda in grid.        
    warningTable<- data.frame(
                    IMIN, IMIN == nl, IMIN==1,
                    gcv.grid$lambda[IMIN],
                    gcv.grid$trA[IMIN],
     row.names = c("GCV","GCV.model", "GCV.one", "RMSE", "pure error", "REML")
                     ) 
    warning<- (warningTable[,2]|warningTable[,3])&
                      (!is.na(warningTable[,1]))
    indRefine<- (!warningTable[,2]) & (!warningTable[,3]) & 
                        (!is.na(warningTable[,1]))   
    warningTable<- cbind( warning, indRefine, warningTable ) 
    names( warningTable)<- c("Warning","Refine","indexMIN", "leftEndpoint", "rightEndpoint",
                             "lambda","effdf")                                    
    # now optimze the search producing refined optima
    if (verbose) 
        print(gcv.grid)
    # setup output matrix for refined values
    lambda.est <- matrix(NA, ncol = 6, nrow = 6, dimnames = list(
         c("GCV", "GCV.model", "GCV.one", "RMSE", "pure error", "REML"), 
         c("lambda", "trA", "GCV", "shat","-lnLike Prof" , "converge")))
    # fill in grid search estimates
      for( k in 1:6){
      	if( !is.na(IMIN[k])){
      		lambda.est[k,1]<-  gcv.grid$lambda[IMIN[k]]
      	}
      }   
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
    #  6- Maxmize the restricted maxmimum likelihood (REML)
    #   standard GCV w/o replicates  
    if( verbose){
    	print( warningTable)
    }  
    if(indRefine[1]){    
        starts <- lambda.grid[IMIN[1] + c(-1,0,1)]
        out <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=Krig.fgcv, f.extra = info, tol = tol)               
        lambda.est[1,1]<-  out$x
        lambda.est[1,6]<-  out$iter
        }
    if( indRefine[2]) {
        starts <- lambda.grid[IMIN[2] + c(-1,0,1)]
        out <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=Krig.fgcv.model, f.extra = info, tol = tol)               
        lambda.est[2,1]<-  out$x 
        lambda.est[2,6]<-  out$iter     
    }
    if( indRefine[3]) {
        starts <- lambda.grid[IMIN[3] + c(-1,0,1)]
        out <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=Krig.fgcv.one, f.extra = info, tol = tol) 
         lambda.est[3, 1] <-out$x 
          lambda.est[3,6]<-  out$iter                 
        }
    if ( indRefine[6] ){
        starts <- lambda.grid[IMIN[6] + c(-1,0,1)]
        out <- golden.section.search(ax=starts[1],bx=starts[2],cx=starts[3],
                           f=Krig.flplike, f.extra = info, tol = tol)
        lambda.est[6,1]<-  out$x
        lambda.est[6,6]<-  out$iter
     }             
    if (  indRefine[4] ) {
            guess<- gcv.grid$lambda[IMIN[4]]
            lambda.rmse <- find.upcross(Krig.fs2hat, info,
                          upcross.level = rmse^2, 
                          guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
        }  
    #
    # matching estimate of sigma from reps.
    if (  indRefine[5] ) { 	    
            guess <- gcv.grid$lambda[IMIN[5]]     
            lambda.pure.error <- find.upcross(Krig.fs2hat, info, 
                    upcross.level = shat.pure.error^2, guess = guess, 
                    tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
    }
    #
    # OK done with all six methods
    # NOTE that not all may 
    # fill in return matrix with all the right stuff
    # fill in REML results
    lam.ml <- lambda.est[6, 1]
    lambda.est[6, 2] <- Krig.ftrace(lam.ml, D)
    lambda.est[6, 3] <- Krig.fgcv(lam.ml, info)
    lambda.est[6, 4] <- sqrt(Krig.fs2hat(lam.ml, info))
    lambda.est[6, 5] <- Krig.flplike(lam.ml, info)
    # fill in GCV results
    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            lambda.est[k, 2] <- Krig.ftrace(lam, D)
            if (k == 1 | k > 3) {
                lambda.est[k, 3] <- Krig.fgcv(lam, info)
                lambda.est[k, 5] <- Krig.flplike(lam, info)
            }
            if (k == 2) {
                lambda.est[k, 3] <- Krig.fgcv.model(lam, info)
            }
            if (k == 3) {
                lambda.est[k, 3] <- Krig.fgcv.one(lam, info)
            }
            lambda.est[k, 4] <- sqrt(Krig.fs2hat(lam, info))
        }
    }
    # Note that the estimate by default is
    # REML == restricted maximum likelihood. 
    if( give.warnings & any(warningTable$Warning)){
    	cat("Methods at endpoints of grid search:", fill=TRUE)
    	print(warningTable[warningTable$Warning,])
    }
    list(gcv.grid = gcv.grid, lambda.est = lambda.est, warningTable=warningTable)
}
