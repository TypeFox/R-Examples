
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

LKrig.MLE <- function(x, y, ..., LKinfo, use.cholesky=NULL, par.grid = NULL, 
    lambda.profile = TRUE, verbose = FALSE, lowerBoundLogLambda=-16,
                      nTasks=1, taskID=1, tol=.005) {
    LKrig.args <- c(list(x = x, y = y), list(...))
    # at this point LKinfo has the correct value for the number of multiresolution levels
    par.grid <- LKrig.make.par.grid(par.grid = par.grid, LKinfo = LKinfo) 
    NG <- length(par.grid$alpha)
    # adjust the for loop from 1:NG if this
    # function is called through Rmpi (i.e. if nTasks!=1)
    # idea is that task with ID will work on
    # NG1 through NG2 -- a fraction of the total number in par.grid
     if( NG < nTasks ){
        stop("Too many tasks for the number of parameter values")
      }
    indexTemp<- round(seq( 0,NG,,nTasks+1))
    NG1<- indexTemp[taskID] + 1
    NG2<- indexTemp[taskID+1]
    #
    out <- matrix(NA, nrow = NG, ncol = 9,
                  dimnames = list(NULL, c("EffDf", "lnProfLike", "GCV", 
                  "sigma.MLE", "rho.MLE", "llambda.MLE", "lnLike", "counts value", "grad")))
    optim.counts <- rep(NA, 2)
    # evaluate parameters but do an optimzation over lambda
    lnProfileLike.max <- -1e+20
    lnLike.eval<- NULL
    llambda.opt<- NA
    for (k in NG1:NG2) {
    # if starting value is missing use the optimum from previous fit
    # this only really makes sense if other parameters have some sort of continuity from k-1 to k.
        llambda.start<- ifelse (is.na( par.grid$llambda[k]), llambda.opt,  par.grid$llambda[k] )  
     # first fit to get cholesky symbolic decomposition
        LKinfo.temp<- do.call("LKinfoUpdate",c( list(LKinfo=LKinfo),list(
                                   a.wght = (par.grid$a.wght[[k]]),
                                    alpha =  (par.grid$alpha[[k]]),
                                       nu = par.grid$nu[k],
                                   lambda = exp(llambda.start)
                                    ) )
                              )      
      if( verbose){
      	cat("LKrig.MLE: initial LKinfo object:", fill=TRUE)
      	print( LKinfo.temp)
      }                             
    # for first pass save symbolic decomposition for M
        if( k == NG1 ){ use.cholesky <- NULL}
           obj <- do.call("LKrigFindLambda",
                      c( LKrig.args,
                        list(       LKinfo = LKinfo.temp,
                            lambda.profile = lambda.profile,
                              use.cholesky = use.cholesky,
                                       tol = tol,
                                   verbose = verbose
                            )
                        )
                        )
          if( verbose){
          	cat("LKrig.MLE: find lambda grid value", k,  fill=TRUE)
          	print( obj)
          }              
     # Note: if lambda.profile == FALSE the starting lambda is passed through as the "optimal" one
        llambda.opt<- obj$summary["llambda.opt"]
     # compare to current largest likelihood and update the LKinfo.MLE list if bigger.
        if ( obj$summary["lnProfLike"] > lnProfileLike.max ) {
              lnProfileLike.max <- obj$summary["lnProfLike"] 
              LKinfo.MLE <- obj$LKinfo
              lambda.MLE<- exp( llambda.opt)
        }
     # save summary results from this set of parameters.
       out[k, ] <- obj$summary
       lnLike.eval <- c( lnLike.eval, list(obj$lnLike.eval))
     #   
       if (verbose) { print( c(k,out[k,])) }
    }
    return(list(
                 summary = out,
                par.grid = par.grid,
                  LKinfo = LKinfo, 
              LKinfo.MLE = LKinfo.MLE,
             lnLike.eval = lnLike.eval,
              lambda.MLE = lambda.MLE,
                    call = match.call(),
                  taskID = taskID
                )
           )
}

