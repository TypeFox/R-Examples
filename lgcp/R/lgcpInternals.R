##' at function
##' @param t change in time parameter, see Brix and Diggle (2001)
##' @param mu mean
##' @param theta parameter beta in Brix and Diggle
##' @return ...
at <- function(t,mu,theta){
    return(mu*(1-exp(-t*theta)))
}



##' bt.scalar function
##' @param t change in time, see Brix and Diggle (2001)
##' @param theta parameter beta in Brix and Diggle
##' @return ...
bt.scalar <- function(t,theta){
    return(exp(-t*theta))
}



##' d.func function
##' @param mat1il matrix 1
##' @param mat2jk matrix 2
##' @param i index matrix 1 number 1
##' @param j index matrix 2 number 1
##' @param l index matrix 1 number 2
##' @param k index matrix 2 number 2
##' @return ...
d.func <- function(mat1il,mat2jk,i,j,l,k){
# warning global variables
  return(sqrt(mat1il[i,l]^2+mat2jk[j,k]^2))
}






##' gu function
##' @param u distance
##' @param sigma variance parameter, see Brix and Diggle (2001)
##' @param phi scale parameter, see Brix and Diggle (2001)
##' @param model correlation type, see ?CovarianceFct
##' @param additionalparameters vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @return this is just a wrapper for CovarianceFct
gu <- function(u,sigma,phi,model,additionalparameters){
        return(CovarianceFct(x=u,param=c(mean=0,variance=sigma^2,nugget=0,scale=phi,additionalparameters),model=model))

}








##' lgcpPredictSpatialINLA function
##'
##' -------------------------------------------------------
##' !IMPORTANT! after library(lgcp) this will be a dummy function. 
##' In order to use, type getlgcpPredictSpatialINLA() at the console. This will download and install the true function.
##' -------------------------------------------------------
##'
##' The function \code{lgcpPredictSpatialINLA} performs spatial prediction for log-Gaussian Cox Processes using the integrated nested Laplace approximation.
##'
##' The following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s)}{\mathcal Y(s)} be a spatial Gaussian process and \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space. 
##' Cases occur at spatial positions \eqn{x \in W}{x \in W} 
##'  according to an inhomogeneous spatial Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x)}{R(x)},
##'   The number of cases, \eqn{X_{S}}{X_{S}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S} \sim \mbox{Poisson}\left\{\int_S R(s)ds\right\}}{X_{S} \sim \mbox{Poisson}\left\{\int_S R(s)ds\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005) (but ignoring temporal variation), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\exp\{\mathcal Y(s,t)\}.}{R(s,t) = \lambda(s)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1.}{\int_W\lambda(s)d s=1.}
##'
##'
##' Before calling this function, the user must decide on the parameters, spatial covariance model, spatial discretisation,
##' fixed spatial (\eqn{\lambda(s)}{\lambda(s)}) component and whether or not any output is
##' required. Note there is no autorotate option for this function.
##'
##' @param sd a spatial point pattern object, see ?ppp
##' @param ns size of neighbourhood to use for GMRF approximation ns=1 corresponds to 3^2-1=8 eight neighbours around each point, ns=2 corresponds to 5^2-1=24 neighbours etc ...
##' @param model.parameters values for parameters, see ?lgcppars
##' @param spatial.covmodel correlation type see ?CovarianceFct 
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @param cellwidth width of grid cells on which to do MALA (grid cells are square). Note EITHER gridsize OR cellwidthe must be specified.
##' @param gridsize size of output grid required. Note EITHER gridsize OR cellwidthe must be specified.
##' @param spatial.intensity the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param optimverbose logical whether to print optimisation details of covariance matching step
##' @param inlaverbose loogical whether to print the inla fitting procedure to the console
##' @param generic0hyper optional hyperparameter list specification for "generic0" INLA model. default is list(theta=list(initial=0,fixed=TRUE)), which effectively treats the precision matrix as known.
##' @param strategy inla strategy
##' @param method optimisation method to be used in function matchcovariance, default is "Nelder-Mead". See ?matchcovariance 
##' @return the results of fitting the model in an object of class \code{lgcpPredict}
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle (2013). Journal of Statistical Software, 52(4), 1-40. URL http://www.jstatsoft.org/v52/i04/
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{lgcpPredict} \link{KinhomAverage}, \link{ginhomAverage}, \link{lambdaEst}, \link{muEst}, \link{spatialparsEst}, \link{thetaEst},  
##' \link{spatialAtRisk}, \link{temporalAtRisk}, \link{lgcppars}, \link{CovarianceFct}, \link{mcmcpars}, \link{setoutput} 
##' \link{print.lgcpPredict}, \link{xvals.lgcpPredict}, \link{yvals.lgcpPredict}, \link{plot.lgcpPredict}, \link{meanfield.lgcpPredict},
##' \link{rr.lgcpPredict}, \link{serr.lgcpPredict}, \link{intens.lgcpPredict},   
##' \link{varfield.lgcpPredict}, \link{gridfun.lgcpPredict}, \link{gridav.lgcpPredict}, \link{hvals.lgcpPredict}, \link{window.lgcpPredict},
##' \link{mcmctrace.lgcpPredict}, \link{plotExceed.lgcpPredict}, \link{quantile.lgcpPredict}, \link{identify.lgcpPredict}, \link{expectation.lgcpPredict},
##' \link{extract.lgcpPredict}, \link{showGrid.lgcpPredict},
##' @export 
    
lgcpPredictSpatialINLA <- function( sd,
                                    ns,
            					    model.parameters=lgcppars(),
            					    spatial.covmodel="exponential",
            					    covpars=c(),
            					    cellwidth=NULL,
            					    gridsize=NULL,
            					    spatial.intensity,				
            					    ext=2,
            					    optimverbose=FALSE,
            					    inlaverbose=TRUE,
            					    generic0hyper=list(theta=list(initial=0,fixed=TRUE)),
            					    strategy="simplified.laplace",
            					    method="Nelder-Mead"){

    cat("  lgcpPredictSpatialINLA not installed yet ...\n")
    cat("  type getlgcpPredictSpatialINLA() to get it.\n")


    # note this is a dummy function. In order to use, type getlgcpPredictSpatialINLA() at the console.

}

            					    
