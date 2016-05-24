##' Modification of the function  \code{\link[KrigInv]{integration_design}} from the package \code{\link[KrigInv]{KrigInv}} to 
##' be usable for SUR-based optimization. Handles two or three objectives.
##' Available important sampling schemes: none so far.
##' @title Function to build integration points (for the SUR criterion)
##' @param SURcontrol Optional list specifying the procedure to build the integration points and weights. 
##'        Many options are possible; see 'Details'.
##' @param d The dimension of the input set. If not provided \code{d} is set equal to the length of \code{lower}.
##' @param lower Vector containing the lower bounds of the design space.
##' @param upper Vector containing the upper bounds of the design space.
##' @param model A list of kriging models of \code{km} class.
##' @param min.prob This argument applies only when importance sampling distributions are chosen. 
##'       For numerical reasons we give a minimum probability for a point to
##'       belong to the importance sample. This avoids probabilities equal to zero and importance sampling
##'       weights equal to infinity. In an importance sample of \code{M} points, the maximum weight becomes 
##'       \code{1/min.prob * 1/M}.
##' @references
##' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}.
##' @seealso \code{\link[GPareto]{GParetoptim}} \code{\link[GPareto]{crit_SUR}} \code{\link[KrigInv]{integration_design}}
##' @details 
##' Options for the \code{SURcontrol} list are :
##' \itemize{
##'   \item A) If nothing is specified, \code{100 * d} points are chosen using the Sobol sequence;
##'   \item B) One can directly set the field \code{integration.points} (\code{p * d} matrix) for prespecified integration points. 
##'        In this case these integration points and the corresponding vector \code{integration.weights} will be used 
##'        for all the iterations of the algorithm;
##'   \item C) If the field \code{integration.points} is not set then the integration points are renewed at each iteration. 
##'        In that case one can control the number of integration points \code{n.points} (default: \code{100*d}) and a specific 
##'        distribution \code{distrib}. Possible values for distrib are: "\code{sobol}", "\code{MC}" and "\code{SUR}"
##'         (default: "\code{sobol}"):
##'         \itemize{
##'           \item C.1) The choice "\code{sobol}" corresponds to integration points chosen with the Sobol sequence in dimension \code{d} (uniform weight);
##'           \item C.2) The choice "\code{MC}" corresponds to points chosen randomly, uniformly on the domain;
##'           \item C.3) The choice "\code{SUR}" corresponds to importance sampling distributions (unequal weights). \cr
##'         When important sampling procedures are chosen, \code{n.points} points are chosen using importance sampling among a discrete 
##'         set of \code{n.candidates} points (default: \code{n.points*10}) which are distributed according to a distribution \code{init.distrib} 
##'         (default: "\code{sobol}"). Possible values for \code{init.distrib} are the space filling distributions "\code{sobol}" and "\code{MC}" 
##'         or an user defined distribution "\code{spec}". The "\code{sobol}" and "\code{MC}" choices correspond to quasi random and random points 
##'         in the domain. If the "\code{spec}" value is chosen the user must fill in manually the field \code{init.distrib.spec} to specify 
##'         himself a \code{n.candidates * d} matrix of points in dimension \code{d}.
##'         }
##' }
##' 
##'        
##'        
##' @return 
##' A list with components:
##' \itemize{
##' \item{\code{integration.points}}{ \code{p x d} matrix of p points used for the numerical calculation of integrals}
##' \item{\code{integration.weights}}{ a vector of size \code{p} corresponding to the weight of each point. If all the points are equally 
##' weighted, \code{integration.weights} is set to \code{NULL}}
##' }
##' @importFrom randtoolbox sobol
##' @export

integration_design_optim <- function(SURcontrol=NULL,d=NULL,lower,upper,model=NULL,min.prob=0.001){
  #####################################################################################  
  # Generic function to build integration points for some criterion
  # Modification of the function integration_design from the package KrigInv to 
  # be usable for SUR-based optimization. Handles one constraint and 2/3 objectives
  # Available important sampling schemes: none so far
  #####################################################################################  
  
  result <- NULL
  if(is.null(d)) d <- length(lower)
  if (length(lower) != length(upper) ){
    print("Error in integration_Parameters: 'lower' and 'upper' must have the same length")
    return(NULL)
  }
  
  #Trivial case 1
  if(is.null(SURcontrol)){
    #nothing has been specified, thus we use default values
    n.int.points<-d*100
    integration.points <- lower+sobol(n=n.int.points,dim=d)*(upper-lower)
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(!is.null(model)) colnames(integration.points)<- colnames(model[[1]]@X)
    
    result$integration.points <- integration.points
    result$integration.weights <- NULL
    return(result)
  }
  
  #Trivial case 2
  if(!is.null(SURcontrol$integration.points)){
    #integration points pre-specified
    #     if(!is.null(model) && d>1) colnames(SURcontrol$integration.points) <- colnames(model@X)
    result$integration.points  <- SURcontrol$integration.points
    result$integration.weights <- SURcontrol$integration.weights
    return(result)
  }
  
  #non trivial cases:
  if(is.null(SURcontrol$n.points)) SURcontrol$n.points <- d*100
  if(is.null(SURcontrol$distrib)) SURcontrol$distrib <- "MC"
  
  if(SURcontrol$distrib=="sobol"){
    integration.points <- lower+sobol(n=SURcontrol$n.points,dim=d)*(upper-lower)
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
    result$integration.points <- integration.points
    result$integration.weights<-NULL
    return(result)
  }
  
  if(SURcontrol$distrib=="MC"){
    integration.points <- lower+matrix(runif(d*SURcontrol$n.points),ncol=d)*(upper-lower)
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
    result$integration.points <- integration.points
    result$integration.weights<-NULL
    return(result)
  }
  
  
  if(SURcontrol$distrib=="SUR"){
    if(is.null(SURcontrol$n.candidates)) SURcontrol$n.candidates <- SURcontrol$n.points*10
    if(is.null(SURcontrol$init.distrib)) SURcontrol$init.distrib <- "MC"
    
    #generation of the initial candidates points:
    if(SURcontrol$init.distrib=="sobol") initial.integration.points <- t(lower+t(sobol(n=SURcontrol$n.candidates,dim=d))*(upper-lower))
    if(SURcontrol$init.distrib=="MC") initial.integration.points <- t(lower+t(matrix(runif(d*SURcontrol$n.candidates),ncol=d))*(upper-lower))
    if(SURcontrol$init.distrib=="spec") initial.integration.points <- SURcontrol$init.distrib.spec
    
    if(d==1) initial.integration.points<-matrix(initial.integration.points,ncol=1)
    
    #prediction on these initial candidate points
    if(is.null(model)){
      print("Error in integration_Parameters: for the 'SUR' importance sampling distribution, 
            you must set the argument 'model'")
      return(NULL)
    }
    
    #--------------------------------------------------------------
    if(SURcontrol$distrib=="SUR"){
      Tau.n <- prob.of.non.domination(model=model, integration.points=initial.integration.points) 
    }
    #--------------------------------------------------------------
    Tau.n.sum <- sum(Tau.n)
    if(Tau.n.sum==0) Tau.n.sum <- 1
    prob.n <- pmax(Tau.n/Tau.n.sum,min.prob/SURcontrol$n.candidates)
    prob.n <- prob.n/sum(prob.n)
    weight.n <- 1/(prob.n*SURcontrol$n.candidates*SURcontrol$n.points)
    
    prob.n.copy <- c(0,prob.n)
    prob.n.cum  <- cumsum(prob.n.copy)
    
    my.indices <- findInterval(runif(SURcontrol$n.points),prob.n.cum,all.inside=TRUE)
    integration.points <- initial.integration.points[my.indices,]
    integration.weights <- weight.n[my.indices]
    
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(SURcontrol$n.points==1) integration.points <- matrix(integration.points,ncol=d)
    
    if(!is.null(model)){
      if (length(model) > 1){
        colnames(integration.points) <- colnames(model[[1]]@X)
      } else {
        colnames(integration.points)<- colnames(model@X) 
      }
    }
    result$integration.points  <- integration.points
    result$integration.weights <- integration.weights
    return(result)
  }
  }
