#' Constructor for a Mixed Membership Model Object
#' 
#' Constructor for a \code{mixedMemModel} object which can be used for analysis 
#' in the \code{mixedMem} package.
#' 
#' The function returns an object of \code{mixedMemModel} class. This object contains dimensions of the model,
#' the observed data, and the model parameters. Once a \code{mixedMemModel} object is created,
#' the specified model can be fit for the data using the \code{mmVarFit} function. If the inputs are inconsistent (ie if dimensions do not match,
#' or if observations and distribution types are not compatible,\code{mixedMemModel} will throw an error. For additional details on usage, and model
#' assumptions, see the corresponding vignette "Fitting Mixed Membership Models Using \code{mixedMem}".
#' 
#' @param Total the number of individuals in the sample.
#' @param J the number of variables observed on each individual.
#' @param Rj a vector of length J specifying the number of repeated measurements
#'  for each variable.
#' @param Nijr an array with dimension (\code{Total}, \code{J}, \code{max(Rj)}) indicating the number
#'  of ranking levels for each variable and each replication. For multinomial and Bernoulli
#'  variables, \code{Nijr}[i,j,r] = 1. For rank variables, \code{Nijr}[i,j,r] indicates the
#'  number of alternatives ranked.
#' @param K the number of sub-populations.
#' @param Vj a vector of length \code{J} specifying the number of possible candidates
#'  for each variable. For a Bernoulli variable \code{Vj}[j] = 1. 
#' @param alpha a vector of length \code{K} which is the parameter for Dirichlet
#'  membership distribution.
#' @param theta an array with dimensions (\code{J},\code{K},\code{max(Vj)}) which governs the variable
#'  distributions. The parameter \code{theta}[j,k,] governs the distribution of variable J for a complete member of sub-population k. For instance, if variable j is a Bernoulli variable, theta[j,k,1] is the probability of success; if
#'  variable j is a multinomial variable, \code{theta}[j,k, 1:Vj[j]] is the probability for each of the \code{Vj}[j] categories; if variable j
#'  is a rank variable, \code{theta}[j,k, 1:Vj[j]] are the support parameters for each of the \code{Vj}[j] alternatives. Since the dimension of the relevant parameters
#'  can differ across variables, any unused elements of the array should be set to 0, while all other elements should be positive.
#' @param phi an array with dimensions (\code{Total},\code{K}) which specifies the variational
#'  parameters for the membership vectors, \eqn{\lambda}. The default group membership initialization is uniform across all groups (phi[i,k] = 1/K for all k).
#'  The default initialization is highly recommended.
#' @param delta an array with dimensions (\code{Total},\code{J},\code{max(Rj)}, \code{max(Nijr)}, \code{K}) which specifies
#'  the variational parameters for the context vectors Z. The default initialization is
#'   uniform across all sub-populations (\code{delta}[i, j, r, n, k] = 1/K for all k).
#' @param dist a vector of length \code{J} specifying variable types. Options are
#'  "bernoulli", "multinomial" or "rank" corresponing to the distributions
#'   of the observed variables.
#' @param obs an array with dimensions (\code{Total},\code{J},\code{max(Rj)}, \code{max(Nijr)}) corresponding to 
#' the observed data. For Bernoulli random variables, the data consist of 0/1's. 
#' For multinomial or rank data the data consist of integers 0,1,\ldots,(Vj[j]-1).
#' @param fixedObs an array with dimensions (1, \code{J}, \code{max(Rj)}, \code{max(Nijr)}) corresonding to the observed
#' responses for a fixed group in the extended GoM model from Erosheva et al (2007)
#' @param P scalar between 0 and 1 corresponding to initial value for the proportion of individuals in the fixed group
#' @param beta scalar between 0 and 1 corresponding to the initial value of beta, the conditional probability of being in the fixed group
#' for an individual who's observed responses match the fixed group.   
#' @return returns an object of class \code{mixedMemModel}
#' @examples
#' set.seed(123)
#' Total <- 50 # 50 Individuals
#' J <- 3 # 3 different variables
#' # distributions of each variable
#' dist <- c("multinomial", "bernoulli", "rank") 
#' # 100 repeated measures of the multinomial, 5 repeated measures of the
#' # Bernoulli, 1 repeated measure of the rank
#' Rj <- c(100, 5, 1) 
#' 
#' K <- 4 # 4 sub-populations
#' alpha <- rep(.5, K) #hyperparameter for dirichlet distribution
#' 
#' # Number of categories/alternatives for each variable. For the Bernoulli, Vj = 1
#' Vj <- c(10, 1, 4) 
#' 
#' 
#' theta <- array(0, dim = c(J, K, max(Vj)))
#' # Parameters governing multinomial
#' theta[1,,] <- gtools::rdirichlet(K, rep(.3, Vj[1]))
#' #parameters governing Bernoulli
#' theta[2,,] <- cbind(rbeta(K, 1,1), matrix(0, nrow = K, ncol = Vj[1]-1))
#' theta[3,,] <- cbind(gtools::rdirichlet(K, rep(.3, Vj[3])),
#'  matrix(0, nrow = K, ncol = Vj[1]-Vj[3]))
#' 
#' # Items selected for each observation. For Multinomial and Bernoulli, this is always 1
#' # For rank data, this will be the number of alternatives ranked
#' Nijr = array(0, dim = c(Total, J, max(Rj)))
#' Nijr[,1,c(1:Rj[1])] = 1 # N_ijr is 1 for multinomial variables
#' Nijr[,2,c(1:Rj[2])] = 1 # N_ijr is 1 for Bernoulli variables
#' Nijr[,3,c(1:Rj[3])] = sample.int(Vj[3], size = Total, replace = TRUE)
#' 
#' # generate random sample of observations
#' sampleMixedMem <- rmixedMem(Total, J, Rj, Nijr, K, Vj,
#' dist, theta, alpha)
#' 
#' ## Initialize a mixedMemModel object
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj,
#'  Nijr= Nijr, K = K, Vj = Vj,dist = dist, obs = sampleMixedMem$obs,
#'   alpha = alpha, theta = theta)
#' # Look at Summary of the initialized model
#' summary(test_model)
#' # Plot the current values for theta
#' plot(test_model) 
#' @export

mixedMemModel = function(Total, J, Rj, Nijr, K, Vj, alpha, theta, phi = NULL, delta = NULL, dist, obs,
                         fixedObs = NULL, P = NULL, beta = NULL)
{
  # Checks if model defaults are used and fills in defaults
  if(is.null(delta))
  {
    delta = array(0, dim = c(Total,J, max(Rj), max(Nijr), K))
    for(i in 1:Total)
    {
      for(j in 1:J)
      {
        for(r in 1:Rj[j])
        {
          for(n in 1:Nijr[i,j,r])
          {
            delta[i,j,r,n,] = rep(1/K,K)
          }
        }
      }
    }
  }
  
  if(is.null(phi))
  {
    phi = array(1/K, dim = c(Total,K))
  }
  #put objects in a list
  if(is.null(fixedObs)){
    model_obj = list(Total, J, Rj, Nijr, K, Vj,
                     alpha + 1 -1, theta + 1 -1, phi + 1 -1, delta + 1 -1,
                     dist, obs);
    names(model_obj) = c("Total", "J", "Rj", "Nijr", "K", "Vj", "alpha","theta", "phi",
                         "delta", "dist" ,"obs")
  } else {
    model_obj = list(Total, J, Rj, Nijr, K, Vj,
                     alpha + 1 -1, theta + 1 -1, phi + 1 -1, delta + 1 -1,
                     dist, obs, fixedObs,
                     P+ 1 -1, beta + 1 -1);
    names(model_obj) = c("Total", "J", "Rj", "Nijr", "K", "Vj", "alpha","theta", "phi",
                         "delta", "dist" ,"obs", "fixedObs", "P", "beta")    
  }
  class(model_obj) = "mixedMemModel"
  
  dimnames(model_obj$theta) <- list(paste("Var", c(1:J)),
                                    paste("Group", c(1:K)),
                                    paste("Cand", c(0:(max(Vj)-1))))
  
  names(model_obj$alpha) <- paste("Group", c(0:(K-1)))
  
  names(model_obj$Vj) <- paste("Var", c(1:J))
  names(model_obj$Rj) <- paste("Var", c(1:J))
  names(model_obj$dist) <- paste("Var", c(1:J))
  
  dimnames(model_obj$phi) <- list(paste("Ind", c(1:Total)),
                                  paste("Group", c(0:(K-1))))
  
  dimnames(model_obj$delta) <- list(paste("Ind", c(1:Total)),
                                    paste("Var", c(1:J)),
                                    paste("Rep", c(1:max(Rj))),
                                    paste("Rank", c(1:max(Nijr))),
                                    paste("Group", c(0:(K-1))))
  
  dimnames(model_obj$Nijr) <- list(paste("Ind", c(1:Total)),
                                   paste("Var", c(1:J)),
                                   paste("Rep", c(1:max(Rj))))
  
  dimnames(model_obj$obs) <- list(paste("Ind", c(1:Total)),
                                  paste("Var", c(1:J)),
                                  paste("Rep", c(1:max(Rj))),
                                  paste("Rank", c(1:max(Nijr))))
  #check for valid model parameters
  checkModel(model_obj)
  return(model_obj)
}


#' Summary of a Mixed Membership Model
#' 
#' Generic S3 summary function for \code{mixedMemModel} class. 
#' 
#' \code{summary} provides a summary of the given \code{mixedMemModel} object. The function 
#' prints the ELBO, the dimensions of the model and each variable type.
#'  
#' @param object the mixedMemModel object to be summarized
#' @param ... additional parameters
#' @seealso mixedMemModel
#' @export
summary.mixedMemModel = function(object,...)
{
  cat("== Summary for Mixed Membership Model ==\n")
  cat(paste("Total: ", object$Total, "\t\t K: ",object$K, "\t\t ELBO: ",round(computeELBO(object),2),"\n\n" ,sep = ""))
  df = data.frame(paste("  ",c(1:object$J),sep = ""), object$dist, object$Rj,
                  object$Vj)
  colnames(df) = c("Variable  ", "Variable Type    ", "Replicates    ", "Categories  ")
  print(df, row.names = FALSE, right = FALSE)
}

#' Plot a Mixed Membership Model
#' 
#' Generic S3 function to produce visual representation of a \code{mixedMemModel} object. This function calls either the \code{vizTheta} or
#' the \code{vizMem} function. 
#' 
#' @param x the \code{mixedMemModel} object to be plotted.
#' @param type a string which indicates which estimated parameters to plot; valid options are "theta" or "membership".
#'  \code{vizTheta} is called when the type is "theta" and \code{vizMem} is called when the type is "membership". 
#' @param compare an array or matrix for comparison. When type = "theta", \code{compare} should be an array the same size as x$theta.
#' When type = "membership", \code{compare} should be a matrix the same size as x$phi.
#' @param main the main figure title.
#' @param varNames a vector of strings corresponding to names for each variable if plot type is theta.
#' @param groupNames a vector of strings corresponding to labels for each sub-population.
#' @param nrow the number of rows for the grid of plots.
#' @param ncol the number of columns for the grid of plots. If plot type is "theta", this
#' must be K, if plot type is  "membership", this must be a positive integer.
#' @param indices a vector of integers. If the plot type is  "membership", this indicates which individuals to plot.
#' When plot type is  "theta", this indicates which variables to plot. 
#' @param fitNames a vector of strings corresponding to labels for each fit.
#' @param ... additional parameters to be passed.
#' @seealso mixedMemModel, vizTheta, vizMem
#' @export
plot.mixedMemModel = function(x, type = "theta" , compare = NULL,
                              main = NULL,
                              varNames = NULL,
                              groupNames = NULL,
                              nrow = NULL, ncol = NULL, indices = NULL, fitNames = NULL,...){
  if(type =="theta") {
    if(is.null(main)){
      main = "Estimated Theta"
    }
    vizTheta(x, compare = compare, main = main, varNames = varNames,
             groupNames = groupNames, nrow = nrow, fitNames = fitNames, indices = indices)
  } else if (type == "membership") {
    if(is.null(main)){
      main = "Estimated Memberships"
    }
    vizMem(x,compare = compare, main = main, nrow = nrow, ncol = ncol,
           indices = indices, fitNames = fitNames, groupNames = groupNames)
  }
}
