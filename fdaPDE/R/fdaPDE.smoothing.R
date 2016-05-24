#' Spatial regression with differential regularization: stationary and isotropic case (Laplacian)
#' 
#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates \code{x} and \code{y} of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of the algorithm. This usually ensures a much faster computation.
#' @return A list with the following variables:
#' \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#' \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#' \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#' the smoothing parameter is equal to \code{lambda[j]}.}
#' \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; isotropic and stationary case. In particular, the regularizing term involves the Laplacian of the spatial field. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.basis(locations = NULL, observations, FEMbasis, lambda, 
#'        covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
#' @seealso \code{\link{smooth.FEM.PDE.basis}}, \code{\link{smooth.FEM.PDE.sv.basis}}
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
#' @examples
#' library(fdaPDE)
#' ## Load the Meuse data and a domain boundary for these data
#' data(MeuseData)
#' data(MeuseBorder)
#' ## Create a triangular mesh for these data with the provided boundary and plot it
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' plot(mesh)
#' ## Create the Finite Element basis 
#' FEMbasis = create.FEM.basis(mesh)
#' ## Estimate zync field without using covariates, setting the smoothing parameter to 10^3.5
#' data = log(MeuseData[,"zinc"])
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, 
#'                              FEMbasis = FEMbasis, lambda = lambda)
#' ## Plot the estimated spatial field 
#' plot(ZincMeuse$fit.FEM)
#' # Now repeat the analysis using as covariates the square root of the log-distance 
#' # from river \code{sqrt(dist.log(m))} and the altitude \code{elev}
#' desmat = matrix(1,nrow=nrow(MeuseData),ncol=2)
#' desmat[,1] = sqrt(MeuseData[,"dist.log(m)"])
#' desmat[,2] = MeuseData[,"elev"]
#' ZincMeuseCovar = smooth.FEM.basis(observations = data, covariates = desmat, 
#'                                    FEMbasis = FEMbasis, lambda = lambda)
#' # Plot of the non parametric part (f) of the regression model y_i = beta_1 x_i1 + beta_2 x_i2 + f
#' plot(ZincMeuseCovar$fit.FEM)
#' # Print covariates' regression coefficients
#' print(ZincMeuseCovar$beta)


smooth.FEM.basis<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
 
  ##################### Checking parameters, sizes and conversion ################################
  checkSmoothingParameters(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL)
  
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  
  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL)
  ################## End checking parameters, sizes and conversion #############################
  
  
  bigsol = NULL
  if(CPP_CODE == FALSE)
  {
    print('R Code Execution')
    if(!is.null(BC))
    {
      stop('If you want to use Dirichlet boundary conditions, please set CPP_CODE = TRUE')
    }
    
    bigsol = R_smooth.FEM.basis(locations, observations, FEMbasis, lambda, covariates, GCV)   
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.basis(locations, observations, FEMbasis, lambda, covariates, BC, GCV)
  }
  
  numnodes = nrow(FEMbasis$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)  
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates, CPP_CODE)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edf = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}

#' Spatial regression with differential regularization: anysotropic case (elliptic PDE)
#' 
#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates \code{x} and \code{y} of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the parameters of the elliptic PDE in the regularizing term: \code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K; \code{b}, a vector of length 2 of advection coefficients. This induces a 
#' smoothing only in the direction specified by the vector \code{b}; \code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of the algorithm. This usually ensures a much faster computation.
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#'          the smoothing parameter is equal to \code{lambda[j]}.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic case. In particular, the regularizing term involves a second order elliptic PDE, that models the space-variation of the phenomenon. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.PDE.basis(locations = NULL, observations, FEMbasis, 
#'        lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, 
#'        CPP_CODE = TRUE)
#' @seealso \code{\link{smooth.FEM.basis}}, \code{\link{smooth.FEM.PDE.sv.basis}}
#' @references Azzimonti, L., Sangalli, L.M., Secchi, P., Domanin, M., and Nobile, F., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. DOI. 10.1080/01621459.2014.946036. \cr
#'  Azzimonti, L., Nobile, F., Sangalli, L.M., and Secchi, P., 2014. Mixed Finite Elements for Spatial Regression with PDE Penalization. SIAM/ASA Journal on Uncertainty Quantification, 2(1), pp.305-335. 
#' @examples 
#' # Load the mesh and plot it
#' data(mesh.2D.simple)
#' plot(mesh.2D.simple)
#' # Create a vector of noisy samples of an underlying spatial field, 
#' # located over the nodes of the mesh
#' observations = sin(pi*mesh.2D.simple$nodes[,1]) + rnorm(n = nrow(mesh.2D.simple$nodes), sd = 0.1)
#' # Create the FEM basis object
#' FEMbasis = create.FEM.basis(mesh.2D.simple)
#' 
#' # Set a vector of smoothing coefficients
#' lambda = c(10^-4, 1, 10^4)
#' 
#' # Set the anysotropic smoothing matrix K
#' PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
#' # Estimate one field for each smoothing parameter and plot these
#' FEM_CPP_PDE = smooth.FEM.PDE.basis(observations = observations, 
#'                                    FEMbasis = FEMbasis, lambda = lambda, 
#'                                    PDE_parameters = PDE_parameters_anys)
#' plot(FEM_CPP_PDE$fit.FEM)
#' 
#' # Evaluate solution in three points
#' eval.FEM(FEM_CPP_PDE$fit.FEM, locations = rbind(c(0,0),c(0.5,0),c(-2,-2)))
smooth.FEM.PDE.basis<-function(locations = NULL, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  ##################### Checking parameters, sizes and conversion ################################
  checkSmoothingParameters(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = PDE_parameters, PDE_parameters_func = NULL)
  
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  
  if(!is.null(PDE_parameters))
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }
  
  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = PDE_parameters, PDE_parameters_func = NULL)
  ################## End checking parameters, sizes and conversion #############################
  
  
  bigsol = NULL  
  
  if(CPP_CODE == FALSE)
  {
    print('Function implemented only in C++, turn CPP_CODE = TRUE')  
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.basis(locations, observations, FEMbasis, lambda, PDE_parameters, covariates, BC, GCV)
  }
  
  numnodes = nrow(FEMbasis$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)  
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates, CPP_CODE)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edf = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}


#' Spatial regression with differential regularization: anysotropic and non-stationary case (elliptic PDE with space-varying coefficients)
#' 
#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates \code{x} and \code{y} of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the space-varying parameters of the elliptic PDE in the regularizing term: \code{K}, a function that for each spatial location in the spatial domain 
#' (indicated by the vector of the 2 spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a local preferential direction that corresponds to the first eigenvector of the diffusion matrix K.The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' an array with dimensions 2-by-2-by-#points.\code{b}, a function that for each spatial location in the spatial domain returns 
#' a vector of length 2 of transport coefficients. This induces a local smoothing only in the direction specified by the vector \code{b}. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a matrix with dimensions 2-by-#points; \code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{c} induces a shrinkage of the surface to zero. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points; \code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of the algorithm. This usually ensures a much faster computation.
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#'          the smoothing parameter is equal to \code{lambda[j]}.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic and non-stationary case. In particular, the regularizing term involves a second order elliptic PDE with space-varying coefficients, that models the space-variation of the phenomenon. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.PDE.sv.basis(locations = NULL, observations, FEMbasis, 
#'  lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, 
#'  CPP_CODE = TRUE)
#' @seealso \code{\link{smooth.FEM.basis}}, \code{\link{smooth.FEM.PDE.basis}}
#' @references Azzimonti, L., Sangalli, L.M., Secchi, P., Domanin, M., and Nobile, F., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. DOI. 10.1080/01621459.2014.946036. \cr
#'  Azzimonti, L., Nobile, F., Sangalli, L.M., and Secchi, P., 2014. Mixed Finite Elements for Spatial Regression with PDE Penalization. SIAM/ASA Journal on Uncertainty Quantification, 2(1), pp.305-335. 
#' @examples 
#' # Loading the mesh
#' data(mesh.2D.rectangular)
#' # Create the FEM basis object
#' FEMbasis = create.FEM.basis(mesh.2D.rectangular)
#' # Create a vector of noisy samples of an underlying spatial field, 
#' # located over the nodes of the mesh
#' observations = sin(0.2*pi*mesh.2D.rectangular$nodes[,1]) + 
#' rnorm(n = nrow(mesh.2D.rectangular$nodes), sd = 0.1)
#' # Set the smoothing coefficient
#' lambda = c(10^-2)
#' #Set the space vriant coefficients of the penalizying PDE
#' K_func<-function(points)
#' {
#' mat<-c(0.01,0,0,1)
#' output = array(0, c(2, 2, nrow(points)))
#' for (i in 1:nrow(points))
#'    output[,,i] = 0.5*mat %*% t(points[i,1]^2)
#' output
#' }
#' b_func<-function(points)
#' {
#' output = array(0, c(2, nrow(points)))
#' for (i in 1:nrow(points))
#'    output[,i] = 0
#' output
#' }
#' 
#' c_func<-function(points)
#' {
#' rep(c(0), nrow(points))
#' }
#' 
#' u_func<-function(points)
#' {
#' rep(c(0), nrow(points))
#' }
#' # Assemble the parameters in one object
#' PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
#' # Estimate the underlying spatial field and plot these
#' FEM_CPP_PDE = smooth.FEM.PDE.sv.basis(observations = observations, 
#'              FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters)
#' plot(FEM_CPP_PDE$fit.FEM)
smooth.FEM.PDE.sv.basis<-function(locations = NULL, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  ##################### Checking parameters, sizes and conversion ################################
  checkSmoothingParameters(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = NULL, PDE_parameters_func = PDE_parameters)
  
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  
  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = NULL, PDE_parameters_func = PDE_parameters)
  ################## End checking parameters, sizes and conversion #############################
  
  bigsol = NULL  
  
  if(CPP_CODE == FALSE)
  {
    print('Function implemented only in C++, turn CPP_CODE = TRUE')  
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.sv.basis(locations, observations, FEMbasis, lambda, PDE_parameters, covariates, BC, GCV)
  }
  
  numnodes = nrow(FEMbasis$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)  
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates, CPP_CODE)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edf = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}


getBetaCoefficients<-function(locations, observations, fit.FEM, covariates, CPP_CODE = FALSE)
{
  loc_nodes = NULL
  fnhat = NULL
  betahat = NULL
  
  if(!is.null(covariates))
  {
    if(is.null(locations))
    {
      loc_nodes = (1:length(observations))[!is.na(observations)]
      fnhat = as.matrix(fit.FEM$coeff[loc_nodes,])
    }else{
      loc_nodes = 1:length(observations)
      fnhat = eval.FEM(FEM = fit.FEM, locations = locations, CPP_CODE)
    }
    ## #row number of covariates, #col number of functions
    betahat = matrix(0, nrow = ncol(covariates), ncol = ncol(fnhat))
    for(i in 1:ncol(fnhat))
      betahat[,i] = as.vector(lm.fit(covariates,as.vector(observations-fnhat[,i]))$coefficients)
  }
  
 return(betahat)
}


getGCV<-function(locations, observations, fit.FEM, covariates = NULL, edf)
{
  loc_nodes = NULL
  fnhat = NULL
  
  edf = as.vector(edf)
    
  if(is.null(locations))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    fnhat = as.matrix(fit.FEM$coeff[loc_nodes,])
  }else{
    loc_nodes = 1:length(observations)
    fnhat = eval.FEM(FEM = fit.FEM, locations = locations, CPP_CODE = FALSE)
  }
  
  zhat = NULL
  zhat = matrix(nrow = length(loc_nodes), ncol = length(edf))
  if(!is.null(covariates))
  {
    desmatprod = ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    for ( i in 1:length(edf))
    {
      betahat  = desmatprod %*% (observations-fnhat[,i])
      zhat[,i]     = covariates %*% betahat + fnhat[,i]
    }
  }else{
    zhat = fnhat
  }
  
  np = length(loc_nodes)
  
  stderr2 = numeric(length(edf))
  GCV       = numeric(length(edf))
  
  zhat <- as.matrix(zhat)
  
  if(any(np - edf <= 0))
  {
    warning("Some values of 'edf' are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.")  
  }
  
  for (i in 1:length(edf))
  {
    stderr2[i] = t(observations[loc_nodes] - zhat[,i]) %*% (observations[loc_nodes] - zhat[,i]) / ( np - edf[i] )
    GCV[i] = ( np / ( np - edf[i] )) * stderr2[i]
  }
  
  # NA if stderr2 is negative
  stderr = vector('numeric', length(stderr2));
  stderr[stderr2>=0] = sqrt(stderr2[stderr2>=0]);
  stderr[stderr2<0] = NaN;
  
  return(list(stderr = stderr, GCV = GCV))
}
