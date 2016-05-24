#dyn.load("../Release/fdaPDE.so")

CPP_smooth.FEM.basis<-function(locations, observations, FEMbasis, lambda, covariates = NULL, BC = NULL, GCV)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  ##TO BE CHANGED SOON: LOW PERFORMANCES, IMPLIES COPY OF PARAMETERS
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates = as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(lambda)<- "double"
  storage.mode(BC$BC_indices)<- "integer"
  storage.mode(BC$BC_values)<-"double"
  
  GCV = as.integer(GCV)
  storage.mode(GCV)<-"integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, observations, FEMbasis$mesh, 
                  FEMbasis$order, lambda, covariates,
                  BC$BC_indices, BC$BC_values, GCV,
                  package = "fdaPDE")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}

CPP_smooth.FEM.PDE.basis<-function(locations, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation  
  ##TO BE CHANGED SOON: LOW PERFORMANCES, IMPLIES COPY OF PARAMETERS
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  ## Set propr type for correct C++ reading
  
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates = as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(lambda)<- "double"
  storage.mode(BC$BC_indices)<- "integer"
  storage.mode(BC$BC_values)<-"double"
  storage.mode(GCV)<-"integer"
  
  storage.mode(PDE_parameters$K)<-"double"
  storage.mode(PDE_parameters$b)<-"double"
  storage.mode(PDE_parameters$c)<-"double"
  
  ## Call C++ function
  bigsol <- .Call("regression_PDE", locations, observations, FEMbasis$mesh, 
                  FEMbasis$order, lambda, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariates,
                  BC$BC_indices, BC$BC_values, GCV,
                  package = "fdaPDE")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}

CPP_smooth.FEM.PDE.sv.basis<-function(locations, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV)
{
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  ##TO BE CHANGED SOON: LOW PERFORMANCES, IMPLIES COPY OF PARAMETERS
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  
  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = FEMbasis$mesh, order = FEMbasis$order),ncol = 2)
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$b = (PDE_parameters$b)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)
  
  #print(PDE_param_eval)
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates = as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(lambda)<- "double"
  storage.mode(BC$BC_indices)<- "integer"
  storage.mode(BC$BC_values)<-"double"
  storage.mode(GCV)<-"integer"
  
  storage.mode(PDE_param_eval$K)<-"double"
  storage.mode(PDE_param_eval$b)<-"double"
  storage.mode(PDE_param_eval$c)<-"double"
  storage.mode(PDE_param_eval$u)<-"double"
  
  ## Call C++ function
  bigsol <- .Call("regression_PDE_space_varying", locations, observations, FEMbasis$mesh, 
                  FEMbasis$order, lambda, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u, covariates,
                  BC$BC_indices, BC$BC_values, GCV,
                  package = "fdaPDE")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}

CPP_eval.FEM = function(FEM, locations, redundancy)
{
  FEMbasis = FEM$FEMbasis
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  ##TO BE CHANGED SOON: LOW PERFORMANCES, IMPLIES COPY OF PARAMETERS
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  # Imposing types, this is necessary for correct reading from C++
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff = as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,nrow(locations),ncol(coeff))
  for (i in 1:ncol(coeff))
  {
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations[,1], locations[,2], coeff[,i], 
                   FEMbasis$order, redundancy,
                   package = "fdaPDE")
  }
  #Returning the evaluation matrix
  evalmat
}

CPP_get_evaluations_points = function(mesh, order)
{
  # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
  #
  #        arguments:
  # X         an array of x-coordinates.
  # Y         an array of y-coordinates.
  # FELSPLOBJ a FELspline object
  # FAST      a boolean indicating if the walking algorithm should be apply 
  #        output:
  # EVALMAT   an array of the same size as X and Y containing the value of 
  #           FELSPLOBJ at (X,Y).
  
  ## C++ indexes starts from zero, needed conversion.
  mesh$triangles = mesh$triangles - 1
  mesh$edges = mesh$edges - 1
  mesh$neighbors[mesh$neighbors != -1] = mesh$neighbors[mesh$neighbors != -1] - 1
  
  # Imposing types, this is necessary for correct reading from C++
  storage.mode(mesh$points) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  storage.mode(mesh$edges) <- "integer"
  storage.mode(mesh$neighbors) <- "integer"
  storage.mode(order) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  points <- .Call("get_integration_points",mesh, order,
                  package = "fdaPDE")
  
  #Returning the evaluation matrix
  points
}