#' Compute some properties for each triangular element of the mesh
#' 
#' @param mesh A \code{MESH2D} mesh object representing the triangular mesh. This can be created with  \code{\link{create.MESH.2D}}.
#' @return A list with the following variables:
#' \item{\code{detJ}}{A vector of length #triangles. The ith element contains the determinant of the transformation from the reference triangle to the nodes of the i-th triangle. It's values is also the double of the area of each triangle of the basis.}
#' \item{\code{transf}}{A matrix #triangles-by-2-by-2. \code{transf[i,,]} is the 2-by-2 tranformation matrix that transforms the nodes of the reference triangle to the nodes of the i-th triangle.}
#' \item{\code{metric}}{A matrix #triangles-by-2-by-2. \code{metric[i,,]} is the 2-by-2 matrix \cr 
#' \code{transf[i,,]^{-1}*transf[i,,]^{-T}}. This matrix is usuful for the computation
#' of the integrals over the elements of the mesh.} 
#' @description Only executed when the function \code{create.FEM.basis} is run with the option \code{CPP_CODE} = \code{FALSE}. 
#' It computes some quantities associated to the linear map that transforms the ith triangle in the reference triangular element. 
#' These are used for the computation of the integrals necessary to build the mass and stiffness matrix.
#' @usage R_elementProperties(mesh)

R_elementProperties=function(mesh)
{
  nele = dim(mesh$triangles)[[1]]
  
  detJ   = matrix(0,nele,1)      #  vector of determinant of transformations
  metric = array(0,c(nele,2,2))  #  3-d array of metric matrices
  transf = array(0,c(nele,2,2))
  
  for (i in 1:nele)
  {
    diff1x = mesh$nodes[mesh$triangles[i,2],1] - mesh$nodes[mesh$triangles[i,1],1]
    diff1y = mesh$nodes[mesh$triangles[i,2],2] - mesh$nodes[mesh$triangles[i,1],2]
    diff2x = mesh$nodes[mesh$triangles[i,3],1] - mesh$nodes[mesh$triangles[i,1],1]
    diff2y = mesh$nodes[mesh$triangles[i,3],2] - mesh$nodes[mesh$triangles[i,1],2]
    
    transf[i,,] = rbind(cbind(diff1x,diff2x),c(diff1y,diff2y))
    #  Jacobian or double of the area of triangle
    detJ[i] = diff1x*diff2y - diff2x*diff1y
    
    #  Compute controvariant transformation matrix OSS: This is (tranf)^(-T)
    Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/detJ[i]
    
    #  Compute metric matrix
    metric[i,,] = t(Ael)%*%Ael
  } 
  
  FEStruct <- list(detJ=detJ, metric=metric, transf=transf)
  return(FEStruct)
}

#' Compute the mass matrix
#' 
#' @param FEMbasis A \code{FEM} object representing the Finite Element basis. See \code{\link{create.FEM.basis}}.
#' @return A square matrix with the integrals of all the basis' functions pairwise products.
#' The dimension of the matrix is equal to the number of the nodes of the mesh.
#' @description Only executed when \code{smooth.FEM.basis} is run with the option  \code{CPP_CODE} = \code{FALSE}. It computes the mass matrix. The element (i,j) of this matrix contains the integral over the domain of the product between the ith and kth element 
#' of the Finite Element basis. As common practise in Finite Element Analysis, this quantities are computed iterating over all the mesh triangles. 
#' @usage R_mass(FEMbasis)
#' @seealso \code{\link{R_stiff}}

R_mass=function(FEMbasis)
{
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles
  detJ = FEMbasis$detJ
  order = FEMbasis$order
  
  nele  = dim(triangles)[1]
  nnod  = dim(nodes)[1]
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  K0M = NULL
  
  if (order ==2)
  {   
    #  the integrals of products of basis functions for master element:
    
    K0M = matrix(c( 6, -1, -1, -4,  0,  0,
                    -1,  6, -1,  0, -4,  0,
                    -1, -1,  6,  0,  0, -4,
                    -4,  0,  0, 32, 16, 16,
                    0, -4,  0, 16, 32, 16,
                    0,  0, -4, 16, 16, 32), ncol=6, nrow=6, byrow=T)/360
  }else if (order == 1)
  {  
    #  the integrals of products of basis functions for master element:
    K0M = matrix( c( 2,  1,  1,
                     1,  2,  1,
                     1 , 1,  2), ncol=3, nrow=3, byrow=T) / 24
    
    
  }
  
  # assemble the mass matrix
  K0 = matrix(0,nrow=nnod,ncol=nnod)
  for (el in 1:nele)
  {  
    ind = triangles[el,]
    K0[ind,ind] = K0[ind,ind] + K0M * detJ[el]
  }
  
  K0
}

#' Compute the stiffness matrix
#' 
#' @param FEMbasis A \code{FEMbasis} object representing the basis; See \code{\link{create.FEM.basis}}.
#' @return A square matrix with the integrals of all the basis functions' gradients pairwise dot products.
#' The dimension of the matrix is equal to the number of the nodes of the mesh.
#' @description Only executed when \code{smooth.FEM.basis} is run with the option  \code{CPP_CODE} = \code{FALSE}. It computes the stifness matrix. The element (i,j) of this matrix contains the integral over the domain of the scalar product between the gradient of the ith and kth element 
#' of the Finite Element basis. As common practise in Finite Element Analysis, this quantities are computed iterating over all the mesh triangles. 
#' @usage R_stiff(FEMbasis)
#' @seealso \code{\link{R_mass}}

R_stiff= function(FEMbasis)
{
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles
  
  nele  = dim(triangles)[[1]]
  nnod  = dim(nodes)[[1]]
  detJ     = FEMbasis$detJ
  order = FEMbasis$order
  metric = FEMbasis$metric
  
  
  
  KXX = NULL
  KXY = NULL
  KYY = NULL
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  if (order == 2)
  {   
    #  values of K1 for master elements
    KXX = matrix( c( 3, 1, 0, 0, 0,-4, 
                     1, 3, 0, 0, 0,-4,
                     0, 0, 0, 0, 0, 0,
                     0, 0, 0, 8,-8, 0,
                     0, 0, 0,-8, 8, 0,
                     -4,-4, 0, 0, 0, 8), ncol=6, nrow=6, byrow=T)/6
    
    KXY = matrix(c( 3, 0, 1, 0,-4, 0, 
                    1, 0,-1, 4, 0,-4, 
                    0, 0, 0, 0, 0, 0,
                    0, 0, 4, 4,-4,-4,
                    0, 0,-4,-4, 4, 4,
                    -4, 0, 0,-4, 4, 4), ncol=6, nrow=6, byrow=T)/6
    
    KYY = matrix( c( 3, 0, 1, 0,-4, 0,
                     0, 0, 0, 0, 0, 0,
                     1, 0, 3, 0,-4, 0,
                     0, 0, 0, 8, 0,-8,
                     -4, 0,-4, 0, 8, 0,
                     0, 0, 0,-8, 0, 8), ncol=6, nrow=6, byrow=T)/6
  }else if (order == 1)
  {   
    KXX = matrix( c(  1, -1,  0,
                      -1,  1,  0,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KXY = matrix( c(  1,  0, -1,
                      -1,  0,  1,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KYY = matrix( c(  1,  0, -1,
                      0,  0,  0,
                      -1,  0,  1), ncol=3, nrow=3, byrow=T) /2
  }
  
  #  assemble the stiffness matrix
  K1   = matrix(0,nrow=nnod,ncol=nnod)
  for (el in 1:nele)
  {
    ind    = triangles[el,]
    K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY +
             metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
    K1[ind,ind] = K1[ind,ind] + K1M*detJ[el]
  }
  
  K1
}

#' Spatial regression with differential regularization (fully implemented in R code)
#' 
#' @param observations A #observations vector with the observed data values over the domain. The locations of the observations can be specified with the \code{locations} argument.
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table nodes of the mesh. In this last
#' case, an \code{NA} value in the observations vector indicates that there is no observation associated to the corresponding node.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observations in the vector \code{observations}.
#' @param FEMbasis A F\code{EMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @return A list with the following quantities:
#'    \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'    \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#'    \item{\code{beta}}{If covariates is not \code{NULL}, a vector of length #covariates with the regression coefficients associated with each covariate.}
#'    \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @usage R_smooth.FEM.basis(locations, observations, FEMbasis, lambda, covariates, GCV)
#' @seealso \code{\link{smooth.FEM.basis}}, \code{\link{smooth.FEM.PDE.basis}}, \code{\link{smooth.FEM.PDE.sv.basis}}
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp.681.703.

R_smooth.FEM.basis = function(locations, observations, FEMbasis, lambda, covariates = NULL, GCV)
{
  
  # Stores the number of nodes of the mesh. This corresponds to the number of elements of the FE basis.
  numnodes = nrow(FEMbasis$mesh$nodes)
  if(!is.null(covariates))
  {
    covariates = as.matrix(covariates)
  }
  
  if(!is.null(locations))
  {
    locations <- as.matrix(locations)
  }
  
  #  ---------------------------------------------------------------
  # construct mass matrix K0 
  #  ---------------------------------------------------------------
  
  K0 = R_mass(FEMbasis)
  
  #  ---------------------------------------------------------------
  # construct stiffness matrix K1
  #  ---------------------------------------------------------------
  
  K1 = R_stiff(FEMbasis)
  
  
  #  ---------------------------------------------------------------
  # construct the penalty matrix P with ones on diagonal at data points
  #  ---------------------------------------------------------------
  
  #penalty = numeric(numnodes)
  #penalty[data[,1]] = 1
  #P = diag(penalty,nrow=numnodes)
  
  basismat = NULL
  Lmat = matrix(0,nrow = numnodes, ncol = numnodes)
  H = NULL
  Q = NULL
  
  if(!is.null(locations))
  {
    basismat = R_eval.FEM.basis(FEMbasis, locations)
  } 
  
  if(!is.null(covariates))
  {
    if(!is.null(locations))
    {
      H = covariates %*% ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    }else{
      loc_nodes = (1:length(observations))[!is.na(observations)]
      H = covariates[loc_nodes,] %*% ( solve( t(covariates[loc_nodes,]) %*% covariates[loc_nodes,] ) ) %*% t(covariates[loc_nodes,])
    }
    Q=matrix(0,dim(H)[1], dim(H)[2])
    Q = diag(1,dim(H)[1])- H
  }
  
  if(!is.null(locations) && is.null(covariates))
  {
    Lmat = t(basismat) %*% basismat
  }
  if(!is.null(locations) && !is.null(covariates))
  {
    Lmat = t(basismat) %*% Q %*% basismat
  }
  if(is.null(locations) && is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    Lmat[loc_nodes,loc_nodes]=diag(1,length(observations[loc_nodes]))
  }
  if(is.null(locations) && !is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    Lmat[loc_nodes,loc_nodes] = Q
  }
  
  #  ---------------------------------------------------------------
  # construct vector b for system Ax=b
  #  ---------------------------------------------------------------
  
  b = matrix(numeric(numnodes*2),ncol=1)
  if(!is.null(locations) && is.null(covariates))
  {
    b[1:numnodes,] = t(basismat) %*% observations
  }
  if(!is.null(locations) && !is.null(covariates))
  {
    b[1:numnodes,] = t(basismat) %*% Q %*% observations
  }
  if(is.null(locations) && is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    b[loc_nodes,] = observations[loc_nodes]
  }
  if(is.null(locations) && !is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    b[loc_nodes,] = Q %*% observations[loc_nodes]
  }
  #print(b)
  
  #  ---------------------------------------------------------------
  # construct matrix A for system Ax=b.
  #  ---------------------------------------------------------------
  
  bigsol = list(solutions = matrix(nrow = 2*numnodes, ncol = length(lambda)), edf = vector(mode="numeric",length = length(lambda)))
  for(i in 1: length(lambda))
  {
    A  = rbind(cbind(Lmat, -lambda[i]*K1), cbind(K1, K0))
    #print(A)
    # solve system
    bigsol[[1]][,i] = solve(A,b)
    if(GCV == TRUE)
    {
      S = solve(Lmat + lambda[i] * K1 %*% solve(K0) %*% K1)
      if(is.null(locations))
      {
        loc_nodes = (1:length(observations))[!is.na(observations)]
        if(is.null(covariates))
        {
          bigsol[[2]][i]  = sum(diag(S[loc_nodes,loc_nodes]))
        }else{
          bigsol[[2]][i]  = ncol(covariates)+sum(diag(S[loc_nodes,loc_nodes]%*%Q))
        }
      }else{
        if(is.null(covariates))
        {
          bigsol[[2]][i]  = sum(diag(basismat%*%S[1:numnodes,1:numnodes]%*%t(basismat)))
        }else{
          bigsol[[2]][i]  = ncol(covariates)+sum(diag(basismat%*%S[1:numnodes,1:numnodes]%*%t(basismat)%*%Q))
        }
      }
    }else{
      bigsol[[2]][i]  = NA
    }
  }
  return(bigsol)
}

#' Evaluate Finite Element bases and their Derivatives at a set of locations
#' 
#' @param locations A 2-columns matrix with the spatial locations where the bases should be evaluated.
#' @param FEMbasis An \code{FEMbasis} object representing the Finite Element basis; See \link{create.FEM.basis}.
#' @param nderivs A vector of lenght 2 specifying the order of the partial derivatives of the bases to be evaluated. The vectors' entries can
#' be 0,1 or 2, where 0 indicates that only the basis functions, and not their derivatives, should be evaluated.
#' @return 
#' A matrix of basis function values. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' basis function evaluated 
#' @description Only executed when the function \code{smooth.FEM.basis} is run with the option \code{CPP_CODE} = \code{FALSE}. It evaluates the Finite Element basis functions and their derivatives up to order 2 at the specified set of locations. 
#' This version of the function is implemented using only R code. It is called by \link{R_smooth.FEM.basis}.
#' @usage R_eval.FEM.basis(FEMbasis, locations, nderivs = matrix(0,1,2))
#' @seealso \code{\link{R_eval.FEM}}

R_eval.FEM.basis <- function(FEMbasis, locations, nderivs = matrix(0,1,2))
{ 
  if(length(nderivs) != 2)
  {
    stop('NDERIVS not of length 2.')
  }
  
  if(sum(nderivs)>2)
  {
    stop('Maximum derivative order is greater than two.')
  }
  
  N = nrow(locations)
  nbasis = FEMbasis$nbasis
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  mesh = FEMbasis$mesh
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  order = FEMbasis$order
  #nodeindex = params$nodeindex
  detJ = FEMbasis$detJ
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[,1],]
  v2 = nodes[triangles[,2],]
  v3 = nodes[triangles[,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = FEMbasis$detJ
  ones3 = matrix(1,3,1)
  modJMat = modJ %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = R_insideIndex(mesh, locations[i,])
  }
  
  evalmat = matrix(0,N,nbasis)
  
  for(i in 1:N)
  {
    indi = ind[i]
    
    if(!is.nan(indi))
    {
      baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
      baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
      baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,triangles[indi,1]] = 2* baryc1^2 - baryc1
          evalmat[i,triangles[indi,2]] = 2* baryc2^2 - baryc2
          evalmat[i,triangles[indi,3]] = 2* baryc3^2 - baryc3
          evalmat[i,triangles[indi,6]] = 4* baryc1 * baryc2
          evalmat[i,triangles[indi,4]] = 4* baryc2 * baryc3
          evalmat[i,triangles[indi,5]] = 4* baryc3 * baryc1
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = (4* baryc1 - 1) * M1[indi,2]
          evalmat[i,triangles[indi,2]] = (4* baryc2 - 1) * M2[indi,2]
          evalmat[i,triangles[indi,3]] = (4* baryc3 - 1) * M3[indi,2]
          evalmat[i,triangles[indi,6]] = (4* baryc2 ) * M1[indi,2] + 4*baryc1 * M2[indi,2]
          evalmat[i,triangles[indi,4]] = (4* baryc3 ) * M2[indi,2] + 4*baryc2 * M3[indi,2]
          evalmat[i,triangles[indi,5]] = (4* baryc1 ) * M3[indi,2] + 4*baryc3 * M1[indi,2]
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = (4*baryc1 - 1)*M1[indi,3]
          evalmat[i,triangles[indi,2]] = (4*baryc2 - 1)*M2[indi,3]
          evalmat[i,triangles[indi,3]] = (4*baryc3 - 1)*M3[indi,3]
          evalmat[i,triangles[indi,6]] = 4*baryc2*M1[indi,3] + 4*baryc1*M2[indi,3]
          evalmat[i,triangles[indi,4]] = 4*baryc3*M2[indi,3] + 4*baryc2*M3[indi,3]
          evalmat[i,triangles[indi,5]] = 4*baryc1*M3[indi,3] + 4*baryc3*M1[indi,3]
        }
        else if(nderivs[1] == 1 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,2]%*%M1[indi,3];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,2]%*%M2[indi,3];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,2]%*%M3[indi,3];
          evalmat[i,triangles[indi,6]] = 4*M2[indi,2]%*%M1[indi,3] + 4*M2[indi,3]%*%M1[indi,2];
          evalmat[i,triangles[indi,4]] = 4*M3[indi,2]%*%M2[indi,3] + 4*M3[indi,3]%*%M2[indi,2];
          evalmat[i,triangles[indi,5]] = 4*M1[indi,2]%*%M3[indi,3] + 4*M1[indi,3]%*%M3[indi,2];
        }
        else if(nderivs[1] == 2 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,2]%*%M1[indi,2];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,2]%*%M2[indi,2];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,2]%*%M3[indi,2];
          evalmat[i,triangles[indi,6]] = 8*M2[indi,2]%*%M1[indi,2];
          evalmat[i,triangles[indi,4]] = 8*M3[indi,2]%*%M2[indi,2];
          evalmat[i,triangles[indi,5]] = 8*M1[indi,2]%*%M3[indi,2];
        }
        else if(nderivs[1] == 0 && nderivs[2] == 2)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,3]%*%M1[indi,3];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,3]%*%M2[indi,3];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,3]%*%M3[indi,3];
          evalmat[i,triangles[indi,6]] = 8*M2[indi,3]%*%M1[indi,3];
          evalmat[i,triangles[indi,4]] = 8*M3[indi,3]%*%M2[indi,3];
          evalmat[i,triangles[indi,5]] = 8*M1[indi,3]%*%M3[indi,3];
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
      else
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,triangles[indi,1]] = baryc1;
          evalmat[i,triangles[indi,2]] = baryc2;
          evalmat[i,triangles[indi,3]] = baryc3;
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = M1[indi[1],2];
          evalmat[i,triangles[indi,2]] = M1[indi[2],2];
          evalmat[i,triangles[indi,3]] = M1[indi[3],2]; 
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = M1[indi[1],3];
          evalmat[i,triangles[indi,2]] = M1[indi[2],3];
          evalmat[i,triangles[indi,3]] = M1[indi[3],3]; 
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
    }
  }
  return(evalmat)
}

#' Evaluate a FEM object at a set of locations
#' 
#' @param FEM A \code{FEM} object to be evaluated
#' @param locations A 2-columns matrix with the spatial locations where the FEM object should be evaluated
#' @return 
#' A matrix of numeric evaluations of the \code{FEM} object. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' function evaluated.
#' @description Only executed when the function \code{smooth.FEM.basis} is run with the option \code{CPP_CODE} = \code{FALSE}. It evaluates a FEM object at the specified set of locations. 
#' @usage R_eval.FEM(FEM, locations)
#' @seealso \code{\link{R_eval.FEM.basis}}

R_eval.FEM <- function(FEM, locations)
{ 
  if (is.vector(locations))
  {
    locations = t(as.matrix(locations))
  }
  
  N = nrow(locations)
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  FEMbasis = FEM$FEMbasis
  mesh = FEMbasis$mesh
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  coeff = FEM$coeff
  nsurf = dim(coeff)[2]
  
  FEMbasis = FEM$FEMbasis
  order = FEMbasis$order
  #nodeindex = params$nodeindex
  detJ = FEMbasis$detJ
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[,1],]
  v2 = nodes[triangles[,2],]
  v3 = nodes[triangles[,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = FEMbasis$detJ
  ones3 = matrix(1,3,1)
  modJMat = modJ %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = R_insideIndex(mesh, as.numeric(locations[i,]))
  }
  
  evalmat = matrix(NA, nrow=N, ncol=nsurf)
  
  for (isurf in 1:nsurf)
  {
    for(i in 1:N)
    {
      indi = ind[i]
      
      if(!is.nan(indi))
      {
        baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
        baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
        baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
        
        if(order == 2)
        {
          c1 = coeff[triangles[indi,1],isurf]
          c2 = coeff[triangles[indi,2],isurf]
          c3 = coeff[triangles[indi,3],isurf]
          c4 = coeff[triangles[indi,6],isurf]
          c5 = coeff[triangles[indi,4],isurf]
          c6 = coeff[triangles[indi,5],isurf]
          
          fval =  c1*(2* baryc1^2 - baryc1) +
            c2*(2* baryc2^2 - baryc2) +
            c3*(2* baryc3^2 - baryc3) +
            c4*(4* baryc1 * baryc2) +
            c5*(4* baryc2 * baryc3) +
            c6*(4* baryc3 * baryc1)
          evalmat[i,isurf] = fval
        }
        else
        {
          c1 = coeff[triangles[indi,1],isurf]
          c2 = coeff[triangles[indi,2],isurf]
          c3 = coeff[triangles[indi,3],isurf]
          fval = c1*baryc1 + c2*baryc2 + c3*baryc3
          evalmat[i,isurf] = fval
        }
      }
    }
  }
  return(evalmat)
}

R_tricoefCal = function(mesh)
{
  #  TRICOEFCAL compute the coefficient matrix TRICOEF
  #  required to test of a point is indside a triangle
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  ntri   = dim(triangles)[[1]]
  
  #  compute coefficients for computing barycentric coordinates if
  #  needed
  
  tricoef = matrix(0, nrow=ntri, ncol=4)
  tricoef[,1] = nodes[triangles[,1],1]-nodes[triangles[,3],1]
  tricoef[,2] = nodes[triangles[,2],1]-nodes[triangles[,3],1]
  tricoef[,3] = nodes[triangles[,1],2]-nodes[triangles[,3],2]
  tricoef[,4] = nodes[triangles[,2],2]-nodes[triangles[,3],2]
  detT = matrix((tricoef[,1]*tricoef[,4] - tricoef[,2]*tricoef[,3]),ncol=1)
  tricoef = tricoef/(detT %*% matrix(1,nrow=1,ncol=4))
  
  return(tricoef)
}

R_insideIndex = function (mesh, location)
{
  #  insideIndex returns the index of the triangle containing the point
  # (X,Y) if such a triangle exists, and NaN otherwise.
  #  TRICOEF may have already been calculated for efficiency,
  #  but if the function is called with four arguments, it is calculated.
  
  
  eps=2.2204e-016
  small = 10000*eps
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  X = location[1]
  Y = location[2]
  
  ntri   = dim(triangles)[[1]]
  indtri   = matrix(1:ntri,ncol=1)
  
  #  compute coefficients for computing barycentric coordinates if needed
  
  tricoef = R_tricoefCal(mesh)
  
  #  compute barycentric coordinates
  r3 = X - nodes[triangles[,3],1]
  s3 = Y - nodes[triangles[,3],2]
  lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
  lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
  lam3 = 1 - lam1 - lam2
  
  #  test these coordinates for a triple that are all between 0 and 1
  int  = (-small <= lam1 & lam1 <= 1+small) & 
    (-small <= lam2 & lam2 <= 1+small) & 
    (-small <= lam3 & lam3 <= 1+small)
  
  #  return the index of this triple, or NaN if it doesn't exist
  indi = indtri[int]
  if (length(indi)<1)
  {
    ind = NA
  }else{
    ind = min(indi)
  }
  
  ind
}

R_eval_local.FEM = function(FEM, locations, element_index)
{
  N = nrow(locations)
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  mesh = FEMbasis$mesh
  nodes = mesh$nodes
  triangles = mesh$triangles
  coeff = FEM$coeff
  nsurf = dim(coeff)[2]
  
  FEMbasis = FEM$FEMbasis
  order = FEMbasis$order
  #nodeindex = params$nodeindex
  detJ = FEMbasis$detJ
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[element_index,1],]
  v2 = nodes[triangles[element_index,2],]
  v3 = nodes[triangles[element_index,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = FEMbasis$detJ[element_index]
  ones3 = matrix(1,3,1)
  #modJMat = modJ %*% t(ones3)
  
  M1 = c(v2[1]*v3[2] - v3[1]*v2[2], v2[2] - v3[2], v3[1] - v2[1])/(modJ)
  M2 = c(v3[1]*v1[2] - v1[1]*v3[2], v3[2] - v1[2], v1[1] - v3[1])/(modJ)
  M3 = c(v1[1]*v2[2] - v2[1]*v1[2], v1[2] - v2[2], v2[1] - v1[1])/(modJ)
  
  evalmat = matrix(NA, nrow=N, ncol=nsurf)
  
  for (isurf in 1:nsurf)
  {
    for(i in 1:N)
    {
      baryc1 = (M1*Pgpts[i,]) %*% ones3
      baryc2 = (M2*Pgpts[i,]) %*% ones3
      baryc3 = (M3*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        c1 = coeff[triangles[element_index,1],isurf]
        c2 = coeff[triangles[element_index,2],isurf]
        c3 = coeff[triangles[element_index,3],isurf]
        c4 = coeff[triangles[element_index,6],isurf]
        c5 = coeff[triangles[element_index,4],isurf]
        c6 = coeff[triangles[element_index,5],isurf]
        
        fval =  c1*(2* baryc1^2 - baryc1) +
          c2*(2* baryc2^2 - baryc2) +
          c3*(2* baryc3^2 - baryc3) +
          c4*(4* baryc1 * baryc2) +
          c5*(4* baryc2 * baryc3) +
          c6*(4* baryc3 * baryc1)
        evalmat[i,isurf] = fval
      }else{
        c1 = coeff[triangles[element_index,1],isurf]
        c2 = coeff[triangles[element_index,2],isurf]
        c3 = coeff[triangles[element_index,3],isurf]
        fval = c1*baryc1 + c2*baryc2 + c3*baryc3
        evalmat[i,isurf] = fval
      }
    }
  }
  return(evalmat)
}

R_plot.ORD1.FEM = function(FEM, ...)  
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
  # vectors X and Y;
  #
  
  
  
  #   if (!is.fd(fdobj))
  #   {
  #     stop('FDOBJ is not an FD object')
  #   }
  
  nodes = FEM$FEMbasis$mesh$nodes
  triangles = FEM$FEMbasis$mesh$triangles
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$basis
  
  mesh = FEMbasis$mesh
  
  heat = heat.colors(100)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    axes3d()
    
    z = coeff[as.vector(t(triangles)),isurf]
    rgl.triangles(x = nodes[as.vector(t(triangles)) ,1], y = nodes[as.vector(t(triangles)) ,2], 
                  z=coeff[as.vector(t(triangles)),isurf], 
                  color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1],...)
    aspect3d(2,2,1)
    rgl.viewpoint(0,-45)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_plot.ORDN.FEM = function(FEM, num_refinements, ...)  
{
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  heat = heat.colors(100)
  
  coeff = FEM$coeff
  
  # num_refinements sets the number od division on each triangle edge to be applied for rifenment
  if(is.null(num_refinements))
  {
    num_refinements = 10
  }
  
  # For the reference triangles we construct a regular mesh
  x = seq(from = 0, to = 1, length.out = num_refinements+1)
  y = seq(from = 0, to = 1, length.out = num_refinements+1)
  points_ref = expand.grid(x,y)
  points_ref = points_ref[points_ref[,1] + points_ref[,2] <= 1,]
  
  
  meshi = create.MESH.2D(nodes = points_ref, order = 1)   
  #plot(meshi)
  
  # locations is the matrix with that will contain the coordinate of the points where the function is 
  # evaluated (1st and 2nd column) and the columns with the evaluation of the ith fucntion on that point
  
  locations = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$nodes), ncol = 2+ncol(coeff))
  triangles = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$triangles), ncol = 3)
  tot = 0
  
  
  for (i in 1:nrow(mesh$triangles))
  {
    # For each traingle we define a fine mesh as the transofrmation of the one constructed for the reference
    pointsi = t(FEMbasis$transf[i,,]%*%t(meshi$nodes) + mesh$nodes[mesh$triangles[i,1],])
    #We evaluate the fine mesh OBS: we know the triangle we are working on no need for point location
    z = R_eval_local.FEM(FEM, locations = pointsi, element_index = i)
    
    #We store the results
    locations[((i-1)*nrow(pointsi)+1):(i*nrow(pointsi)),] = cbind(pointsi,z)
    triangles[((i-1)*nrow(meshi$triangles)+1):(i*nrow(meshi$triangles)),] = meshi$triangles+tot
    tot = tot + nrow(meshi$nodes)
  }
  
  heat = heat.colors(100)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    z = locations[as.vector(t(triangles)), 2 + isurf]
    rgl.triangles(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2], 
                  z = z, 
                  color = heat[round(99*(z-min(z))/(max(z)-min(z)))+1],...)
    aspect3d(2,2,1)
    rgl.viewpoint(0,-45)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_image.ORD1.FEM = function(FEM)  
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
  # vectors X and Y;
  #
  
  #   if (!is.fd(fdobj))
  #   {
  #     stop('FDOBJ is not an FD object')
  #   }
  
  nodes = FEM$FEMbasis$mesh$nodes
  triangles = FEM$FEMbasis$mesh$triangles
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  heat = heat.colors(100)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    #rgl.open()
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    z = coeff[as.vector(t(triangles)),isurf]
    rgl.triangles(x = nodes[as.vector(t(triangles)) ,1], y = nodes[as.vector(t(triangles)) ,2], 
                  z=0, 
                  color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,0)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_image.ORDN.FEM = function(FEM, num_refinements)  
{
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  heat = heat.colors(100)
  
  coeff = FEM$coeff
  
  if(is.null(num_refinements))
  {
    num_refinements = 10
  }
  
  x = seq(from = 0, to = 1, length.out = num_refinements+1)
  y = seq(from = 0, to = 1, length.out = num_refinements+1)
  points_ref = expand.grid(x,y)
  points_ref = points_ref[points_ref[,1] + points_ref[,2] <= 1,]
  
  meshi = create.MESH.2D(nodes = points_ref, order = 1)   
  #plot(meshi)
  
  locations = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$nodes), ncol = 3*ncol(coeff))
  triangles = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$triangles), ncol = 3*ncol(coeff))
  tot = 0
  
  for (i in 1:nrow(mesh$triangles))
  {
    pointsi = t(FEMbasis$transf[i,,]%*%t(meshi$nodes) + mesh$nodes[mesh$triangles[i,1],])
    z = R_eval_local.FEM(FEM, locations = pointsi, element_index = i)
    
    #mesh3 <- addNormals(subdivision3d(tmesh3d(vertices = t(locations), indices = t(triangles), homogeneous = FALSE)),deform = TRUE)
    
    locations[((i-1)*nrow(pointsi)+1):(i*nrow(pointsi)),] = cbind(pointsi,z)
    triangles[((i-1)*nrow(meshi$triangles)+1):(i*nrow(meshi$triangles)),] = meshi$triangles+tot
    tot = tot + nrow(meshi$nodes)
  }
  
  heat = heat.colors(100)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    z = locations[as.vector(t(triangles)), 2 + isurf];
    rgl.triangles(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2], 
                  z=0, 
                  color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,0)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

