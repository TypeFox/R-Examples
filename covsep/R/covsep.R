#' covsep: tests for determining if the covariance structure of
#' 2-dimensional data is separable
#'
#' Functions for testing if the covariance structure of 2-dimensional data (e.g.
#' samples of surfaces X_i = X_i(s,t)) is separable, i.e. if cov(X) = C_1 x C_2. 
#' A complete descriptions of the implemented tests can be found in the paper
#' arXiv:1505.02023. 
#'
#' @section Main functions:
#' The main functions are \itemize{
#' \item \code{\link{clt_test}},
#' \item \code{\link{gaussian_bootstrap_test}}, 
#' \item \code{\link{empirical_bootstrap_test}},
#' \item \code{\link{HS_gaussian_bootstrap_test}},
#' \item \code{\link{HS_empirical_bootstrap_test}}
#' }
#' 
#'
#' @docType package
#' @name covsep
NULL










#' Generate a sample from a Matrix Gaussian distribution
#' 
#' @param N sample size
#' @param C1 row covariance
#' @param C2 column covariance
#' @param M mean matrix
#' @return A \code{N x dim(C1)[1] x dim(C2)[1]} array containing the
#' generated data
#' @examples 
#' Data = rmtnorm(30, C1, C2)
#' @export
# TODO ' optimize!

rmtnorm = function(N, C1, C2, M=matrix(0, nrow(C1), nrow(C2)) ){
    Data<-array(NA, c(N,dim(C1)[1],dim(C2)[1]))

    ## compute square root (using Cholesky) # use 'eigen' if there is an error
    SC1 <- tryCatch({t(chol(C1))},
        error = function(err){
            eigen.tmp <- eigen(C1, symmetric=TRUE)
            sqrt.diag  <- diag(sqrt(pmax(eigen.tmp$values, 0)), dim(C1)[1]) # correct the eigenvalues if they are negative, and take square root
            return( (eigen.tmp$vectors %*% sqrt.diag %*% t(eigen.tmp$vectors)) )
        },
        finally = {
        })

    SC2 <- tryCatch({t(chol(C2))},
        error = function(err){
            eigen.tmp <- eigen(C2, symmetric=TRUE)
            sqrt.diag  <- diag(sqrt(pmax(eigen.tmp$values, 0)), dim(C2)[1]) # correct the eigenvalues if they are negative, and take square root
            return( (eigen.tmp$vectors %*% sqrt.diag %*% t(eigen.tmp$vectors)) )
        },
        finally = {
        })



    # old
    #SC1<-svd1$u%*%diag(sqrt(svd1$d))%*%t(svd1$v)
    #SC2<-svd2$u%*%diag(sqrt(svd2$d))%*%t(svd2$v)
    for (nn in 1:N){ 
        ## this is faster
        Res<-array(stats::rnorm(prod(dim(Data)[2:3]),0,1),dim(Data)[2:3])
        #      for (i in 1:dim(Res)[1]){
        #        for (j in 1:dim(Res)[2]){ 
        #        Res[i,j]<-rnorm(1,0,1)
        #      }
        #      }
        Data[nn,,]<-SC1%*%Res%*%t(SC2)+M
    }
    return(Data)
}


#' Generate surface data 
#' 
#' Generate samples of surface data 
#' @inheritParams rmtnorm
#' @param gamma parameter to specify how much the covariance is separable. 
#' @param distribution distribution of the data
#' 
#' @section Details:
#' \code{gamma} can take values between 0 and 1; \code{gamma=0} corresponds to
#' a separable covariance, \code{gamma=1} corresponds to a non-separable
#' covariance (described in the paper \url{http://arxiv.org/abs/1505.02023}).
#' Values of \code{gamma} between 0 and 1 corresponds to an interpolation between
#' these two covariances
#'
#' \code{ distribution } can take the values 'gaussian' or 'student'
#' @return A \code{N x dim(C1)[1] x dim(C2)[1]} array containing the
#' generated data
#' 
#' @examples 
#' Data = generate_surface_data(30, C1, C2, gamma=0)
#'
#' @references \cite{John A. D. Aston, Davide Pigoli, Shahin Tavakoli, "Tests
#' for separability in nonparametric covariance operators of random surfaces",
#' 2015, under revision}, \url{http://arxiv.org/abs/1505.02023}
#' @export


generate_surface_data  <- function(N, C1, C2, gamma, distribution='gaussian'){

    # Vectorization map to project 4d cov to 2d
    I<-matrix(0,dim(C1)[1]*dim(C2)[1],2)
    for (i in 1:dim(I)[1]){
        I[i,]<-c(ceiling(i/dim(C2)[1]),i-(ceiling(i/dim(C2)[1])-1)*dim(C2)[1])
    }

    # Define the covariance operator for the simulation, based on gam
    CL<-matrix(0,dim(C1)[1]*dim(C2)[1],dim(C1)[1]*dim(C2)[1])
    x<-seq(0,1,len=dim(C1)[1])
    y<-seq(0,1,len=dim(C2)[1])
    for (p in 1:(dim(C1)[1]*dim(C2)[1])){
        for (q in 1:(dim(C1)[1]*dim(C2)[1])){
            CL[p,q]<-gamma*(1/((I[p,2]-I[q,2])^2+1)*exp(-(I[p,1]-I[q,1])^2/((I[p,2]-I[q,2])^2+1)))+(1-gamma)*C1[I[p,1],I[q,1]]*C2[I[p,2],I[q,2]]
        }
    }

    # data simulation
    if( distribution == 'gaussian' ){ 
        X<-mvtnorm::rmvnorm(N,mean=rep(0,dim(CL)[1]),sigma=CL)
    }else if (distribution == 'student'){   
    # correlation operator for the multivariate t simulations
        CLC<-stats::cov2cor(CL)
        X<-mvtnorm::rmvt(N,sigma=CLC,df=6)
    } else {
        stop(paste0('unknown distribution ', distribution))
    }

    # inverse of the vectorization:

    Y<-array(0,c(N,dim(C1)[1],dim(C2)[1]))
    for (i in 1:dim(X)[2]){
        Y[,ceiling(i/dim(Y)[3]),i-(ceiling(i/dim(Y)[3])-1)*dim(Y)[3]]<-X[,i]
    }
    return(Y)
} 

############## Description of the Datasets

#' A covariance matrix
#' 
#' Marginal covariance matrix C1 used for simulations in the paper \url{http://arxiv.org/abs/1505.02023}
#' 
#' @section Details:
#' This is a 32x32 real-valued covariance matrix.  
#' 
#' @examples
#' data(C1)
#' str(C1) 
#' 
#' @references \cite{John A. D. Aston, Davide Pigoli, Shahin Tavakoli, "Tests
#' for separability in nonparametric covariance operators of random surfaces",
#' 2015, under revision}, \url{http://arxiv.org/abs/1505.02023}

"C1"



#' A covariance matrix
#' 
#' Marginal covariance matrix C2 used for simulations in paper
#' \url{http://arxiv.org/abs/1505.02023}
#' 
#' @section Details: This is a 7x7 real-valued covariance matrix. 
#' @examples 
#' data(C2)
#' str(C2)
#' 
#' @references \cite{John A. D. Aston, Davide Pigoli, Shahin Tavakoli, "Tests
#' for separability in nonparametric covariance operators of random surfaces",
#' 2015, under revision}, \url{http://arxiv.org/abs/1505.02023}

"C2"



#' A data set of surfaces
#' 
#' Dataset of 50 surfaces simulated from a Gaussian process with a separable covariance structure. \code{SurfacesData[i,,]} corresponds to the i-th surface observed on a 32x7 uniform grid.
#' @section Details:
#'  This is a \code{50 x 32 x 7} array.
#' @examples
#' data(SurfacesData)
#' image(SurfacesData[1,,]) # color image of the first surface in the dataset
"SurfacesData"




#' renormalize a matrix normal random matrix to have iid entries
#' @param X a matrix normal random matrix with mean zero
#' @inheritParams rmtnorm
#' @param type the type of renormalization to do.  Possible options are 'no',
#' 'diag' or 'full' (see details section).
#' @section Details:
#' \code{ type } can take the values \describe{
#'     \item{'diag'}{each entry of \code{X} is renormalized by its marginal standard deviation}
#'     \item{'full'}{\code{X} is renormalized by its root inverse covariance}
#' }
#' @return A matrix with renormalized entries
#' @examples 
#' Data  <-  rmtnorm(30, C1, C2)
#' ans  <-  renormalize_mtnorm(Data[1,,], C1, C2)
#' @export
### needs to be optimized

renormalize_mtnorm  =  function(X, C1, C2, type='full')
{
    ## check compatibility of Data
#    C1  <- as.matrix(C1)
#    C2  <- as.matrix(C2)
#    if( !is.matrix(X) ){ # X should be a vector of length nrow(C2)
#        if( 1 != ncol(C1)) stop("dimension missmatch between X and C1")
#        if( length(X) != ncol(C2) ) stop("dimension missmatch between X and C2")
#    } else
    if( !(is.matrix(X) && is.matrix(C1) && is.matrix(C2)) ) stop("renormalize_mtnorm's arguments must be passed as matrices")
    if( nrow(X) != ncol(C1)) stop("dimension missmatch between X and C1")
    if( ncol(X) != ncol(C2)) stop("dimension missmatch between X and C2")

    if( type == 'full' ){
    svd1 = svd(C1)
    svd2 = svd(C2)
    ans =  svd1$u %*% diag((1/sqrt(svd1$d)), nrow(C1)) %*% t(svd1$v) %*% X %*%  svd2$u %*% diag( (1/sqrt(svd2$d)), nrow(C2) ) %*% t(svd2$v) 
    } else if ( type == 'diag' ){
    ans =  diag(1/sqrt(diag(C1)), nrow(C1))  %*% X %*%   diag(1/sqrt(diag(C2)), nrow(C2)) 
    } else {
        stop( paste0( 'unknown type ', type) )
    } 
    return(ans)
}

#' estimates marginal covariances (e.g. row and column covariances) of bi-dimensional sample
#' @param Data a (non-empty) \code{N x d1 x d2} array of data values. The first
#' direction indices the \eqn{N} observations, each consisting of a \code{d1 x d2}
#' discretization of the surface, so that \code{Data[i,,]} corresponds to the
#' i-th observed surface.
#' @return A list containing the row covariance (\code{C1}) and column covariance (\code{C2})
#' @examples
#' Data  <-  rmtnorm(30, C1, C2)
#' marginal.cov  <- marginal_covariances(Data)
#' @export

marginal_covariances  = function(Data){

    ## size of the data
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]

    Data=sweep(Data, c(2,3), apply(Data, c(2,3), mean)) # remove the mean from the data
    
    ## compute marginal covariance estimate
    tmp.mat=matrix(c(aperm(Data, c(1,3,2))), ncol=d1)
    C1=stats::cov(tmp.mat)*(N*d2-1)/N
    rm(tmp.mat)

    ## compute marginal covariance estimate
    tmp.mat=matrix(c(Data), ncol=d2)
    C2=stats::cov(tmp.mat)*(N*d1-1)/N
    rm(tmp.mat)

    C1<-C1/sqrt(sum(diag(C1)))
    C2<-C2/sqrt(sum(diag(C2)))
    return(list(C1=C1, C2=C2))

}

#' Compute the projection of the rescaled difference between the sample covariance
#' and its separable approximation onto the separable eigenfunctions
#' @inheritParams marginal_covariances
#' @param l1 number of eigenfunctions to be used in the first (row) dimension for the projection
#' @param l2 number of eigenfunctions to be used in the second (column) dimension for the projection
#' @param with.asymptotic.variances logical variable; if TRUE, the function outputs the estimate asymptotic variances of the projected differences
#' @return A list with
#' \describe{
#'  \item{T.N}{The projected differences}
#' \item{sigma.left}{The row covariances of \code{T.N}}
#' \item{sigma.right}{The column covariances of \code{T.N}}
#' }
#' @section Details:
#' The function computes the projection of the rescaled difference between the sample covariance
#' and its separable approximation onto the separable eigenfunctions \code{u_i
#' x v_j : i = 1, \ldots, l1; j = 1, \ldots, l2}.
#'
#' @examples 
#' Data  <-  rmtnorm(30, C1, C2)
#' ans <- projected_differences(Data, l1=1, l2=2)
#' @export

projected_differences = function(Data, l1=1, l2=1, with.asymptotic.variances=TRUE)
{
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]

    Data=sweep(Data, c(2,3), apply(Data, c(2,3), mean)) # remove the mean from the data

    ## compute marginal covariances
    marginal.covariances = marginal_covariances(Data)

    ## and their eigenstructure
    svd1 <- svd(marginal.covariances$C1)
    svd2 <- svd(marginal.covariances$C2)
    lambda=svd1$d
    gamma=svd2$d


    # computation of the matrix of the centering terms for the pair of eigenfunctions (l1,l2)
    shift.stat<-matrix(NA,l1,l2) 

    for (rowi in 1:l1){ 
        for (coli in 1:l2){
            # estimation of the projection of the covariance operator on the
            # eigendirections rowi and coli
            tmp <-numeric(1)
            for (i in 1:N){
                tmp <-tmp +(as.double(t(svd1$u[,rowi])%*%(Data[i,,])%*%svd2$u[,coli]))^2
            }
            tmp <-tmp/N
            shift.stat[rowi,coli]<-tmp -lambda[rowi]*gamma[coli]
        }
    }
    shift.stat = sqrt(N)*shift.stat # rescale it

    ## values by default
    sigma.left=NULL
    sigma.right=NULL

    if(with.asymptotic.variances){
        sigma.left  <- sqrt(2)*outer(lambda[1:l1], lambda[1:l1]) * (diag(l1)*(sum(lambda)^2)+ matrix(sum(lambda^2), l1, l1) - sum(lambda)*outer(lambda[1:l1], lambda[1:l1], '+') )/(sum(lambda)*sum(gamma))
        #
        sigma.right  <- sqrt(2)*outer(gamma[1:l2], gamma[1:l2]) * (diag(l2)*(sum(gamma)^2)+ matrix(sum(gamma^2), l2, l2) - sum(gamma)*outer(gamma[1:l2], gamma[1:l2], '+') )/(sum(lambda)*sum(gamma))
    }

    ans=list(T.N = shift.stat, sigma.left=sigma.left, sigma.right=sigma.right)
    return(ans)
}


#' Test for separability of covariance operators for Gaussian process.
#'
#' This function performs the asymptotic test  for the separability of the covariance operator for a random surface generated from a Gaussian process (described in the paper \url{http://arxiv.org/abs/1505.02023}).
#'
#' @inheritParams marginal_covariances
#' @param L1 an integer or vector of integers in \eqn{1:p} indicating the
#' eigenfunctions in the first direction to be used for the test. 
#' @param L2 an integer or vector of integers in \eqn{1:q} indicating the
#' eigenfunctions in the second direction to be used for the test. 
#' 
#' @section Details:
#' 
#' If L1 and L2 are vectors, they need to be of the same length.
#' 
#' The function tests for separability using the projection of the covariance
#' operator in the separable eigenfunctions \code{u_i tensor v_j : i = 1, \ldots, l1;
#' j = 1, \ldots, l2}, for each pair \code{(l1,l2) = (L1[k], L2[k])}, for \code{k = 1:length(L1)}.
#' 
#' The test works by using asymptotics, and is
#' only valid if the data is assumed to be Gaussian.
#' 
#' The surface data needs to be measured or resampled on a common regular grid
#' or on common basis functions.
#'
#' @return The p-value of the test for each pair \code{(l1,l2) = (L1[k], L2[k])}, for \code{k = 1:length(L1)}.
#' @seealso \code{\link{empirical_bootstrap_test}}, \code{\link{gaussian_bootstrap_test}}
#' 
#' @examples 
#' data(SurfacesData)
#' clt_test(SurfacesData, L1=c(1,2), L2=c(1,4))
#'
#' @references \cite{John A. D. Aston, Davide Pigoli, Shahin Tavakoli, "Tests
#' for separability in nonparametric covariance operators of random surfaces",
#' 2015, under revision}, \url{http://arxiv.org/abs/1505.02023}
#' @export

clt_test  <- function(Data, L1, L2)
{
    #browser()
    ## compute projected differences
    proj.diff <- projected_differences(Data, max(L1), max(L2), with.asymptotic.variances=TRUE)
    ## compute chi-2 statistic and its p-value
    pvals = array(NA, length(L1))
    for(k in 1:length(L1)){
        l1  <-  L1[k]
        l2  <- L2[k]
    pvals[k]=stats::pchisq( sum(renormalize_mtnorm(   proj.diff$T.N[1:l1, 1:l2, drop=FALSE] ,
                                               proj.diff$sigma.left[1:l1, 1:l1, drop=FALSE],
                                               proj.diff$sigma.right[1:l2, 1:l2, drop=FALSE], type='full')^2),
                       df=l1*l2, lower.tail=F)
    }
    return(pvals)
}


#'  Projection-based Gaussian (parametric) bootstrap test for separability of covariance structure
#' 
#' This function performs the test for the separability of covariance structure
#' of a random surface generated from a Gaussian process, based on the
#' parametric bootstrap procedure described in the paper
#' \url{http://arxiv.org/abs/1505.02023}
#' @inheritParams clt_test
#' @param studentize parameter to specify which type of studentization is performed. Possible options are 'no', 'diag' or 'full' (see details section).
#' @param B number of bootstrap replicates to be used.
#' @param verbose logical parameter for printing progress
#' 
#' @section Details:
#' 
#' This function performs the test of separability 
#' of the covariance structure for a random surface (introduced in the paper
#' \url{http://arxiv.org/abs/1505.02023}), when generated from a Gaussian
#' process. The sample surfaces need to be measured on a common regular grid. The test
#' consider a subspace formed by the tensor product of eigenfunctions of the separable
#' covariances. It is possible to specify the number of eigenfunctions to be considered
#' in each direction. 
#' 
#' If L1 and L2 are vectors, they need to be of the same length.
#' 
#' The function tests for separability using the projection of the covariance
#' operator in the separable eigenfunctions \code{u_i x v_j : i = 1, \ldots, l1;
#' j = 1, \ldots, l2}, for each pair (l1,l2) = (L1[k], L2[k]), for k = 1:length(L1).
#' 
#' \code{studentize} can take the values \describe{
#'     \item{'no'}{no studentization is performed}
#'     \item{'diag'}{each projection coordinate is renormalized by an estimate of its standard deviation}
#'     \item{'full'}{the projection coordinates are renormalized by an estimate of their joint covariance}
#' }
#'
#' @return The p-value of the test for each pair \code{(l1,l2) = (L1[k], L2[k])}, for \code{k = 1:length(L1)}.
#' 
#' @seealso \code{\link{empirical_bootstrap_test}}, \code{\link{clt_test}}
#' 
#' @examples
#' data(SurfacesData)
#' \dontrun{gaussian_bootstrap_test(SurfacesData,L1=1,L2=1,B=100,
#' studentize='full')} #' pvalue of the test
#'
#' @references \cite{John A. D. Aston, Davide Pigoli, Shahin Tavakoli, "Tests
#' for separability in nonparametric covariance operators of random surfaces",
#' 2015, under revision}, \url{http://arxiv.org/abs/1505.02023}
#' @export


gaussian_bootstrap_test  <- function(Data, L1, L2, studentize='full', B=100, verbose=FALSE)
{
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]

    ## compute value of test statistic and save C1, C2
    marginal.cov <- marginal_covariances(Data)
    proj.diff <- projected_differences(Data, max(L1), max(L2), with.asymptotic.variances=TRUE)
    stat = array(NA, length(L1))
    for (k in 1:length(L1)){ 
        l1 = L1[k]
        l2 = L2[k]
        if( studentize == 'no'){
            stat[k] = sum(proj.diff$T.N[1:l1, 1:l2]^2)
        } else if (studentize == 'diag' || studentize == 'full'){
            stat[k] = sum(renormalize_mtnorm(  proj.diff$T.N[1:l1, 1:l2, drop=FALSE] ,
                                             proj.diff$sigma.left[1:l1, 1:l1, drop=FALSE],
                                             proj.diff$sigma.right[1:l2, 1:l2, drop=FALSE], 
                                             type=studentize)^2)
        } else {
            stop( paste0( 'unknown value for studentize: ', studentize) )
        }
    }
    ## Bootstrap runs
    boot.stat = array(NA, c(B, length(L1)) )
    for (b in 1:B){
        ## print bootstrap  progress
        if(verbose){
            if((b %% floor(1+B/200))==0){ cat("GB",round(100*b/B,1), "%--", sep="") }
        }
        ## generate sample (TODO : replace with more efficient code?)
        Data.boot <- rmtnorm(N, marginal.cov$C1,  marginal.cov$C2) 
        ## compute projected differences
        proj.diff.boot <- projected_differences(Data.boot, max(L1), max(L2), with.asymptotic.variances=TRUE)
        ## compute renormalized test statistic
        for ( k in 1:length(L1) ){
            l1 = L1[k]
            l2 = L2[k]
            if( studentize == 'no'){
                boot.stat[b,k] = sum(proj.diff.boot$T.N[1:l1, 1:l2]^2)
            } else if (studentize == 'diag' || studentize == 'full'){
                boot.stat[b,k] = sum(renormalize_mtnorm(
                                                        proj.diff.boot$T.N[1:l1,
                                                                           1:l2,
                                                                           drop=FALSE]
                                                        ,
                                                        proj.diff.boot$sigma.left[1:l1,
                                                                                  1:l1,
                                                                                  drop=FALSE],
                                                        proj.diff.boot$sigma.right[1:l2,
                                                                                   1:l2,
                                                                                   drop=FALSE],
                                                        type=studentize)^2)
            } else {
                stop( paste0( 'unknown value for studentize: ', studentize) )
            }
        }
        rm(proj.diff.boot)
    }
    ## compute and  return pvalue
    pvalues<-rep(NA,dim(boot.stat)[2])
    for (k in 1:length(L1)){ 
        pvalues[k]<-mean(boot.stat[,k]>stat[k])
    }
    return(pvalues)
}

#'   Projection-based empirical bootstrap test for separability of covariance structure
#' 
#' This function performs the test for the separability of covariance structure
#' of a random surface based on the empirical bootstrap procedure described in
#' the paper \url{http://arxiv.org/abs/1505.02023}.
#'
#' @inheritParams gaussian_bootstrap_test
#'
#' @section Details:
#' 
#' This function performs the test of separability 
#' of the covariance structure for a random surface (introduced in the paper
#' \url{http://arxiv.org/abs/1505.02023}), when generated from a Gaussian
#' process. The sample surfaces need to be measured on a common regular grid. The test
#' consider a subspace formed by the tensor product of eigenfunctions of the separable
#' covariances. It is possible to specify the number of eigenfunctions to be considered
#' in each direction. 
#' 
#' If L1 and L2 are vectors, they need to be of the same length.
#' 
#' The function tests for separability using the projection of the covariance
#' operator in the separable eigenfunctions \code{u_i x v_j : i = 1, \ldots, l1;
#' j = 1, \ldots, l2}, for each pair (l1,l2) = (L1[k], L2[k]), for k = 1:length(L1).
#' 
#' \code{studentize} can take the values \describe{
#'     \item{'no'}{no studentization is performed}
#'     \item{'diag'}{each projection coordinate is renormalized by an estimate of its standard deviation}
#'     \item{'full'}{the projection coordinates are renormalized by an estimate of their joint covariance}
#' }
#'
#' @return The p-value of the test for each pair \code{(l1,l2) = (L1[k], L2[k])}, for \code{k = 1:length(L1)}.
#' 
#' @seealso \code{\link{gaussian_bootstrap_test}}, \code{\link{clt_test}}
#' 
#' @examples 
#' data(SurfacesData)
#' \dontrun{empirical_bootstrap_test(SurfacesData,L1=1,L2=1, B=100, studentize='full')}
#' @references \cite{John A. D. Aston, Davide Pigoli, Shahin Tavakoli, "Tests
#' for separability in nonparametric covariance operators of random surfaces",
#' 2015, under revision}, \url{http://arxiv.org/abs/1505.02023}
#' @export

empirical_bootstrap_test  <- function(Data, L1, L2, studentize, B=100, verbose=FALSE)
{
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]

    ## compute value of test statistic and save C1, C2
    marginal.cov <- marginal_covariances(Data)
    proj.diff <- projected_differences(Data, max(L1), max(L2), with.asymptotic.variances=TRUE)
    stat = array(NA, length(L1))
    for (k in 1:length(L1)){ 
        l1 = L1[k]
        l2 = L2[k]
        if( studentize == 'no'){
            stat[k] = sum(proj.diff$T.N[1:l1, 1:l2]^2)
        } else if (studentize == 'diag' || studentize == 'full'){
            stat[k] = sum(renormalize_mtnorm(   proj.diff$T.N[1:l1, 1:l2, drop=FALSE],
                                             proj.diff$sigma.left[1:l1, 1:l1, drop=FALSE],
                                             proj.diff$sigma.right[1:l2, 1:l2, drop=FALSE], type=studentize)^2)
        } else {
            stop( paste0( 'unknown value for studentize: ', studentize) )
        }
    }
    ## computation of the bootstrap distribution 
    if(verbose) cat("*** Empirical Bootstrap ***\n")
    boot.stat<-matrix(NA,B,length(L1))
    for (b in 1:B){
        ## print bootstrap  progress
        if(verbose){
        if((b %% floor(1+B/200))==0){ cat("EB",round(100*b/B,1), "%--", sep="") }
        }
        # resample from the empirical density
        Data.boot<-Data[sample.int(N,N,replace=TRUE),,]
        proj.diff.boot = projected_differences(Data.boot, max(L1), max(L2), with.asymptotic.variances=TRUE)
        for ( k in 1:length(L1) ){
            l1 = L1[k]
            l2 = L2[k]
            if( studentize == 'no'){
                boot.stat[b,k]  <-  sum((proj.diff.boot$T.N[1:l1, 1:l2] -
                                         proj.diff$T.N[1:l1, 1:l2])^2 )
            } else if (studentize == 'diag' || studentize == 'full'){
                boot.stat[b,k] = sum(renormalize_mtnorm(   proj.diff.boot$T.N[1:l1, 1:l2, drop=FALSE] - proj.diff$T.N[1:l1, 1:l2, drop=FALSE]  ,
                                                        proj.diff.boot$sigma.left[1:l1, 1:l1, drop=FALSE],
                                                        proj.diff.boot$sigma.right[1:l2, 1:l2, drop=FALSE], type=studentize)^2)
            } else {
                stop( paste0( 'unknown value for studentize: ', studentize) )
            }
        }
        rm(proj.diff.boot)
    }
    ## compute and  return pvalue
    pvalues<-rep(NA,dim(boot.stat)[2])
    for (k in 1:length(L1)){ 
        pvalues[k]<-mean(boot.stat[,k]>stat[k])
    }
    return(pvalues)
}


#' compute the difference between the full sample covariance and its separable
#' approximation
#' @inheritParams marginal_covariances
#' @return A \code{d1 x d2 x d1 x d2} array, where \code{d1 = nrow(Data)} and
#' \code{d2 = ncol(Data)}.
#' @section Details:
#' This is an internal function.
## to be optimized
difference_fullcov <- function(Data){
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]
    ## remove the mean from the data ### no need, should be passed only after removing the mean
    Data=sweep(Data, c(2,3), apply(Data, c(2,3), mean))

    marginal.cov  <- marginal_covariances(Data)

    ## compute Difference 'D'
    tmp.mat=matrix(Data, nrow=N)
    Difference.D = array(stats::cov(tmp.mat)*(N-1)/N, c(d1,d2,d1,d2)) - aperm(outer(marginal.cov$C1,marginal.cov$C2), c(1,3,2,4))
    return((Difference.D))
}


#' Gaussian (parametric) bootstrap test
#' for separability of covariance structure using Hilbert--Schmidt distance
#' @inheritParams gaussian_bootstrap_test
#'
#' @section Details:
#' 
#' This function performs the test of separability 
#' of the covariance structure for a random surface (introduced in the paper
#' \url{http://arxiv.org/abs/1505.02023}), when generated from a Gaussian
#' process. The sample surfaces need to be measured on a common regular grid. The test
#' considers the Hilbert--Schmidt distance between the sample covariance and its separable approximation.
#' 
#'
#' @return The p-value of the test.
#' 

#' @examples 
#' data(SurfacesData)
#' \dontrun{HS_gaussian_bootstrap_test(SurfacesData, B = 100)}
#' @export

HS_gaussian_bootstrap_test  <- function(Data, B=100, verbose=FALSE)
{
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]

    ## compute value of test statistic and save C1, C2
    marginal.cov <- marginal_covariances(Data)
    stat = sum( difference_fullcov(Data)^2 )
#
    boot.stat = array(NA, B)
    ## Bootstrap runs
    for (b in 1:B){
        ## print bootstrap  progress
        if(verbose){
        if((b %% floor(1+B/200))==0){ cat("GB",round(100*b/B,1), "%--", sep="") }
        }
        ## generate sample (TODO : replace with more efficient code?)
        Data.boot <- rmtnorm(N, marginal.cov$C1,  marginal.cov$C2) 
        ## compute test statistic
        boot.stat[b]  <- sum( difference_fullcov(Data.boot)^2 )
    }
#
    ## compute and  return pvalue
    pvalue<-mean(boot.stat>stat)
    return(pvalue)
}

#' Empirical bootstrap test for
#' separability of covariance structure using Hilbert--Schmidt distance
#' @inheritParams gaussian_bootstrap_test
#'
#' @section Details:
#' 
#' This function performs the test of separability 
#' of the covariance structure for a random surface (introduced in the paper
#' \url{http://arxiv.org/abs/1505.02023}), when generated from a Gaussian
#' process. The sample surfaces need to be measured on a common regular grid. The test
#' considers the Hilbert--Schmidt distance between the sample covariance and its separable approximation.
#'
#' @return The p-value of the test.
#' 
#' @examples 
#' data(SurfacesData)
#' \dontrun{HS_empirical_bootstrap_test(SurfacesData, B = 100)}
#' @export

HS_empirical_bootstrap_test  <- function(Data, B=100, verbose=FALSE)
{
    N<-dim(Data)[1]
    d1 <- dim(Data)[2]
    d2 <- dim(Data)[3]
    ## compute value of test statistic and save C1, C2
    marginal.cov <- marginal_covariances(Data)
    stat = sum( difference_fullcov(Data)^2 )
#
    boot.stat = array(NA, B)
    ## Bootstrap runs
    for (b in 1:B){
        ## print bootstrap  progress
        if(verbose){
            if((b %% floor(1+B/200))==0){ cat("GB",round(100*b/B,1), "%--", sep="") }
        }
        # resample from the empirical density
        Data.boot<-Data[sample.int(N,N,replace=TRUE),,]
        ## compute test statistic
        boot.stat[b]  <- sum( (difference_fullcov(Data.boot) - difference_fullcov(Data))^2 )
    }
#
    ## compute and  return pvalue
    pvalue<-mean( boot.stat > stat )
    return(pvalue)
}

