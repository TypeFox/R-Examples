#' Split connectivity matrix into subpopulations
#'
#' This function tries to optimally split a given subpopulation into
#' two smaller subpopulations.
#'
#' @param indices vector of indices of sites in a subpopulation
#' @param conn.mat a square connectivity matrix.  This matrix has
#' typically been normalized and made symmetric prior to using this
#' function.
#' @param beta controls degree of splitting of connectivity matrix,
#' with larger values generating more subpopulations.
#' @param tries how many times to restart the optimization algorithm. Defaults to 5.
#' @param threshold controls when to stop each "try".  Defaults to 1e-10.
#' @param alpha controls rate of convergence to solution
#' @param maxit Maximum number of iterations to perform per "try".
#'
#' @return List with one or two elements, each containing a vector of
#' indices of sites in a subpopulations
#'
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson,
#' P. R. 2012. Identification of subpopulations from connectivity
#' matrices. Ecography, 35: 1004-1016.
#' 
#' @seealso See also \code{\link{optimalSplitConnMat}},
#' \code{\link{recSplitConnMat}},
#' \code{\link{subpopsVectorToList}}
#' 
#' @author
#' David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
splitConnMat <- function(indices,conn.mat,beta,tries=5,
                    threshold=1e-10,alpha=0.1,maxit=500) {
    # makes a submatrix of the total connectivity matrix only
    # involving the sites in the index list
    ppp = conn.mat[indices, indices, drop=F]

    # al appears to be just 1/beta
    al = 1/beta
    
    n = dim(ppp)[2];
    s = matrix(1,nrow=n)
    
    eNoSplit = -t(s) %*% ppp %*% s + al * sum(s)^2;

    best = 1e10
    
    for (kk in 1:tries) {
        s = matrix( runif(n,min=-1,max=1), ncol = 1 )
        ds = s + 1
        sTot = sum(s)
        sOld = sign(s)

        for (it in 1:maxit) {
            if (sqrt(sum((s-ds)^2)) <= threshold)
                break

            # the following three lines implements EQ. 8 in the paper
            v = ppp %*% s - al * sum(s)            
            ds = sqrt(abs(v)) * sign(v)
            s = alpha * s + (1-alpha) * ds
        }

        if (it == maxit) warning("Reached max iterations")

        s = sign(s)
        e = -t(s) %*% ppp %*% s + al * sum(s)^2;
        # calclulates the value of the cost function

        if (e < best) {
            sBest = s
            best = e
        }
    }

    part = indices[ sBest == -1 ]
    notPart = setdiff( indices, part )

    if ( (best < eNoSplit) & (length(part)>0) & (length(notPart)>0) ) {
        return(list(part,notPart))
    }

    return(list(indices))
}

#' Recursively subdivides a set of subpoplations
#' 
#' This funtion recursively splits each subpopulation of a list of
#' subpopulations until none of the subpopulations can be split
#' further to improve the minimization.
#'
#' @param subpops.lst A list whose elements are vectors of indices for each subpopulation.  See \code{\link{subpopsVectorToList}}.
#' @param conn.mat A square connectivity matrix.  This matrix has
#' typically been normalized and made symmetric prior to using this
#' function.
#' @param beta Controls degree of splitting of connectivity matrix,
#' with larger values generating more subpopulations.
#' @param \dots further arguments to be passed to \code{\link{splitConnMat}}
#' 
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson,
#' P. R. 2012. Identification of subpopulations from connectivity
#' matrices. Ecography, 35: 1004-1016.
#' 
#' @seealso See also \code{\link{optimalSplitConnMat}},
#' \code{\link{splitConnMat}},
#' \code{\link{subpopsVectorToList}}
#' 
#' @author
#' David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
recSplitConnMat <- function(subpops.lst, conn.mat, beta, ...) {
    old = 0
    while (length(subpops.lst) > old) {
        old = length(subpops.lst)
        subpops.lst = unlist( sapply( subpops.lst, splitConnMat,
            conn.mat=conn.mat, beta=beta, ..., simplify=F),
            recursive=F )
    }
    return(subpops.lst)
}

#' Merge subpopulations
#'
#' This function tries to merge random subopoulations, checking if the
#' result is a better soluton to the minimization problem.
#'
#' @param subpops.lst A list whose elements are vectors of indices for each subpopulation.  See \code{\link{subpopsVectorToList}}.
#' @param conn.mat A square connectivity matrix.  This matrix has
#' typically been normalized and made symmetric prior to using this
#' function.
#' @param beta Controls degree of splitting of connectivity matrix,
#' with larger values generating more subpopulations.
#'
#' @return List of the same format as subpops.lst, but with
#' potentially fewer subpopulations.
#'
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson,
#' P. R. 2012. Identification of subpopulations from connectivity
#' matrices. Ecography, 35: 1004-1016.
#' 
#' @seealso See also \code{\link{optimalSplitConnMat}},
#' 
#' @author
#' David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
mergeSubpops <- function ( subpops.lst,  conn.mat, beta ) {
    nIt = length(subpops.lst)^2
    al = 1/beta
    
    for (it in 1:nIt) {
        if (length(subpops.lst) < 2) break # Must have at least 2 to merge
        
        ii = sort(sample( 1:length(subpops.lst), 2 ))
        li = c( subpops.lst[[ ii[1] ]], subpops.lst[[ ii[2] ]] )
        pTest = conn.mat[li,li]
        
        s = matrix(1,nrow=length(li))
        eTogether = -t(s) %*% pTest %*% s + al * sum(s)^2
        s[1:length(subpops.lst[[ ii[1] ]])] = -1
        eSplit = -t(s) %*% pTest %*% s + al * sum(s)^2

        if (eTogether<eSplit) {
            subpops.lst[[ ii[1] ]] = li
            subpops.lst = subpops.lst[ -ii[2] ]
        }
    }

    return(subpops.lst)
}

#' Reduced connectivity matrix according to a set of subpopulations
#' 
#' Reduces a connectivity matrix based on a set of subpopulations.  If there are
#' N subpopulations, then the reduced matrix will have dimensions NxN.  The 
#' reduced matrix will be ordered according to the order of subpopulations in 
#' \code{subpops.lst}.
#' 
#' @param subpops.lst A list whose elements are vectors of indices for each 
#'   subpopulation.  If a vector of integers is given, then 
#'   \code{\link{subpopsVectorToList}} is applied to convert it to a list of 
#'   subpopulations.
#' @param conn.mat A square connectivity matrix.
#'   
#' @return A reduced connectivity matrix.  The sum of all elements of this
#'   reduced connectivity matrix will be equal to the sum of all elements of the
#'   original connectivity matrix.
#'   
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson, P. R. 2012. 
#'   Identification of subpopulations from connectivity matrices. Ecography, 35:
#'   1004-1016.
#'   
#' @seealso See also \code{\link{qualitySubpops}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @example tests/test.optimalSplitConnMat.R
#' @export
reducedConnMat <- function( subpops.lst, conn.mat ) {
  if (!is.list(subpops.lst))
    subpops.lst = subpopsVectorToList(subpops.lst)
  
  pii = matrix( 0.0, nrow=dim(conn.mat)[2], ncol=length(subpops.lst) )
  
  for (kk in 1:length(subpops.lst)) {
    pii[ subpops.lst[[kk]], kk ] = 1
  }
  
  # Note I use p instead of t(p), as was in Jacobi code.
  # This is because my matrices are oriented columns to rows
  #pt = t(pii) %*% p %*% pii %*% solve( t(pii) %*% pii )
  
  # Somewhat quicker method
  pt =  t(pii) %*% conn.mat %*% pii #%*% diag( 1 / sapply(subpops.lst,length) )
  # No need to normalize by number of sites in cluster of larval
  # origin because we are just going to normalize each colum
 
  return(pt)
}

#' Quality measure for subpopulation division
#' 
#' A measure of the leakage between subpopulations for a given division of the 
#' connectivity matrix into subpopulations.  This statistic is equal to 1 - 
#' mean(RLR) of the reduced connectivity matrix, where RLR=relative local 
#' retention (\code{\link{relativeLocalRetention}}), i.e., the fraction of 
#' settling individuals that originated at their site of settlement.
#' 
#' @inheritParams reducedConnMat
#'   
#' @return The quality statistic.
#'   
#'   A smaller value of the quality statistic indicates less leakage.
#'   
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson, P. R. 2012. 
#'   Identification of subpopulations from connectivity matrices. Ecography, 35:
#'   1004-1016.
#'   
#' @seealso See also \code{\link{optimalSplitConnMat}}, 
#'   \code{\link{subpopsVectorToList}}, \code{\link{relativeLocalRetention}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
#' @include retentionStats.R
qualitySubpops <- function( subpops.lst, conn.mat )  
    (1 - mean(relativeLocalRetention(reducedConnMat(subpops.lst,conn.mat))))

#' Compute vector of beta values
#' 
#' Helper function to compute a set of beta values using formula used in Jacobi
#' et al. (2012).
#' 
#' @param n numerator of formula from Jacobi et al. (2012).  Normally will be
#'   the number of columns in the connectivity matrix if one normalizes the
#'   columns (otherwise, it would typically be \code{N^2 / sum(conn.mat)}, where
#'   \code{N} is the number of columns of \code{conn.mat}.
#' @param steps number of beta values to return.  Defaults to 10.
#' @param cycles how many cycles of \code{2*pi} to do.
#' @param coeff coefficient in front of sine function
#' @param pwr exponent in denominator
#'   
#' @return vector of beta values
#'   
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson, P. R. 2012.
#'   Identification of subpopulations from connectivity matrices. Ecography, 35:
#'   1004-1016.
#'   
#' @seealso See also \code{\link{optimalSplitConnMat}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
betasVectorDefault <- function(n,steps=10,cycles=3/4,
                                           coeff=0.8,pwr=3.0)
  n/(1+coeff*sin(seq(0,cycles*2*pi,length.out=steps)))^pwr

#' Iteratively, optimally split a connectivity matrix
#' 
#' Algorithm for iteratively determining subpopulations of
#' highly-connected sites.  Uses an iterative method described in
#' Jacobi et al. (2012)
#'
#' @param conn.mat A square connectivity matrix.
#' @param normalize.cols A boolean indicating if columns of conn.mat
#' should be normalized by the sum of all elements in the column.
#' Defaults to TRUE.
#' @param make.symmetric A string indicating how to force conn.mat to
#' be symmetric.  "mean" (the default) will replace C_ij by (C_ij +
#' C_ji)/2.  "max" will replace C_ij by the maximum of C_ij and C_ji.
#' @param remove.diagonal A boolean indicating if the diagonal
#' elements of conn.mat should be removed before determining the
#' subpopulations.  Defaults to FALSE.
#' @param cycles Number of times to pass over values in betas.
#' @param betas Vector of beta values to try.  If not given, will
#' default to \code{\link{betasVectorDefault}(dim(conn.mat)[2],steps)}.
#' @param steps Number of beta values to produce using
#' betasVectorDefault.  Ignored if betas argument is explicitly
#' given.
#' @param \dots further arguments to be passed to \code{\link{splitConnMat}}
#'
#' @return A list with the following elements:
#' \item{betas}{Vector of all beta values tested}
#'
#' \item{num.subpops}{Vector of number of subpopulations found for
#' each value of beta}
#'
#' \item{qualities}{Vector of the quality statistic for each
#' subpopulation division}
#' 
#' \item{subpops}{A matrix with dimensions dim(conn.mat)[2] X
#' length(betas) indicating which subpopulation each site belongs to}
#' 
#' \item{best.splits}{A list indicating for each number of
#' subpopulations, which column of subpops contains the division with
#' the lowest quality statistic.  E.g.,
#' \code{best.splits[["4"]]$index} contains the column index of the
#' optimal division of the connectivity matrix into 4 subpopulations.}
#' 
#' @seealso See also \code{\link{splitConnMat}},
#' \code{\link{recSplitConnMat}}, \code{\link{mergeSubpops}},
#' \code{\link{qualitySubpops}}
#' @author
#' David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @example tests/test.optimalSplitConnMat.R
#' @export
#' 
#' @note In Jacobi et al. (2012) paper, the connectivity matrix is
#' oriented so that \eqn{C_ij} is dispersal from i to j, whereas in this R
#' package, the connectivity matrix is oriented so that \eqn{C_ij} is
#' dispersal from j to i.  This choice of orientation is arbitrary,
#' but one must always be consistent.  From j to i is more common in
#' population dynamics because it works well with matrix
#' multiplication (e.g., \code{settlers = conn.mat \%*\% eggs}).
#'
#' @references Jacobi, M. N., André, C., Döös, K., and Jonsson,
#'  P. R. 2012. Identification of subpopulations from connectivity
#'  matrices. Ecography, 35: 1004-1016.
#' 
optimalSplitConnMat <-
    function(conn.mat, normalize.cols=TRUE,
             make.symmetric="mean", remove.diagonal=FALSE, 
             cycles = 2, 
             betas=betasVectorDefault(
               ifelse(normalize.cols,dim(conn.mat)[2],
                      prod(dim(conn.mat))/sum(conn.mat)),steps),
             steps=10,
             ... ) {
    if (class(conn.mat) != "matrix")
        stop("Input conn.mat must be a matrix.")
    
    pp <- conn.mat

    # Make larval loss uniform over space
    if (normalize.cols)
        for (kk in 1:dim(pp)[2]) {
            ss = sum(pp[,kk])
            if (ss>0) pp[,kk] = pp[,kk] / ss
        }

    # Force symmetric if not already the case
    mymax = function(x) { ii = x < t(x); x[ii] = t(x)[ii]; return(x) }
    pp = switch(make.symmetric,
        mean = (pp + t(pp)) / 2.0,
        max = mymax(pp),
        stop("Bad max.symmetric string"))

    # Not sure if this is obligatory or optional
    # I believe it should be optional
    if (remove.diagonal)
        diag(pp) <- 0

    subpops = matrix(NA,nrow=dim(conn.mat)[2],ncol=cycles*length(betas))
    num.subpops = rep(NA,cycles*length(betas))
    qualities = num.subpops
    for (kk in 1:cycles) {
        print( paste( "Starting cycle",kk ) )

        # Initialize with one big cluster
        ta = list(1:dim(conn.mat)[2])

        # Loop over betas
        for (ll in 1:length(betas)) {
            beta = betas[ll]

            print( paste( "beta =", beta ) )
            
            taOld = list()
            while (!setequal(ta,taOld)) {
                taOld = ta
                ta = recSplitConnMat(ta, pp, beta, ...)
                #if (length(unlist(ta))!=dim(pp)[2]) stop("recSplitConnMat error")
                ta = mergeSubpops(ta,  pp, beta )
                #if (length(unlist(ta))!=dim(pp)[2]) stop("mergeSubpops error")
            }

            # Store results
            nn = (kk-1)*length(betas)+ll
            qualities[nn] = qualitySubpops(ta,conn.mat)
            num.subpops[nn] = length(ta)
            for (mm in 1:length(ta))
                subpops[ta[[mm]],nn] = mm
        }
    }

    # For each number of subpops, find the configuration with the
    # best quality measure
    imin <- function(x) { which( min(x) == x )[1] }
    best.splits = by( data.frame(quality=qualities,
        index=1:length(qualities)),
        INDICES=num.subpops,
        FUN=function(d) d[imin(d$quality),] )
    
    return(list(betas=rep(betas,cycles),
                subpops=subpops,qualities=qualities,
                num.subpops=num.subpops,best.splits=best.splits))
}

#' Convert subpopulation vector to a list of indices
#' 
#' A helper function to convert a vector of subpopulation identifications into a
#' list appropriate for \code{\link{recSplitConnMat}},
#' \code{\link{qualitySubpops}}, etc.
#' 
#' Note that subpopulations list will be ordered according to the numerical
#' order of the subpopulation indices in the original matrix, which will not
#' necessarily be that of the spatial order of sites in the original
#' connectivity matrix.
#' 
#' @param x vector of subpopulation identifications
#'   
#' @return A list where each element is a vector of indices for a given 
#'   subpopulation.
#'   
#' @seealso See also \code{\link{recSplitConnMat}}, \code{\link{qualitySubpops}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
subpopsVectorToList <- function(x) {
    xx = sort(unique(x))

    ta = list()
    for (yy in xx)
        ta[[length(ta)+1]] = which( x == yy )
    return(ta)
}
