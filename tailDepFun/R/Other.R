#' Empirical stable tail dependence function
#'
#' Returns the stable tail dependence function in dimension \code{d}, evaluated in a point \code{cst}.
#'
#' @param ranks A \code{n} x \code{d} matrix, where each column is a permutation of the integers \code{1:n}, representing the ranks computed from a sample of size \code{n}.
#' @param k An integer between 1 and \eqn{n - 1}; the threshold parameter in the definition of the empirical stable tail dependence function.
#' @param cst The value in which the tail dependence function is evaluated: defaults to \code{rep(1,d)}, i.e., the extremal coefficient.
#' @return A scalar between \eqn{\max(x_1,\ldots,x_d)} and \eqn{x_1 + \cdots + x_d}.
#' @seealso \code{\link{stdfEmpCorr}}
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @export
#' @examples
#' ## Simulate data from the Gumbel copula and compute the extremal coefficient in dimension four.
#' set.seed(2)
#' cop <- copula::gumbelCopula(param = 2, dim = 4)
#' data <- copula::rCopula(n = 1000, copula = cop)
#' stdfEmp(apply(data,2,rank), k = 50)
stdfEmp <- function(ranks, k, cst = rep(1,ncol(ranks))) {
    d <- ncol(ranks)
    n <- nrow(ranks)
    if(length(cst) != ncol(ranks)){
        warning("The length of cst should be equal to ncol(ranks)")
    }
    res<-.C("ecGeneral",as.integer(as.vector(ranks)),as.integer(d),
                  as.double(cst),as.integer(n),as.integer(k),result=double(1),
            PACKAGE = "tailDepFun")$result
    return(res)
}

#' Bias-corrected empirical stable tail dependence function
#'
#' Returns the bias-corrected stable tail dependence function in dimension \code{d}, evaluated in a point \code{cst}.
#'
#' @param ranks A \code{n} x \code{d} matrix, where each column is a permutation of the integers \code{1:n}, representing the ranks computed from a sample of size \code{n}.
#' @param k An integer between 1 and \eqn{n - 1}; the threshold parameter in the definition of the empirical stable tail dependence function.
#' @param cst The value in which the tail dependence function is evaluated: defaults to \code{rep(1,d)}.
#' @param tau The parameter of the power kernel. Defaults to 5.
#' @param k1 An integer between 1 and \eqn{n}; defaults to \eqn{n} - 10.
#' @return A scalar between \eqn{\max(x_1,\ldots,x_d)} and \eqn{x_1 + \cdots + x_d}.
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @references Beirlant, J., Escobar-Bach, M., Goegebeur, Y., and Guillou, A. (2016). Bias-corrected estimation of stable tail dependence function. Journal of Multivariate Analysis, 143, 453-466.
#' @seealso \code{\link{stdfEmp}}
#' @details The values for \code{k1} and \code{tau} are chosen as recommended in Beirlant et al. (2016). This function might be slow for large \code{n}.
#' @export
#' @examples
#' ## Simulate data from the Gumbel copula
#' set.seed(2)
#' cop <- copula::gumbelCopula(param = 2, dim = 4)
#' data <- copula::rCopula(n = 1000, copula = cop)
#' stdfEmpCorr(apply(data,2,rank), k = 50)
stdfEmpCorr<-function(ranks,k,cst = rep(1,ncol(ranks)),tau = 5,k1 = (nrow(ranks) - 10)){
    a <- r <- 0.4
    rhot<-rhoL(ranks,k1,cst,tau,a,r)
    temp1<-sapply(c(1:k)/(k+1), function(i) Kernel(i,tau)*(i^(-rhot)))
    temp2<-sapply(c(1:k)/(k+1), function(i) Kernel(i,tau))
    alphat<-alphaL(ranks,k1,cst,-rhot)
    res<-tildeL(ranks,k,tau,cst) - ((k1/k)^rhot)*alphat*(sum(temp1)/k)
    finalval<-res/(sum(temp2)/k)
    return(ifelse(finalval < max(cst), max(cst), ifelse(finalval > sum(cst), sum(cst), finalval)))
}

#' Integrated empirical stable tail dependence function
#'
#' Analytical implementation of the integral of the bivariate stable tail dependence function over the unit square.
#'
#' @param ranks A \code{n} x 2 matrix, where each column is a permutation of the integers \code{1:n}, representing the ranks computed from a sample of size \code{n}.
#' @param k An integer between 1 and \eqn{n - 1}; the threshold parameter in the definition of the empirical stable tail dependence function.#' @return A scalar.
#' @export
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @examples
#' ranks <- cbind(sample(1:20), sample(1:20))
#' stdfEmpInt(ranks, k = 5)
stdfEmpInt <- function(ranks, k) {
    n <- nrow(ranks)
    res<-.C("stdf",as.integer(ranks[, 1]),as.integer(ranks[, 2]),as.integer(n),
            as.integer(k),result=double(1), PACKAGE = "tailDepFun")$result
    return(res)
}

#' Selects a grid of indices.
#'
#' Returns a regular grid of indices in which to evaluate the stable tail dependence function.
#'
#' @param cst A vector containing the values used to construct the grid. Must contain 0.
#' @param d An integer, representing the dimension.
#' @param nonzero An vector containing integers between \eqn{2} and \eqn{d}, representing the number of non-zero elements in every row of the grid. Defaults to \eqn{2}.
#' @param locations A \code{d} x \eqn{2} matrix containing the Cartesian coordinates of \eqn{d} points in the plane. Used for the Brown-Resnick process only. If not \code{NULL}, then \code{cst} must be \code{c(0,1)} and \code{nonzero} must be \code{2}.
#' @param maxDistance If \code{locations} is not \code{NULL}, pairs of locations with distance not larger than \code{maxDistance} will be selected.
#' @return A matrix with \code{q} rows and \code{d} columns, where every row represents a vector in which we will evaluate the stable tail dependence function (for the weighted least squares estimator) or where every row indicates which pairs of variables to use (for the M-estimator)
#' @export
#' @examples
#' selectGrid(cst = c(0,0.5,1), d = 3, nonzero = c(2,3))
#' locations <- cbind(rep(1:3, each = 3), rep(1:3,3))
#' selectGrid(c(0,1), d = 9, locations = locations, maxDistance = 1.5)
selectGrid <- function(cst, d, nonzero = 2, locations = NULL, maxDistance = 10^6){
    if((! 0 %in% cst)){
        warning("cst needs to contain zero.")
    } else if(! is.whole(d)){
        warning("d should be an integer")
    } else if((!is.vector(nonzero)) || !all(is.whole(nonzero))){
        warning("nonzero should be a vector of integers")
    }
    if(length(cst) > 2 && !is.null(locations)){
        warning("If locations is not NULL, cst is set to c(0,1)")
        cst <- c(0,1)
    }
    if((!all(cst %in% c(0,1))) && !is.null(locations)){
        warning("If locations is not NULL, cst is set to c(0,1)")
        cst <- c(0,1)
    }
    if(!is.null(locations) && (length(nonzero) > 1 || (length(nonzero)==1 && nonzero !=2))){
        warning("If locations is not NULL, nonzero is set to 2")
        nonzero <- 2
    }
    candidate <- lapply(nonzero, function(i) t(utils::combn(d, i)))
    q <- sum(unlist(lapply(candidate, function(i) nrow(i))))
    qtot <- lapply(candidate, function(i) nrow(i)*((length(cst)-1)^ncol(i)))
    grid <- matrix(0,ncol=d,nrow=sum(unlist(qtot)))
    count <- 1
    for(i in 1:length(nonzero)){
        for(j in 1:nrow(candidate[[i]])){
            temp <- as.matrix(expand.grid(lapply(c(1:nonzero[i]), function(I) cst[cst!=0])))
            colnames(temp) <- NULL
            grid[count:(count + nrow(temp) - 1),candidate[[i]][j,]] <- temp
            count <- count + nrow(temp)
        }
    }
    if(is.null(locations)){
        return(grid)
    } else{
        if(d != nrow(locations)){
            warning("dimension is set to nrow(locations)")
            d <- nrow(locations)
        }
        result <- apply(grid, 1, function(i){
            temp <- locations[i == TRUE,][2,] - locations[i == TRUE,][1,]
            if(sqrt(temp[1]^2 + temp[2]^2) <= maxDistance + 1e-10){
                return(TRUE)
            } else{
                return(FALSE)
            }
        })
        return(grid[result,])
    }
}

#' Wind speeds in the Netherlands.
#'
#' Daily maximal speeds of wind gusts, measured in 0.1 m/s. The data are observed at
#' 22 inland weather stations in the Netherlands. Only the summer months are presented
#' here (June, July, August). Also included are the Euclidian coordinates of the 22
#' weather stations, where a distance of 1 corresponds to 100 kilometers.
#'
#' @format dataKNMI$data is a matrix with 672 rows and 22 columns, dataKNMI$loc is a matrix with 22 rows
#' and 2 columns.
#' @source KNMI
#' @name dataKNMI
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @examples
#' data(dataKNMI)
#' n <- nrow(dataKNMI$data)
#' locations <- dataKNMI$loc
#' x <- apply(dataKNMI$data, 2, function(i) n/(n + 0.5 - rank(i)))
#' indices <- selectGrid(cst = c(0,1), d = 22, locations = locations, maxDistance = 0.5)
#' EstimationBR(x, locations, indices, k = 60, method = "Mestimator", isotropic = TRUE,
#' covMat = FALSE)$theta
NULL

#' EUROSTOXX50 weekly negative log-returns.
#'
#' The first three columns represent the weekly negative log-returns of the index prices
#' of the EUROSTOXX50 and of its subindices correspoding to the supersectors chemicals and insurance.
#' The fourth and fifth columns represent the weekly negative log-returns of the index prices of
#' the DAX and the CAC40 indices.
#' The sixth to tenth columnds represent the weekly negative log-returns of the stock prices of Bayer,
#' BASF, Allianz, AXA, and Airliquide respectively.
#'
#' @format \code{dataEUROSTOXX} is a matrix with 711 rows and 10 columns.
#' @source Yahoo Finance
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @name dataEUROSTOXX
#' @examples
#' data(dataEUROSTOXX)
#' ## Transform data to unit Pareto margins
#' n <- nrow(dataEUROSTOXX)
#' x <- apply(dataEUROSTOXX, 2, function(i) n/(n + 0.5 - rank(i)))
#' ## Define indices in which we evaluate the estimator
#' indices <- selectGrid(c(0,0.5,1), d = 10, nonzero = c(2,3))
#' start <- c(0.67,0.8,0.77,0.91,0.41,0.47,0.25,0.7,0.72,0.19,0.37,0.7,0.09,0.58)
#' ## Estimate the parameters. Lasts up to ten minutes.
#' ## EstimationMaxLinear(x, indices, k = 40, method = "WLS", startingValue = start,
#' ## covMat = FALSE, EURO = TRUE)

NULL

#' tailDepFun
#'
#' The package \code{tailDepFun} provides functions implementing two rank-based minimal distance estimation
#' methods for parametric tail dependence models for distributions attracted to a max-stable law.
#' The estimators, referred to as the pairwise M-estimator and the weighted least squares estimator, are
#' described in Einmahl et al. (2016a) and Einmahl et al. (2016b). Extensive examples to illustrate the use
#' of the package can be found in the accompanying vignette.
#'
#' Currently, this package allows for estimation of the Brown-Resnick process, the Gumbel (or logistic) model
#' and max-linear models (possibly on a directed acyclic graph). The main functions of this package are
#' \code{\link{EstimationBR}}, \code{\link{EstimationGumbel}} and \code{\link{EstimationMaxLinear}},
#' but several other functions are exported as well: \code{\link{stdfEmpInt}}
#' returns the integral of the bivariate empirical stable tail dependence function over the unit square, and
#' \code{\link{stdfEmp}} and \code{\link{stdfEmpCorr}} return the (bias-corrected) empirical stable tail dependence
#' function. The functions \code{\link{AsymVarBR}}, \code{\link{AsymVarGumbel}}, \code{\link{AsymVarMaxLinear}}
#' return the asymptotic covariance matrices of the estimators. An auxiliary function to select a regular
#' grid of indices in which to evaluate the stable tail dependence function is exported as well,
#' \code{\link{selectGrid}}. Finally, two datasets are available: \code{\link{dataKNMI}} (Einmahl et al., 2016a)
#' and \code{\link{dataEUROSTOXX}} (Einmahl et al., 2016b).
#'
#' @name tailDepFun
#' @docType package
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A., and Segers, J. (2016a). An Mestimator of spatial tail dependence. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(1), 275-298.
#' @references Einmahl, J.H.J., Kiriliouk, A., and Segers, J. (2016b). A continuous updating weighted least squares estimator of tail dependence in high dimensions. See http://arxiv.org/abs/1601.04826.
#' @examples
#' ## get a list of all help files of user-visible functions in the package
#' help(package = tailDepFun)
NULL

