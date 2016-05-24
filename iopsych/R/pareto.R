# These functions calaculate pareto tradeoffs


# Internal ---------------------------------------------------------------------


#' An internal function for generating criterion weights for XX pareto plots.
#' 
#' @param pts How many different pairs of criterion weights to calculate.
#' @return A matrix of of criterion weights with one column per criteiron. 
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.wtPair
.wtPair <- function(pts) {
    wt_one <- seq(from=0, to=1, by=(1/(pts-1)))
    wt_two <- (1-wt_one)
    return(cbind(wt_one, wt_two))
}

#' Computes data needed for a XX Pareto plot.
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param pts The number of points used. Determines accuracy.
#' @return Points along the pareto optimal surface and the predictor weights
#'         used to compute them. 
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.paretoXX
.paretoXX <- function(r_mat, x_col, y_col, pts=100) {
    #Terms
    rxx <- r_mat[x_col, x_col]
    rxy1 <- r_mat[y_col[1], x_col]
    rxy2 <- r_mat[y_col[2], x_col]
    #Partially Evaluate Regression
    betaFun <- .rmatBetaPE(rxx)
    fn <- function(wt) {
        wts <-c(wt, (1-wt))
        rxy <- (fuseVec(r_mat=r_mat, a=y_col, wt_a=wts))[x_col]
        return(betaFun(rxy=rxy))
    }
    #Find beta weights for each pair of criteria weights
    wt_one <- seq(from=0, to=1, by=(1/(pts-1)))
    beta_mat <- lapply(wt_one, fn)
    #Find Multiple Correlation
    mr_1 <- .solveWtVec(wt=beta_mat, rxx=rxx, rxy=rxy1)
    mr_2 <- .solveWtVec(wt=beta_mat, rxx=rxx, rxy=rxy2)
    mr_mat <- cbind(mr_1, mr_2)
    #Format output
    beta_out <- t(do.call("cbind", beta_mat))
    out <- list(beta_out, wt_one, mr_mat)
    names(out) <- c("betas", "wt_one", "multiple_r")
    return(out)
}


# External --------------------------------------------------------------------


#' Computes data needed for a XX Pareto plot.
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param pts The number of points used. Determines accuracy.
#' @return  \itemize{
#'             \item{betas}{A matrix of beta weights for each criteria weight}
#'             \item{wt_one}{The weight given to the first criterion}
#'             \item{multiple_r}{The correlation between the predictor and 
#'                                criterion composites}
#'              }
#' @author Allen Goebl and Jeff Jones
#' @examples
#' # Setup Data
#' data(dls2007)
#' r_mat <- dls2007[1:6, 2:7]
#'
#' #Run Model
#' XX1 <- paretoXX(r_mat=r_mat, x_col=1:4, y_col=5:6)
#'
#' # Plot Multiple correlations
#' plot(c(0,1), c(.3,.5), type="n", xlab="C1 Wt", ylab="mr") 
#' lines(XX1$wt_one, (XX1$R2)[,1])
#' lines(XX1$wt_one, (XX1$R2)[,2])
#' @export
paretoXX <- function(r_mat, x_col, y_col, pts=100) {
    #Check Input
    .isCorMat(r_mat)
    if(length(y_col) != 2) {stop("y_col should be a vector of length 2")}
    if(max(c(x_col)) > ncol(r_mat)) {stop("x_col is out of bounds")}
    if(max(c(y_col)) > ncol(r_mat)) {stop("y_col is out of bounds")}
    #Call Function
    out <- .paretoXX(r_mat=r_mat, x_col=x_col, y_col=y_col, pts=100)
    #Format Output
    var_names <- colnames(r_mat)[x_col]
    colnames(out[["betas"]]) <- var_names
    return(out)
}


#Pareto XY --------------------------------------------------------------------


#' Computes data needed for a XY Pareto plot.
#'
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param d_vec A vector of d scores.
#' @param gen The number of iterations used by the algorithim.
#' @param pop The population or number of cases used by the algorithim.
#' @param pred_lower The minimum weight allowed for each predictor. 
#' @param pred_upper The maximum weight allowed for each predictor.
#' @return  \itemize{
#'             \item{betas}{A matrix of beta weights for each criteria weight}
#'             \item{mr_d}{A matrix of multiple correlations or d values 
#'                         corresponding to each row of beta weights.}
#'             \item{pareto_optimal}{A vector indicating whether each value is
#'                                   pareto optimal}
#'          }
#' @author Allen Goebl Jeff Jones
#' @examples
#' data(dls2007)
#' dat <- dls2007
#' r_mat <- dat[1:6, 2:7]
#' x_col <- 1:4 
#' y_col <- 5:6
#' d_vec <- dat[1:4, 1]
#'
#' paretoXY(r_mat=r_mat, x_col=1:4, y_col=5, d_vec=d_vec, pred_lower=c(0,0,0,0))
#' paretoXY(r_mat=r_mat, x_col=1:4, y_col=c(5,6))
#' @export
paretoXY <- function(r_mat, x_col, y_col, d_vec, gen=100, pop=100, 
                     pred_lower=rep(-2, length(x_col)),
                     pred_upper=rep( 2, length(x_col))) {
    #Terms
    rxx <- r_mat[x_col, x_col]
    rxy <- r_mat[y_col, x_col]
    rxy_list <- lapply(y_col, function(x) (r_mat[x, x_col]))
    if(!missing(d_vec)) { rxy_list$d <- unname(-d_vec) }
    #Set Objective function
    obj <- function(a) -.solveWtExp(wt=a, rxx=rxx, rxy_list=rxy_list)
    #Optimize Objective
    out <- mco::nsga2(fn           = obj, 
                      idim         = length(x_col), 
                      odim         = length(rxy_list), 
                      generations  = 100, 
                      popsize      = 500,
                      lower.bounds = pred_lower,
                      upper.bounds = pred_upper,
                      vectorized   = TRUE)
    names(out) <- c("betas", "mr_d", "pareto_optimal")
    return(out)
}
