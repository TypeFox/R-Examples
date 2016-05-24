#' RDD coefficient prediction
#'
#' Function to predict the RDD coefficient in presence of covariate (without covariates, returns the same than \code{\link{rdd_coef}})
#' @param object A RDD regression object
#' @param covdata New data.frame specifying the values of the covariates, can have multiple rows. 
#' @param se.fit A switch indicating if standard errors are required.
#' @param vcov. Specific covariance function (see package sandwich ), by default uses the \code{\link{vcov}}
#' @param newdata Another data on which to evaluate the x/D variables. Useful in very few cases. 
#' @param stat The statistic to use if there are multiple predictions, 'identity' just returns the single values, 'mean' averages them 
#' @param weights Eventual weights for the averaging of the predicted values. 
#' @details The function \code{rdd_pred} does a simple prediction of the RDD effect
#'  \deqn{RDDeffect= \mu(x, z, D=1) - \mu(x, z, D=0)}
#' When there are no covariates (and z is irrelevant in the equation above), this amounts exactly to the usual RDD coefficient, 
#' shown in the outputs, or obtained with \code{\link{rdd_coef}}. If there were covariates, and if these covariates were estimated using the 
#' \dQuote{include} \emph{strategy} and with different coefficients left and right to the cutoff (i.e.
#' had argument \emph{slope} = \dQuote{separate}), than the RDD effect is also dependent on the value of the covariate(s). 
#' \code{rdd_pred} allows to set the value of the covariate(s) at which to evaluate the RDD effect, by providing a data.frame with
#' the values for the covariates. Note that the effect can be evaluated at multiple points, if you provide multiple rows of \code{covdata}. 
#'
#' In pressence of covariate-specific RDD effect, one may wish to estimate an average effect. This can be done by setting the argument \code{stat='mean'}. 
#' Weights can additionally be added, with the argument \code{weights}, to obtain a weighted-average of the predictions. Note however that in most cases, 
#' this will be equivalent to provide covariates at their (weighted) mean value, which will be much faster also!
#'
#' Standard errors, obtained setting the argument \code{se.fit=TRUE}, are computed using following formula:
#'  \deqn{x_i \Omega x_i^{'}}
#' where \eqn{\Omega} is the estimated variance-covariance matrix ( by default \eqn{\sigma^2(X^{'}X)^{-1}} using \code{\link{vcov}}) and 
#' \eqn{x_i} is the input data (a mix of covdata and input data). If one wishes individual predictions, standard errors are simply obtained 
#' as the square of that diagonal matrix, whereas for mean/sum, covariances are taken into account. 
#' @return Returns the predicted value(s), and, if se.fit=TRUE, their standard errors. 
#' @export
#' @references Froehlich (2007) Regression discontinuity design with covariates, IZA discussion paper 3024
#' @examples
#' # Load data, add (artificial) covariates:
#' data(house)
#' n_Lee <- nrow(house)
#' z1 <- runif(n_Lee)
#' house_rdd <- rdd_data(y=y, x=x, data=house, covar=z1, cutpoint=0)
#' 
#' # estimation without covariates: rdd_pred is the same than rdd_coef:
#' reg_para <- rdd_reg_lm(rdd_object=house_rdd)
#' 
#' rdd_pred(reg_para)
#' rdd_coef(reg_para, allInfo=TRUE)
#' 
#' # estimation with covariates: 
#' reg_para_cov <- rdd_reg_lm(rdd_object=house_rdd,
#'                           covariates='z1',
#'                           covar.opt=list(slope='separate') )
#'
#' # should obtain same result as with RDestimate                             
#' rdd_pred(reg_para_cov, covdata=data.frame(z1=0)) 
#'   
#' # evaluate at mean of z1 (as comes from uniform)
#' rdd_pred(reg_para_cov, covdata=data.frame(z1=0.5))

rdd_pred <- function(object, covdata, se.fit = TRUE, vcov. = NULL, newdata, stat = c("identity", "sum", "mean"), weights) {
    
    stat <- match.arg(stat)
    
    if (!missing(weights)) {
        if (missing(covdata)) 
            stop("Arg 'weights' only useful with arg 'covdata'")
        if (stat == "identity") 
            stop("Argument 'weights' not useful when arg: stat='identity'")
        if (stat == "sum") {
            warning("Providing weights for a sum makes little sense?!")
        }
        if (length(weights) != NROW(covdata)) 
            stop("Weights should be of the same length than covdata")
    }
    
    x_call <- getCall(object)
    hasCo <- hasCovar(object)
    
    if (is.null(x_call$covar.opt)) {
        covar.slope <- "same"
        covar.strat <- "include"
    } else {
        covar.slope <- ifelse(is.null(x_call$covar.opt$slope), "same", x_call$covar.opt$slope)
        covar.strat <- ifelse(is.null(x_call$covar.opt$strategy), "include", x_call$covar.opt$strategy)
    }
    
    
    ## get original data structure:
    mf <- model.frame(object)[1:2, -1]
    if (any(grepl("\\(weights\\)", colnames(mf)))) 
        mf <- mf[, -grep("\\(weights\\)", colnames(mf))]
    
    ## Fill orig struc with 0/1
    if (missing(newdata)) {
        which.D <- grep("^D$", colnames(mf))
        mf[, which.D] <- c(0, 1)  ## set coeff of interest
        mf[, -which.D] <- 0  ## remove others (not absolutely necessary actually)
        newdata <- mf
    }
    
    ## Merge covdata with newdata:
    
    if (!missing(covdata)) {
        if (covar.strat == "residual") 
            stop("Do not provide 'covdata' if covariates were use with 'residual' strategy")
        if (covar.slope == "separate") {
            Nrow_cov <- nrow(covdata)
            if (Nrow_cov > 1) 
                newdata <- newdata[c(1, rep(2, Nrow_cov)), ]
            if (!is.null(rownames(covdata))) {
                if ("1" %in% rownames(covdata)) 
                  rownames(newdata)[1] <- "0"
                rownames(newdata)[-1] <- rownames(covdata)
            } else {
                rownames(newdata) <- c(0, seq_len(Nrow_cov))
            }
            colnames_cov <- colnames(covdata)
            ind <- seq(from = 2, by = 2, length.out = Nrow_cov)
            if (!all(colnames_cov %in% colnames(newdata))) 
                stop("Arg 'covdata' contains colnames not in the data")
            newdata[2:nrow(newdata), paste(colnames(covdata), "D", sep = ":")] <- covdata
        }
    }
    
    multiN <- nrow(newdata) > 2
    
    ## Merge and check no NAs
    X_i <- as.matrix(cbind(1, newdata))
    if (any(is.na(X_i))) {
        warning("data contains NA. Were removed")
        X_i <- X_i[-apply(X_i, 1, function(x) any(is.na(x))), ]
    }
    
    ## Set up variance matrix: X_i (X'X)^{-1} X_i'
    if (is.null(vcov.)) 
        vcov. <- vcov(object)
    X_inv <- vcov.
    mat <- X_i %*% X_inv %*% t(X_i)
    
    ## preds:
    
    if (!multiN) {
        pred_point <- drop(diff(X_i %*% rdd_coef(object, allCo = TRUE)))
        if (se.fit) 
            pred_se <- sqrt(sum(c(diag(mat), -2 * mat[1, 2])))
    } else {
        d <- X_i %*% coef(object)
        
        
        Mat_SUM <- cbind(1, diag(nrow(d) - 1))
        Mat_DIAG <- matrix(diag(mat), ncol = 1)
        if (missing(weights)) {
            MAT_SmallSum <- matrix(c(-(nrow(d) - 1), rep(1, nrow(d) - 1)), nrow = 1)  ## create vector: [- n-1, 1, 1, 1....]
        } else {
            MAT_SmallSum <- matrix(c(-1, weights), nrow = 1)  ## create vector: [- 1, w_1, w_2, w_n]
        }
        
        if (stat == "identity") {
            Mat_DIFF <- Mat_SUM
            Mat_DIFF[, 1] <- -1
            pred_point <- drop(Mat_DIFF %*% d)
            if (se.fit) 
                pred_se <- drop(sqrt(Mat_SUM %*% Mat_DIAG - 2 * mat[1, 2:ncol(mat)]))
        } else {
            if (stat == "mean" & missing(weights)) 
                MAT_SmallSum <- MAT_SmallSum/Nrow_cov
            pred_point <- drop(MAT_SmallSum %*% d)
            if (se.fit) 
                pred_se <- drop(sqrt(MAT_SmallSum %*% mat %*% t(MAT_SmallSum)))
        }
    }
    
    
    ## result:
    if (se.fit) {
        res <- list()
        res$fit <- pred_point
        res$se.fit <- pred_se
    } else {
        res <- pred_point
    }
    res
} 
