#' Prepare the data
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, 'SID' and co-varaiables(if necessary).
#' @param SID A string indicating the name of the column with the study ID.
#' @param formula.se An optional object of class "formula": A symbolic description of a linear model to be fitted to mean E(x) of sensitivity in the logit scale.
#' the default (when no covariates are included) symbolic description is SID ~ 1 corresponds to the model formula E(x) = mu = exp(a)/(1 + exp(a)) where a is the intercept.
#' When the covariates are categorical and the relative measures are needed it is important to remove the interecept from the model to obtain meaningful parameters. EG for
#' a covariate 'Test' with two levels(A and B) and relative sensitivity of B versus A is needed, then the correct formula is SID ~ Test - 1. See **** .
#' Further information on interpretation of parameters in logistic regression see ***
#' @param formula.sp An optional object of class "formula": A symbolic description of a linear model to be fitted to specificity data.
#' When no covariates are included, the formula is not necessary. By default the covariate information for sensitivity is used.
#' @param formula.omega An optional object of class "formula": A symbolic description of a linear model to be fitted to the copula function.
#' When no covariates are included, the formula is not necessary. By default the covariate information for sensitivity is used.
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
#' @keywords internal

prep.data <- function(data,
                      SID,
                      formula.se,
                      formula.sp,
                      formula.omega){

    if (is.null(data)) stop("Need to supply valid data")

    if (sum(stats::complete.cases(data))!=nrow(data)) print ("Data contains missing values. Only complete cases will be used")

    data <- data[stats::complete.cases(data),]

    if ((sum(grepl('tp', names(data), ignore.case = TRUE)) + sum(grepl('tn', names(data), ignore.case = TRUE)) +
         sum(grepl('fp', names(data), ignore.case = TRUE)) + sum(grepl('fn', names(data), ignore.case = TRUE))) != 4) stop("Data should contain: TP, FN, TN, FP")

    if (is.null(SID)) stop("Need to define study ID")

    data$SID <- data[,grepl(SID, names(data), ignore.case = TRUE)]

    Ns <- nrow(data)

    XSE <- stats::model.matrix(stats::terms(formula.se), data)

    if(ncol(XSE)==1){

        attr(XSE, 'assign') <- NULL
        attr(XSE, 'dimnames') <- NULL

    } else {

        attr(XSE, 'assign') <- NULL
        attr(XSE, 'contrasts') <- NULL
        attr(XSE, 'dimnames') <- NULL
        XSE <- as.matrix(XSE)
    }

    XSP <- stats::model.matrix(stats::terms(formula.sp), data)
    if(ncol(XSP)==1) {
        attr(XSP, 'assign') <- NULL
        attr(XSP, 'dimnames') <- NULL

    } else {

        XSP <- stats::model.matrix(stats::terms(formula.sp), data)

        attr(XSP, 'assign') <- NULL
        attr(XSP, 'contrasts') <- NULL
        attr(XSP, 'dimnames') <- NULL
        XSP <- as.matrix(XSP)
    }

    XOMEGA <- stats::model.matrix(stats::terms(formula.omega), data)
    if(ncol(XOMEGA)==1) {
        XOMEGA <- XSE

    } else {

        XOMEGA <- stats::model.matrix(stats::terms(formula.omega) , data)

        attr(XOMEGA, 'assign') <- NULL
        attr(XOMEGA, 'contrasts') <- NULL
        attr(XOMEGA, 'dimnames') <- NULL
        XOMEGA <- as.matrix(XOMEGA)

    }

    TP <- data[,grepl('tp', names(data), ignore.case = TRUE)]
    TN <- data[,grepl('tn', names(data), ignore.case = TRUE)]
    FP <- data[,grepl('fp', names(data), ignore.case = TRUE)]
    FN <- data[,grepl('fn', names(data), ignore.case = TRUE)]

    if (((!is.numeric(TP)) &  (!is.numeric(FP)) & (!is.numeric(TN)) & (!is.numeric(FN)))) stop("Either (TP, FN, TN, FP) is not numeric")

    Npse <- ncol(XSE)
    Npsp <- ncol(XSP)
    Npomega <- ncol(XOMEGA)

    return(list(data=data, XSE=XSE, XSP=XSP, XOMEGA=XOMEGA, Ns=Ns, Npse=Npse, Npsp=Npsp, Npomega=Npomega))
}

#'  Compute log pointwise predictive density, effective number of parameters and WAIC.
#'
#' @param model stanfit object
#'
#' @author Victoria N Nyaga
#' @keywords internal
#============================== WAIC =====================================================#

waic <- function (model){
    log_sum_exp <- function(x) {
        x_max <- base::max(x)
        x_max + log(base::sum(exp(x - x_max)))
    }
	#Calculate posterior variances from simulation
	colVars <- function (a){
		diff <- a - base::matrix (base::colMeans(a), nrow(a), ncol(a), byrow=TRUE)
		vars <- base::colMeans(diff^2)*base::nrow(a)/(base::nrow(a)-1)
		return (vars)
	}

    log_lik <- rstan::extract(model, "loglik")$loglik
    summands <- apply(log_lik, 2, log_sum_exp)
    correc <- - ncol(log_lik) * log(nrow(log_lik))
    lppd <- sum(summands) + correc
    p_waic_1 <- 2 * (sum(summands - base::colMeans(log_lik)) + correc)
    p_waic_2 <- sum (colVars(log_lik))
    waic_2 <- -2*lppd + 2*p_waic_2
    return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd,
                  p_waic_1=p_waic_1))
}

#============================= DEBYE =================================================#
#'  Compute transform omega to ktau.
#'
#' @param theta correlation parameter(s) from the frank copula function.
#'
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
#' @keywords internal
#




omega.to.ktau <- function(omega){

	debye <- function(omega){

		ft <- function(t){
			t/(exp(t) - 1)
		    }
	(1/omega)*stats::integrate(ft, 0, omega)[[1]]
	    }

        1 + (4*(debye(omega) - 1))/omega
}

#'@examples
#' ceiling(debye(0.0000000001))==1
#
#' data(telomerase)
#'
#' df <- prep.data(data=telomerase,
#'                 SID = "ID",
#'                 formula.se=model1@modelargs$formula.se,
#'                 formula.sp=model1@modelargs$formula.sp,
#'                 formula.omega=model1@modelargs$formula.omega)
#'
#' typeof(df)=="list"


