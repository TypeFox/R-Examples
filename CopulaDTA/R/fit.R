#' Fit copula based bivariate beta-binomial distribution to diagnostic data.
#' @param cdtamodel An object of cdtamodel class from \link{cdtamodel}.
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, 'SID' and co-varaiables(if necessary).
#' @param SID A string indicating the name of the column with the study identifier.
#'@param chains A positive numeric value specifying the number of chains, default is 3.
#'@param iter A positive numeric value specifying the number of iterations per chain. The default is 6000.
#'@param warmup A positive numeric value (<iter) specifying the number of iterations to be discarded(burn-in/warm-up). The default is 1000.
#'@param thin A positive numeric value specifying the interval in which the samples are stored. The default is 10.
#'@param cores A positive numeric values specifying the number of cores to use to execute parallel sampling. When the hardware has more at least 4 cores,
#'the default is 3 cores and otherwise 1 core.
#'@param ... Other optional parameters as specified in \link[rstan]{stan}.
#'@return An object of cdtafit class.
#'

#'@examples
#' \dontrun{
#' fit1 <- fit(model1,
#'                 SID='ID',
#'                 data=telomerase,
#'                 iter=2000,
#'                 warmup=1000,
#'                 thin=1,
#'                 seed=3)
#'
#'
#
#' fit2 <- fit(model2,
#'                 SID='StudyID',
#'                 data=ascus,
#'                 iter=2000,
#'                 warmup=1000,
#'                 thin=1,
#'                 seed=3)
#' }
#'
#'@references {Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Inc.}
#'@references {Clayton DG (1978). A model for Association in Bivariate Life Tables and its Application in
#'Epidemiological Studies of Familial Tendency in Chronic Disease Incidence. Biometrika,65(1), 141-151.}
#'@references {Frank MJ (1979). On The Simultaneous Associativity of F(x, y) and x + y - F(x, y). Aequationes Mathematicae, pp. 194-226.}
#'@references {Farlie DGJ (1960). The Performance of Some Correlation Coefficients for a General Bivariate
#'Distribution. Biometrika, 47, 307-323.}
#'@references {Gumbel EJ (1960). Bivariate Exponential Distributions. Journal of the American Statistical Association, 55, 698-707.}
#'@references {Meyer C (2013). The Bivariate Normal Copula. Communications in Statistics - Theory and Methods, 42(13), 2402-2422.}
#'@references {Morgenstern D (1956). Einfache Beispiele Zweidimensionaler Verteilungen. Mitteilungsblatt furMathematische Statistik, 8, 23 - 235.}
#'@references {Sklar A (1959). Fonctions de Repartition a n Dimensions et Leurs Marges. Publications de l'Institut de Statistique de L'Universite de Paris, 8, 229-231.}
#'@export
#'@importFrom rstan sampling
#'@importFrom rstan stan_model
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
fit.cdtamodel <- function(cdtamodel,
                data,
                SID,
                cores=3,
                chains=3,
                iter=6000,
                warmup=1000,
                thin=10,
                ...) {
#=================================Specify the data=========================================#

    df <- prep.data(data=data,
                    SID = SID,
                    formula.se=cdtamodel@modelargs$formula.se,
                    formula.sp=cdtamodel@modelargs$formula.sp,
                    formula.omega=cdtamodel@modelargs$formula.omega)

    datalist <- list(
        Ns = df$Ns,
        Npse = df$Npse,
        Npsp = df$Npsp,
        Npomega = df$Npomega,
        tp = df$data$TP,
        dis = df$data$TP + df$data$FN,
        tn = df$data$TN,
        nondis = df$data$TN + df$data$FP,
        xse = df$XSE,
        xsp = df$XSP,
        xomega = df$XOMEGA)
#=================================Run the model=========================================#

    #Check available cores, if more than 4, use
    if(parallel::detectCores() < 3)  cores <- 1

    stanmodel <- rstan::stan_model(model_code = cdtamodel@modelcode)


    mod <- rstan::sampling(object = stanmodel,
                           data=datalist,
                           warmup=warmup,
                           thin=thin,
                           chains=chains,
                           cores=cores,
                           iter=iter,
                           ...)

out <- new("cdtafit",
           data=data,
           SID=SID,
           copula=cdtamodel@copula,
           modelargs=cdtamodel@modelargs,
           fit=mod)

return(out)

}

