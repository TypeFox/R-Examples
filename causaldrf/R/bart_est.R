##' The BART estimator
##'
##' This function estimates the ADRF using Bayesian additive regression trees (BART).
##'
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param outcome_formula is the formula used for fitting the outcome surface.
##' gps is one of the independent variables to use in the outcome_formula. ie.
##' \code{Y ~ treat + X.1 + X.2 + ...} or a variation of this.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param grid_val contains the treatment values to be evaluated.
##' @param ... additional arguments to be passed to the bart() outcome function.
##'
##' @details
##' BART is a prediction model that is applicable to many settings, one of which
##'  is causal inference problems.  It is a sum of trees fit, but the influence
##'   of each tree is held back by a regularization prior so that each tree only
##'    contributes a small amount to the overall fit.  Priors are put on the
##'    parameters to avoid overfitting the data and so that no single tree has
##'     a significant influence on the model fit.
##'      For more details see Chipman (2010).
##'
##' BART does not require fitting a treatment model.  Instead, it fits a
##' response surface to the whole dataset and if the response surface is
##' correctly specified, then the causal effect estimate is unbiased.
##' Although most of the focus on BART is for the binary treatment setting,
##'  Hill (2011) also mentions an extension to the continuous or
##'   multidose treatment setting.  When using BART in this continuous treatment
##'    setting, Hill (2011) compares the outcomes of units with
##'    treatment level \eqn{T_i = t} to their outcomes had \eqn{T_i = 0}.
##'    This method infers the treatment effect of units had they not received
##'    treatment compared to their actual observed treatment.  The comparison
##'    is between \eqn{Y_i(0)| (I = 1, T_i = t)} and \eqn{Y_i(t)| (I = 1, T_i = t)}
##'    where \eqn{I = 1} means that the unit is part of the treatment group.
##'    The causal effect is comparing the predicted outcome of units that
##'    received treatment with what their predicted outcome would have been
##'    had they received zero treatment.
##'
##'    This method performs well in simulation studies.
##'    One drawback from BART is the amount of computing time needed.
##'
##' @return \code{bart_est} returns an object of class "causaldrf_simple",
##' a list that contains the following components:
##' \item{param}{parameter estimates for a bart fit.}
##' \item{out_mod}{the result of the bart fit.}
##' \item{call}{the matched call.}
##'
##' @seealso \code{\link{nw_est}}, \code{\link{iw_est}}, \code{\link{hi_est}}, \code{\link{gam_est}},
##' \code{\link{add_spl_est}}, etc. for other estimates.
##'
##' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
##' for use in the different estimates.
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' Hill, Jennifer L. (2011). Bayesian nonparametric modeling for causal
##' inference. \emph{Journal of Computational and Graphical Statistics}
##' \bold{20.1} (2011).
##'
##' Chipman, Hugh A and George, Edward I and McCulloch, Robert E and others (2010).
##' BART: Bayesian additive regression trees.
##' \emph{The Annals of Applied Statistics}
##' \bold{4.1}, 266--298.
##'
##' @examples
##'
##' ## Example from Schafer (2015).  bart takes a few minutes to run (depending on computer).
##'
##' example_data <- sim_data
##'
##' \dontrun{
##' # This estimate takes a long time to run...
##' bart_list <- bart_est(Y = Y,
##'           treat = T,
##'           outcome_formula = Y ~ T + B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'           data = example_data,
##'           grid_val = seq(8, 16, by = 1))
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'     example_data$Y[sample_index],
##'     xlab = "T",
##'     ylab = "Y",
##'     main = "bart estimate")
##'
##' lines(seq(8, 16, by = 1),
##'       bart_list$param,
##'       lty = 2,
##'       lwd = 2,
##'       col = "blue")
##'
##' legend('bottomright',
##'         "bart estimate",
##'         lty=2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y',
##'         cex=1)
##' }
##'
##' rm(example_data, bart_list, sample_index)
##'
##'
##' @usage
##' bart_est(Y,
##'          treat,
##'          outcome_formula,
##'          data,
##'          grid_val,
##'          ...)
##'
##' @export
##'
##'


bart_est <- function(Y,
                     treat,
                     outcome_formula,
                     data,
                     grid_val,
                     ...){

  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # outcome_formula is the formula for the outcome surface
  # data will contain all the data: X, treat, and Y
  # grid_val is the set of grid points on T

  # This function returns a list containing the estimated ADRF,
  # and the object from the bart fit.

  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("outcome_formula" %in% names(tempcall))) stop("No outcome_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("grid_val" %in% names(tempcall)))  stop("No grid_val specified")

  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(Y = data[, as.character(tempcall$Y)])
  tempdat$treat <- data[,as.character(tempcall$treat)]

  # utils::str(m_frame <- model.frame(tempcall$outcome_formula, data))
  m_frame <- model.frame(tempcall$outcome_formula, data)
  covar_matrix_temp <- model.matrix(eval(tempcall$outcome_formula), m_frame)
  covar_mat <- covar_matrix_temp[, -1]



#-------------------------------------------------------------------------------

# create dataset
use=cbind(tempdat$Y, covar_mat)

# xt is everything except for the Y vector
# xt=as.matrix(na.omit(use)[,-1])   # not sure if I need to make xp1 as.matrix?
xt=na.omit(use)[,-1]

nt=nrow(xt)

y=as.numeric(na.omit(use)[,1])

# length_grid is the number of gridpoints
length_grid <- length(grid_val)
xp1=use[,-1]
# xp1=as.matrix(use[,-1])  # not sure if I need to make xp1 as.matrix?


# create list for dealing with arbitrary gridpoints
# LL is a list that contains "length_grid" number of objects.
# each object in the list is a copy of xp1, but with the treatment being set
# to one gridpoint value for the each object.
# each object will be used to predict the average outcome at one gridpoint.
# A mean Y value is calculated from each object in the list.
# mean_bart contains the average Y value at each grid point.
LL <- list(0)
(for (i in 1:length_grid){
  LL <- c(LL, as.name(paste("xp", i, sep = "_")))})
LL[[1]] <- NULL

# this loop sets each object in the list to have treatment value equal to a grid point value.
for (i in 1:length_grid){
LL[[i]] <- xp1
LL[[i]][,1] <- grid_val[i]
}

# this do.call takes a list and makes it into a single matrix
xp <- do.call("rbind", LL)

# calculate the bart components
bart.tot <- BayesTree::bart(x.train=xt,   y.train=y,  x.test=xp, ...)

#### results
index_start <- seq(1, length_grid*nt, by = nt)
index_end <- seq(nt, length_grid*nt, by = nt)


mean_bart <- numeric(length_grid)
for (i in 1:length(mean_bart)){
  mean_bart[i] <- mean(bart.tot$yhat.test.mean[index_start[i]:index_end[i]])

}


z_object <- list(param = mean_bart,
                 out_mod = bart.tot,
                 call = tempcall)

class(z_object) <- "causaldrf_simple"
z_object

}

