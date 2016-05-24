#' reservoir: Tools for Analysis, Design, and Operation of Water Supply Storages
#'
#' Measure single reservoir performance using resilience, reliability, and vulnerability metrics; compute storage-yield-reliability relationships; determine no-fail Rippl storage with sequent peak analysis; optimize release decisions for water supply, hydropower, and multi-objective reservoirs using deterministic and stochastic dynamic programming; evaluate inflow characteristics.
#' 
#' @section Analysis and design functions:
#' The \code{\link{Rippl}} function executes the sequent peak algorithm [Thomas and Burden, 1963] to determine the no-fail storage [Rippl, 1883] for given inflow and release time series.
#' The \code{\link{storage}} function gives the design storage for a specified time-based reliability and yield. Similarly, the \code{\link{yield}} function computes the reliability yield given the storage capacity.
#' The \code{\link{simRes}} function simulates a reservoir under standard operating policy, or using an optimised policy produced by \code{\link{sdp_supply}}.
#' The \code{\link{rrv}} function returns three reliability measures, resilience, and dimensionless vulnerability for given storage, inflow time series, and target release [McMahon et al, 2006]. Users can assume Standard Operating Policy, or can apply the output of \code{\link{sdp_supply}} to determine the RRV metrics under different operating objectives.
#' The \code{\link{Hurst}} function estimates the Hurst coefficient [Hurst, 1951] for an annualized inflow time series, using the method proposed by Pfaff [2008].
#' @section Optimization functions:
#' The Dynamic Programming functions find the optimal sequence of releases for a given reservoir. The Stochastic Dynamic Programming functions find the optimal release policy for a given reservoir, based on storage, within-year time period and, optionally, current-period inflow. 
#' For single-objective water supply reservoirs, users may specify a loss exponent parameter for supply deficits and then optimize reservoir release decisions to minimize summed penalty costs over the operating horizon. This can be done using \code{\link{dp_supply}} or \code{\link{sdp_supply}}. There is also an option to simulate the output of \code{\link{sdp_supply}} using the \code{\link{rrv}} function to validate the policy under alternative inflows or analyze reservoir performance under different operating objectives.
#' The optimal operating policy for hydropower operations can be found using \code{\link{dp_hydro}} or \code{\link{sdp_hydro}}. The operating target is to maximise total energy output over the duration of the input time series of inflows.
#' The \code{\link{dp_multi}} and \code{\link{sdp_multi}} functions allow users to optimize for three weighted objectives representing water supply deficit, flood control, and amenity.
#' @section Storage-depth-area relationships:
#' All reservoir analysis and optimization functions, with the exception of \code{\link{Rippl}}, \code{\link{storage}}, and \code{\link{yield}}, allow the user to account for evaporation losses from the reservoir surface. The package incorporates two storage-depth-area relationships for adjusting the surface area (and therefore evaporation potential) with storage.
#' The simplest is based on the half-pyramid method [Liebe et al, 2005], requiring the user to input the surface area of the reservoir at full capacity via the \code{surface_area} parameter.
#' A more nuanced relationship [Kaveh et al., 2013] is implemeted if the user also provides the maximum depth of the reservoir at full capacity via the \code{max_depth} parameter.
#' Users must use the recommended units when implementing evaporation losses.
#' @section Stochastic generation of synthetic streamflow replicates:
#' The \code{\link{dirtyreps}} function provides quick and dirty generation of stochastic streamflow replicates (seasonal input data, such as monthly or quarterly, only).
#' Two methods are available: the non-parametric kNN bootstrap [Lall and Sharma, 1996] and the parametric periodic Autoregressive Moving Average (PARMA).
#' The PARMA is fitted for p = 1 and q = 1, or PARMA(1,1). Fitting is done numerically by the least-squares method [Salas and Fernandez, 1993].
#' When using the PARMA model, users do not need to transform or deseasonalize the input data as this is done automatically within the algorithm.
#' The kNN bootstrap is non-parametric, so no intial data preparation is required here either.
#' @docType package
#' @name reservoir
#' @examples # 1. Express the distribution of Rippl storage for a known inflow process...
#' layout(1:4)
#' # a) Assume the inflow process follows a lognormal distribution
#' # (meanlog = 0, sdlog = 1):
#' x <- rlnorm(1200)
#' 
#' # b) Convert to a 100-year, monthly time series object beginning Jan 1900
#' x <- ts(x, start = c(1900, 1), frequency = 12)
#' 
#' # c) Begin reservoir analysis... e.g., compute the Rippl storage
#' x_Rippl <- Rippl(x, target = mean(x) * 0.9)
#' no_fail_storage <- x_Rippl$Rippl_storage
#' 
#' # d) Resample x and loop the procedure multiple times to get the
#' # distribution of no-failure storage for the inflow process assuming 
#' # constant release (R) equal to 90 percent of the mean inflow.
#' no_fail_storage <- vector("numeric", 100)
#' for (i in 1:length(no_fail_storage)){
#'   x <- ts(rlnorm(1200), start = c(1900, 1), frequency = 12)
#'   no_fail_storage[i] <- Rippl(x, target = mean(x) * 0.9 ,plot = FALSE)$No_fail_storage
#' }
#' hist(no_fail_storage)
#' 
#' 
#' # 2. Trade off between annual reliability and vulnerability for a given system...
#' layout(1:1)
#' # a) Define the system: inflow time series, storage, and target release.
#' inflow_ts <- resX$Q_Mm3
#' storage_cap <- resX$cap_Mm3
#' demand <- 0.3 * mean(resX$Q_Mm3)
#' 
#' # b) define range of loss exponents to control preference of high reliability
#' # (low loss exponent) or low vulnerability (high loss exponent).
#' loss_exponents <- c(1.0, 1.5, 2)
#' 
#' # c) set up results table
#' pareto_results <- data.frame(matrix(ncol = 2, nrow = length(loss_exponents)))
#' names(pareto_results) <- c("reliability", "vulnerability")
#' row.names(pareto_results) <- loss_exponents
#' 
#' # d) loop the sdp function through all loss exponents and plot results
#' for (loss_f in loss_exponents){
#'  sdp_temp <- sdp_supply(inflow_ts, capacity = storage_cap, target = demand, rep_rrv = TRUE,
#'  S_disc = 100, R_disc = 10, plot = FALSE, loss_exp = loss_f, Markov = TRUE)
#'  pareto_results$reliability[which(row.names(pareto_results)==loss_f)] <- sdp_temp$annual_reliability
#'  pareto_results$vulnerability[which(row.names(pareto_results)==loss_f)] <- sdp_temp$vulnerability
#'  }
#' plot (pareto_results$reliability,pareto_results$vulnerability, type = "b", lty = 3)
#' 
#' @references Hurst, H.E. (1951) Long-term storage capacity of reservoirs, Transactions of the American Society of Civil Engineers 116, 770-808.
#' @references Kaveh, K., H. Hosseinjanzadeh, and K. Hosseini. (2013) A new equation for calculation of reservoir's area-capacity curves, KSCE Journal of Civil Engineering 17(5), 1149-1156.
#' @references Liebe, J., N. Van De Giesen, and Marc Andreini. (2005) Estimation of small reservoir storage capacities in a semi-arid environment: A case study in the Upper East Region of Ghana, Physics and Chemistry of the Earth, 30(6), 448-454.
#' @references Loucks, D.P., van Beek, E., Stedinger, J.R., Dijkman, J.P.M. and Villars, M.T. (2005) Water resources systems planning and management: An introduction to methods, models and applications. Unesco publishing, Paris, France.
#' @references McMahon, T.A., Adeloye, A.J., Zhou, S.L. (2006) Understanding performance measures of reservoirs, Journal of Hydrology 324 (359-382)
#' @references Nicholas E. Graham and Konstantine P. Georgakakos, 2010: Toward Understanding the Value of Climate Information for Multiobjective Reservoir Management under Present and Future Climate and Demand Scenarios. J. Appl. Meteor. Climatol., 49, 557-573.
#' @references Pfaff, B. (2008) Analysis of integrated and cointegrated time series with R, Springer, New York. [p.68]
#' @references Rippl, W. (1883) The capacity of storage reservoirs for water supply, In Proceedings of the Institute of Civil Engineers, 71, 270-278.
#' @references Thomas H.A., Burden R.P. (1963) Operations research in water quality management. Harvard Water Resources Group, Cambridge
#' @references kNN Bootstrap method: Lall, U. and Sharma, A. (1996). A nearest neighbor bootstrap for resampling hydrologic time series. Water Resources Research, 32(3), pp.679-693.
#' @references PARMA method: Salas, J.D. and Fernandez, B. (1993). Models for data generation in hydrology: univariate techniques. In Stochastic Hydrology and its Use in Water Resources Systems Simulation and Optimization (pp. 47-73). Springer Netherlands.

NULL