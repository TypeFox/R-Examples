#' Create surveillance data simulator.
#'
#' \code{create_simulator} creates a pomp object that will run
#' simulations of an SIR or SIS model according to Gillespie's direct
#' method and generate simulated observations of the process.
#'
#' See the vignette "Getting started with spaero" for a description of
#' the model. The "params" argument must include all model
#' parameters. These will become the default parameters for the model
#' object. They can be overriden when the simulation is run via the
#' "params" argument of \code{pomp::simulate}. The case is the same
#' for the "times" argument. The "covar" argument should be a data
#' frame with a column named for each of the time-dependent parameters
#' and a column named time. This data frame describes the time series
#' of each of the time-dependent parameters. In the simulation,
#' interpolation based on this data frame determines the value of
#' these parameters at specific instants in time. The user must ensure
#' that these values result in the parameters remaining non-negative
#' for the course of the simulation.
#'
#' @return A pomp object with which simulations can be run via
#' \code{pomp::simulate}.
#' @param times A numeric vector of increasing times at which the
#' state of the simulation will be sampled.
#' @param t0 The time at which the simulation is started with state
#' variable set to the initial conditions specified via params.
#' @param process_model Character string giving the process
#' model. Allowed values are '"SIR"' and '"SIS"'.
#' @param params A named numeric vector of parameter values and
#' initial conditions.
#' @param covar A data frame containing values of the time-dependent
#' components of the parameters.
#'
#' @seealso \code{\link[pomp]{pomp}} for documentation of pomp objects
#' @useDynLib spaero
#' @export
#' @examples
#'
#' foo <- create_simulator()
#' out <- pomp::simulate(foo, times=seq(0, 20, by=1/26))
#' out <- as(out, "data.frame")
#' head(out)
#'
#' opar <- par(mfrow=c(2, 1))
#' plot((S/N)~time, data=out, type="l")
#' plot(cases~time, data=out, type="l")
#' par(opar)
#'
create_simulator <- function(times=seq(0, 9), t0=min(times),
                             process_model=c("SIR", "SIS"),
                             params=c(gamma=24, mu=1 / 70, d=1 / 70, eta=1e-5,
                                 beta=1e-4, rho=0.1, S_0=1, I_0=0, R_0=0,
                                 N_0=1e5),
                             covar=data.frame(gamma_t=c(0, 0), mu_t=c(0, 0),
                                 d_t=c(0, 0), eta_t=c(0, 0), beta_t=c(0, 0),
                                 time=c(0, 1e6))) {
  process_model <- match.arg(process_model)
  if (!requireNamespace("pomp", quietly = TRUE)) {
    stop(paste("The pomp package is needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
  data <- data.frame(time=times, reports=NA)
  d <- cbind(birth=c(1,1,1,1,0),
             sdeath=c(1,1,1,1,0),
             infection=c(1,1,1,1,0),
             ideath=c(1,1,1,1,0),
             recovery=c(1,1,1,1,0),
             rdeath=c(1,1,1,1,0))
  if (process_model == "SIR") {
    v <- cbind(birth=c(1,0,0,1,0),
               sdeath=c(-1,0,0,-1,0),
               infection=c(-1,1,0,0,0),
               ideath=c(0,-1,0,-1,0),
               recovery=c(0,-1,1,0,1),
               rdeath=c(0,0,-1,-1,0))
  } else {
    v <- cbind(birth=c(1,0,0,1,0),
               sdeath=c(-1,0,0,-1,0),
               infection=c(-1,1,0,0,0),
               ideath=c(0,-1,0,-1,0),
               recovery=c(1,-1,0,0,1),
               rdeath=c(0,0,-1,-1,0))
  }
  rprocess <- pomp::gillespie.sim(rate.fun="_transition_rates",
                                  PACKAGE="spaero", v=v, d=d)
  initializer <- function(params, t0, ...) {
    comp.names <- c("S", "I", "R")
    ic.names <- c("S_0", "I_0", "R_0")
    x0 <- stats::setNames(numeric(5), c("S", "I", "R", "N", "cases"))
    fracs <- params[ic.names]
    x0["N"] <- params["N_0"]
    x0[comp.names] <- round(params["N_0"] * fracs / sum(fracs))
    if(params["rho"] < 0 | params["rho"] > 1) {
      stop("rho must be in [0, 1]", call.=FALSE)
    }
    pos.names <- c("S_0", "I_0", "R_0", "N_0")
    if(any(params[pos.names] < 0)) {
      stop(paste("All", paste(pos.names, collapse=" "), "should be >= 0."),
           call.=FALSE)
    }
    x0
  }
  pomp::pomp(data=data, times="time", t0=t0, params=params, rprocess=rprocess,
             measurement.model=reports~binom(size=cases, prob=rho),
             covar=covar, statenames=c("S", "I", "R", "N", "cases"),
             paramnames=c("gamma", "mu", "d", "eta", "beta", "rho", "S_0",
                 "I_0", "R_0", "N_0"),
             covarnames=c("gamma_t", "mu_t", "d_t", "eta_t", "beta_t"),
             tcovar="time", zeronames="cases", initializer=initializer)
}
