# sim.ltrc() : Function to generate left truncated and right censored data using the Cox's model

#' Generate left-truncated (and right-cencored) data from the Cox model.
#'
#' Various baseline survival functions and truncation distribution are
#' available. Censoring rate can be designated through tuning the parameter
#' \code{Cmax}; \code{Cmas = Inf} means no censoring.
#'
#' @param n the sample size.
#' @param b a numeric vector for true regression coefficients.
#' @param time.dep logical, whether there is the time-dependent covariate (only
#'   one indicator function Zv = I(t >= zeta) is supported); the default is
#'   FALSE.
#' @param Zv.depA logical, whether the time-dependent covariate \code{Zv}
#'   depends on A^* (the only form supported is Zv = I(t >= zeta + A^*)); the
#'   default is FALSE.
#' @param distr.T the baseline survival time (T*) distribution ("exp" or
#'   "weibull").
#' @param shape.T the shape parameter for the Weibull distribution of T*.
#' @param scale.T the scale parameter for the Weibull distributiof of T*.
#' @param meanlog.T the mean for the log-normal distribution of T*.
#' @param sdlog.T the sd for the log-normal distribution of T*.
#' @param distr.A the baseline truncation time (A*) distribution: either of
#'   \code{"weibull"} (the default), \code{"unif"} (Length-Biased Sampling),
#'   \code{"binomial"} or \code{"dunif"}). Note: If distribution name other than
#'   these are provided, \code{"unif"} will be used.
#' @param shape.A the shape parameter for the Weibull distribution of A*.
#' @param scale.A the scale parameter for the Weibull distribution of A*.
#' @param p.A the success probability for the binomial distribution of A*.
#' @param Cmax the upper bound of the uniform distribution of the censoring time
#'   (C).
#' @param fix.seed an optional random seed for simulation.
#' @importFrom stats rweibull runif rbinom rexp plnorm qlnorm quantile
#' @return a list with a data.frame containing the observed survival times
#'   (\code{Ys}), the observed truncation times (\code{As}), the event indicator
#'   (\code{Ds}) and the covariates (\code{Zs}); a vector of certain quantiles
#'   of Ys (\code{taus}); the censoring proportion (\code{PC}) and the
#'   truncation proportiona (\code{PT}).
#' @examples
#' # With time-invariant covariates only
#' sim1 = sim.ltrc(n = 100)
#' head(sim1$dat)
#' # With one time-dependent covariate
#' sim2 = sim.ltrc(n = 100, time.dep = TRUE,
#'          distr.A = "binomial", p.A = 0.8, Cmax = 5)
#' head(sim2$dat)
#' # With one time-dependent covariate with dependence on the truncation time
#' sim3 = sim.ltrc(n = 100, time.dep = TRUE, Zv.depA = TRUE, Cmax = 5)
#' head(sim3$dat)
#' @export
sim.ltrc = function(n=200, b = c(1,1),
                    time.dep = FALSE, Zv.depA = FALSE,
                    distr.T = "weibull",
                    shape.T = 2, scale.T = 1,
                    meanlog.T = 0, sdlog.T = 1,
                    distr.A = "weibull",
                    shape.A = 1, scale.A = 5, p.A = 0.3,
                    Cmax = Inf, fix.seed = NULL){
	PT = NULL
	PC = NULL
	# get the dimention of the covariates (should be << n)
	p = length(b)
	if( !is.null(fix.seed) ) set.seed(fix.seed)
	j = 1
	k = 0
	Ts = NULL
	As = NULL
	Zs = NULL
	while( j <= n ){
		k = k + 1
		As.j = switch(distr.A,
		              "weibull" = rweibull(1, shape.A, scale.A),
  		            "unif" = runif(1,0,100),
  		            "binomial" = rbinom(1,5,p.A),
                  "dunif" = sample(0:5, 1), rweibull(1, shape.A, scale.A))
		if( time.dep ){
		  # Time at which the time-dependent covariate changes from 0 to 1.
		  zeta.j = rexp(1) + ifelse(Zv.depA, As.j, 0)
		  # Continuous time-invariant covariates
		  Zf.j = runif(p-1, -1, 1)
		  # Relative risks: before = RR0.j; after = RR1.j
		  RR0.j = exp(sum(b[-p] * Zf.j))
		  RR1.j = exp(sum(b[-p] * Zf.j) + b[p])
		  eps.j = rexp(1)
		  # Generate Random Samples from the Cox's model.
		  if( distr.T == "lnorm" ){
		    t0 = -plnorm(zeta.j, meanlog.T, sdlog.T,
		                 lower.tail = FALSE, log.p = TRUE)
		    Ts.j = qlnorm(-(min(eps.j, t0 * RR0.j) / RR0.j +
		                    max(eps.j - t0 * RR0.j, 0) / RR1.j),
		                  meanlog.T, sdlog.T,
		                  lower.tail = FALSE, log.p = TRUE)
		  }else if( distr.T == "weibull" ){
		    t0 = (zeta.j / scale.T) ^ shape.T
		    Ts.j = (min(eps.j, t0 * RR0.j) / RR0.j  +
		            max(eps.j - t0 * RR0.j, 0) / RR1.j)^(1 / shape.T) * scale.T
		  }
		}else{
		  Zf.j = runif(p, -1, 1)
		  zeta.j = NULL
		  eps.j = rexp(1)
		  if( distr.T == "lnorm" ){
		    Ts.j = qlnorm(-eps.j * exp(-sum(b * Zf.j)),
		                  meanlog.T, sdlog.T,
		                  lower.tail = FALSE, log.p = TRUE)
		  }else if( distr.T == "weibull" ){
		    Ts.j = (eps.j * exp(sum(-b * Zf.j)))^(1 / shape.T) * scale.T
		  }
		}
		# Keep only untruncated (A,T)'s
		if( As.j < Ts.j ){
			Ts[j] = Ts.j
			As[j] = As.j
      Zs = rbind(Zs, c(Z = Zf.j, zeta = zeta.j))
			j = j + 1
		}else{
			next
		}
	}
	if( Cmax != Inf ){
		Cs = runif(n, 0, Cmax)
	}else{
    # ghost runs to keep the same random seed (=0/>0 censoring)
    Cs = runif(n)
		Cs = rep(Inf,n)
	}
	Cs = As + Cs
	Ds = as.numeric(Ts <= Cs)
	Ys = pmin(Ts,Cs)
	tau = quantile(Ys,seq(0.1,0.9,0.1))
	dat = data.frame(Zs,As,Ys,Ds)
  # Sort the dataset by observed times; PLAC functions assume this.
	dat = dat[order(dat$Ys), ]
  dat = cbind(ID = 1:n, dat)
	rownames(dat) = NULL
	PC = 1 - mean(dat$Ds)
	PT = 1 - n / k
	return(list(dat = dat, PC = PC, PT = PT, tau = tau))
}

#' Perform the paired log-rank test.
#'
#' Perform the paired log-rank test on the truncation times and the residual survival times to check the stationarity assumption (uniform truncation assumption) of the left-truncated right-censored data.
#'
#' @param dat a data.frame of left-truncated right-censored data.
#' @param A.name the name of the truncation time variable in \code{dat}.
#' @param Y.name the name of the survival time variable in \code{dat}.
#' @param D.name the name of the event indicator in \code{dat}.
#' @return a list containing the test statistic and the p-value of the paired log-rant test.
#' @references Jung, S.H. (1999). Rank tests for matched survival data. \emph{Lifetime Data Analysis, 5(1):67-79}.
#' @examples
#' dat = sim.ltrc(n = 100, distr.A = "weibull")$dat
#' plr(dat)
#' @export
plr = function(dat, A.name = "As", Y.name = "Ys", D.name = "Ds"){
  # N-A est. for A and V
  n = nrow(dat)
  names(dat) = gsub(A.name, "As", names(dat))
  names(dat) = gsub(Y.name, "Ys", names(dat))
  names(dat) = gsub(D.name, "Ds", names(dat))
  dat$Vs = dat$Ys - dat$As
  fit.A = survival::survfit(Surv(As, rep(1, n)) ~ 1, data = dat)
  lambda.A = data.frame(et = fit.A$time, hazard.A = fit.A$n.event / fit.A$n.risk)
  Y.A = stepfun(fit.A$time, c(fit.A$n.risk, 0), right = TRUE)
  fit.V = survival::survfit(Surv(Vs,Ds) ~ 1, data = dat)
  lambda.V = data.frame(et = fit.V$time, hazard.V = fit.V$n.event / fit.V$n.risk)
  Y.V = stepfun(fit.V$time, c(fit.V$n.risk, 0), right = TRUE)
  # combine event times
  et = unique(sort(c(fit.A$time, fit.V$time)))
  lambda.diff = merge(lambda.A, lambda.V, by="et", all = TRUE)
  lambda.diff$hazard.A[is.na(lambda.diff$hazard.A)] = 0
  lambda.diff$hazard.V[is.na(lambda.diff$hazard.V)] = 0
  lambda.diff$h.diff = lambda.diff$hazard.A - lambda.diff$hazard.V
  # Hn(t) for log rank test
  Hn = Y.A(et) * Y.V(et) / (Y.A(et) + Y.V(et)) / n
  # Wn = rank statistics
  Wn = sqrt(n) * sum(Hn * lambda.diff$h.diff)
  # variance of Wn
  epsilon = rep(0,n)
  # matrix notation
  X = c(dat$As, dat$Vs)
  delta = c(rep(1, n), dat$Ds)
  foo = function(x,y)x>=y
  Epsilon = outer(X,X,"foo")
  epsilon = delta[1:n] * colSums(Epsilon[(n + 1):(2 * n), 1:n]) /
                         colSums(Epsilon[ , 1:n]) -
            Epsilon[1:n, ] %*% (delta * colSums(Epsilon[(n + 1):(2 * n), ]) /
                               (colSums(Epsilon))^2)-
            delta[(n + 1):(2 * n)] * colSums(Epsilon[1:n, (n + 1):(2 * n)]) /
            colSums(Epsilon[ ,(n + 1):(2 * n)]) +
            Epsilon[(n + 1):(2 * n), ] %*% (delta * colSums(Epsilon[1:n, ]) /
                                           (colSums(Epsilon))^2)
  epsilon = as.vector(epsilon)
  sigma = sqrt(mean(epsilon^2))
  plr.stat = Wn / sigma
  plr.p = (1 - pnorm(abs(plr.stat))) * 2
  return(list(Stat = plr.stat, P = plr.p))
}

