"smuceR" <-
function(y, x = 1:length(y), x0 = 2 * x[1] - x[2], q = thresh.smuceR(length(y)), alpha, r, lengths, family = c("gauss", "gaussvar", "poisson", "binomial"), param, jumpint = confband, confband = FALSE)
{
  family <- match.arg(family)
  if(family == "binomial" & missing(param)) stop("param must be size of binomial distribution")
  
  if(family == "poisson" | family=="gaussvar")
    param <- NULL
  
  dyadic <- FALSE
  if(missing(lengths) && length(y) <= 1000)
    lengths <- length(y):1
  if(missing(lengths) && length(y) > 1000) {
    lengths <- 2^(floor(log2(length(y))):0)
    dyadic <- TRUE
  }
  
  if(missing(alpha) && is.null(q)) stop("q or alpha must be specified")
  if(!missing(alpha)){
    if(missing(r)) {
      if(dyadic)
	q <- quantile(stepR::MRC.asymptotic.dyadic, 1 - alpha)
      else {
	if(length(y)<=1000)
	  q <- quantile(stepR::MRC.1000, 1 - alpha)
	else
	  q <- quantile(stepR::MRC.asymptotic, 1 - alpha)
      }
    } else {
      q <- MRC.quant(1 - alpha, length(y), r, lengths, "sqrt")
    }
  } else {
    if(missing(r)) {
      if(dyadic)
	alpha <- 1 - ecdf(stepR::MRC.asymptotic.dyadic)(q)
      else {
	if(length(y)<=1000)
          alpha <- 1 - ecdf(stepR::MRC.1000)(q)
	else
          alpha <- 1 - ecdf(stepR::MRC.asymptotic)(q)
      }
    } else {
      alpha <- 1 - MRC.pvalue(q, length(y), r, lengths, "sqrt")
    }
  }
  
  if(family=="gauss")
  {
    if(missing(param))
      param <- sdrobnorm(y)
  }
  
  b <- bounds.MRC(y, q = q, penalty = "sqrt", lengths = lengths, family = family, param = param)
  sb <- stepbound(y, b, x = x, x0 = x0, family = family, param = param, jumpint = jumpint, confband = confband)

  attr(sb, "alpha") <- alpha
  attr(sb, "q") <- q
  sb
}


# thresh.smuceR <- function(v) NA
delayedAssign("thresh.smuceR",
with(.thresh.smuceR.data, approxfun(lengths, q_opt, yleft = min(q_opt), yright = max(q_opt))))
