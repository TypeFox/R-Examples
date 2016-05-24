invlink <- function(x, type = c("cloglog", "log", "logit", "probit")){

  switch(type,
    cloglog = 1 - exp(-exp(x)),
    log = exp(x),
    logit = plogis(x),
	  probit = pnorm(x)
  )

}
