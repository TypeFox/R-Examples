link <- function(x, type = c("cloglog", "log", "logit", "probit")){
  
  switch(type,
    cloglog = log(-log(1 - x)),
    log = log(x),
    logit = qlogis(x),
	  probit = qnorm(x)
  )
  
}
