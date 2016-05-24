## Transformation of a variable from the observation scale to the link scale
link <- function(x, type = c("cloglog", "log", "logit")){
  switch(type,
    logit = qlogis(x),
    log = log(x),
    cloglog = log(-log(1 - x)))
    }

## Transformation of a variable from the link scale to the observation scale
invlink <- function(x, type = c("cloglog", "log", "logit")){
  switch(type,
    logit = plogis(x),
    log = exp(x),
    cloglog = 1 - exp( -exp(x)))
    }

