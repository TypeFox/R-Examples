logit <- function(p){
#   logit of p
    log(p/(1-p))
}

expit <- function(theta){
# anti logit function
   1/(1+exp(-theta))
}
