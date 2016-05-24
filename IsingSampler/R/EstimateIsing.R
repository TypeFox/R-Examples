# Wrapper function for:
# - pl: pseudolikelihood
# - uni: univariate logistic regressions
# - bi: bivariate logistic regressions
# - ll: Loglinear model

EstimateIsing <- function(data, responses, beta = 1, method = c('pl', 'uni', 'bi', 'll'),...){

  switch(method[[1]],
        pl = EstimateIsingPL(data, responses, beta, ...)  ,
        uni = EstimateIsingUni(data, responses, beta, ...)  ,
        bi = EstimateIsingBi(data, responses, beta, ...),
        ll = EstimateIsingLL(data, responses, beta, ...)
         )
}