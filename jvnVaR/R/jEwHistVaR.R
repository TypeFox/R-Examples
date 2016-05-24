jEwHistVaR <-
function(s, alpha){
# critizim: each method, this time, hard to choose h
# propose 1 - h > 1/nn, or h < 1- 1/nn
# h shouldn't too large, cause 1 - h too small
# max_cdf = 1 - h^nn shouldn't too small
# propose 1 - h^nn > 1 - alpha
# #or h < alpha ^(1/nn)
# h shouldn't too small, cause (1 - h^(nn-k)) * h^k too small
# 0.95 ^ 100 ~ 0.006; 0.99 ^ 250 ~ 0.08; 0.05^(1/100) ~ 0.97
# basic way: alpha = 0.01 -> smallest in 100.
# #propose h ^ nn > alpha/2
nn <- length(s)
h <- 1-1/nn
x <- sort(s, index.return = TRUE)
h <- min(max(h,(alpha/2)^(1/nn)),1-1/nn, alpha^(1/nn))
weight_proba <- c(1-h, h*rep(1,nn-1))
weight_proba <- cumprod(weight_proba) #weight probability reflect each return
EXW_pmf <- weight_proba[x$ix]
EXW_cdf <- cumsum(EXW_pmf) #cdf of return following exponent weight
if (alpha >= EXW_cdf[1]){ #find quantile
    if (alpha >= EXW_cdf[nn]){
        id <- nn
invF <- x$x[id]
    }else {   
a <- 1:nn
y <- a[EXW_cdf > alpha]
id <- y[1]
temp <- (alpha - EXW_cdf[id-1])/(EXW_cdf[id] - EXW_cdf[id-1])
invF <- x$x[id] * temp + x$x[id-1] * (1 - temp)
    }
}else {
    id <- 1
invF <- x$x[id]
}
return(-invF)
}
