"fun.bimodal.fit.ml" <-
function(data, first.fit, second.fit, prop, param1, param2, selc1, selc2)
{
index1 <- switch(selc1,
rs = 1,
fmkl = 2,
star = 3)
index2 <- switch(selc2,
rs = 1,
fmkl = 2,
star = 3)
result <- optim(c(first.fit[, index1], second.fit[, index2], prop), optim.fun5, data = data, param1 = param1, param2 = param2, control = list(maxit = 5000))
result <- optim(c(result$par), optim.fun5, data = data, param1 = param1, param2 = param2, control = list(maxit = 5000))
return(result)
}

