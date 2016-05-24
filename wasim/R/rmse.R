rmse <-
function (a, b, cond = rep(TRUE, NROW(a))){

cond[is.na(a)] <- FALSE

cond[is.na(b)] <- FALSE

a <- a[cond]

b <- b[cond]

return(sqrt(sum((a-b)^2)/NROW(a)))

}

