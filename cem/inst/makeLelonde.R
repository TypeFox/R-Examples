require(cem)

data(LL)

ans <- c("strongly disagree", "disagree", "neutral", "agree", "strongly agree", "no opinion")

n <- dim(LL)[1]
k <- dim(LL)[2]

set.seed(123)

LeLonde <- data.frame(LL, q1 = sample(ans, n, replace=TRUE) )

idx <- sample(1:n,n*0.1)
sapply(idx, function(i) { LeLonde[i, sample(2:(k+1),1)] <<- NA} )
 

save(LeLonde, file="LeLonde.rda")
