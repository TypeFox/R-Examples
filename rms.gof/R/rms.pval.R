rms.pval <- function(observed, expected, num_sim=1000){

T <- num_sim
if(!is.vector(observed)){
observed <- as.vector(observed)
}
o <- length(observed)
e <- length(expected)
if(o != e){
return(cat("ERROR: Observed and Expected vectors must be same size"))
}

rms <- 0
rarities <- 0
sum_obs <- sum(observed)
sum_exp <- sum(expected)

X <- test.rms(observed,expected)

prob_exp <- expected/sum_exp

for(j in 1:T){
cdvec <- cumsum(prob_exp)
samp_exp <- array(0,c(1,o))
for (i in 1:sum_exp){
t <- runif(1)
ind <- which(cdvec > t)
first <- ind[1]
samp_exp[1,first] <- samp_exp[1,first] + 1
}

rms <- test.rms(observed,samp_exp)
if(rms >X){
rarities <- rarities + 1
}
}

pval <- rarities/T

return(pval)
}
