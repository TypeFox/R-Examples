RunTime <- function(time, U, prob=seq(0,1,.1)) {
#
#  Take the posterior sample of U[1,...nstrata] and compute the percentiles of the run
#  timing. 
#  This uses the quantile() function from the "actuar" package which is designed to compute
#  quantiles of grouped data.
#
#  Input arguments:
#     time - list of sample weeks. It is assumed that there are no salmon prior to the first value
#            in time, and after the last value in time
#     U  - matrix of posterior samples. Each row is a sample from the posterior.
#          Columns correspond to U[1]...U[nstrata]
# 

 timing <- c(min(time):(1+max(time)))
 q.U <- NULL
 #browser()
 for( i in 1:nrow(U)) {  # go through each sample from the posterior    
    U.sample <- U[i,]
    quant <- quantile(grouped.data(Group=timing, Frequency=U.sample), prob=prob)
    q.U <- rbind(q.U, quant)
 }
 q.U

}
