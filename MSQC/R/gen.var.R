gen.var <-
function(x,...) {  # Generalized variance
# Note: must be n > p
#
#
p <- ncol(x) # quality characteristics
m <- nrow(x) # sample
n <- dim(x)[3] # observations or sample size

b1 <- 1 / (n - 1) ^ p * cumprod(n - 1:p)[length(cumprod(n - 1:p))]
b2 <- 1 / (n - 1) ^ (2 * p) * cumprod(n - 1:p)[length(cumprod(n - 1:p))] * 
(cumprod(n - 1:p + 2)[length(cumprod(n - 1:p + 2))] - 
cumprod(n - 1:p)[length(cumprod(n - 1:p))])

S <- covariance(x)
stat <- covariance(x,stat)

LCL <- max(0,(det(S) / b1 * (b1 - 3 * b2 ^ .5))) 
CL <- det(S)
UCL <- det(S) / b1 * (b1 + 3 * b2 ^ .5) 


t3 <- which(stat > UCL || stat < LCL )

if(any(stat > UCL || stat < LCL)){
     cat("The following(s) point(s) fall outside of the control limits" )
      print(t3)}
 
par(mar=c(4,5,3,5))
plot(stat,ylim = c(0,1.1 * max(max(stat),UCL)), main="Generalized Variance Control Chart",
 xlab = "Observation",ylab = expression(det(S)),type = "o") 
mtext(paste(" UCL=", round(UCL, 2)),side=4, at=UCL,las=2)
mtext(paste(" CL=", round(CL, 2)),side=4, at=CL,las=2)
mtext(paste(" LCL=", round(LCL, 2)),side=4, at=LCL,las=2)

points(t3,stat[t3],col = 2) 
segments(0, UCL, m, UCL, col = 2)
segments(0, LCL, m, LCL, col = 2)
segments(0, CL, m, CL, col = 3)

outList = list ("Generalized Variance Control Chart","Upper Control Limit" =  signif(UCL,2),
"Lower Control Limit" =  signif(LCL,2), stat =  signif(stat,2))
 return(outList)

}
