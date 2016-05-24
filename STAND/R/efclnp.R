efclnp <-
function(dd,gam=0.95,L){
#  For random sample size n from any distribution
#  exceedance fraction F= Pr[ x > L]  See Section 3.4
#  (FL,FU) is 100gam-percent CI for F
# USAGE: efclnp(dd=aiha,gam=0.95,L=5)
# ARGUMENTS: dd = matrix dd with x[i] in column 1 and det[i] in col 2
#           gam= confidence level one sided
#           and L= Limit Fo null value(%)
# VALUE: estimate f of F and exact 100*[ 1 - 2*gam ] percent
#       Confidence Interval for the exceedance fraction F= Pr[ x > L]
# DETAILS: see R function  binom.test()
# ASSUMPTION: all non-detects are less than L
nx<- sum( ifelse(dd[,1] > L,1,0) )
n<-dim(dd)[1]
ef<- nx/n
# f0 <- 1 - p
clev<- 1 - (1 - gam)*2
tmp<- binom.test(nx,n,ef,"two.sided",clev)
fcl<- 100*as.numeric( unlist(tmp[4]) )

out<-list("fnp"=100*ef,"fnp.LCL"=fcl[1],"fnp.UCL"=fcl[2],"L"=L,"gamma"=gam)
out

}

