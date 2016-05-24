## demo usage of package rCMA.
##
## The problem solved is the (unconstrained) sphere problem with dimension 2. 
## 
fitFunc <- function(x) {  sum(x*x); }

cma <- cmaNew();
cmaInit(cma,seed=42,dimension=2,initialX=1.5, initialStandardDeviations=0.2);
res1 = cmaOptimDP(cma,fitFunc,iterPrint=30);

plot(res1$fitnessVec,type="l",log="y",col="blue",xlab="Iteration",ylab="Fitness");
str(res1);



