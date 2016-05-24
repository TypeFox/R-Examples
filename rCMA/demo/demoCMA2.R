## demo/demoCMA2.R: demo usage of package rCMA.
##
## After doing the unconstrained sphere (as in demoCMA1.r, for later reference in plot), 
## the constrained sphere problem TR2 is solved. 
## The problem TR2 has in addition to the fitness function 'sphere' the constraint function
## 'above the hyperplane sum_i(x_i) = n', where n is the input space dimension.
## The constrained optimum is at (1,...,1) and it has the value fTarget2=n.
##
## Note that in this second case, the optimimum lies exactly at the boundary 
## of the feasible region: res2$bestX=c(1,...,1).
##
## This script does exactly the same as class CMAExampleConstr in cma_jAll.jar,
## but it allows to define the functions fitFunc and isFeasible on the R side. 
## They can be replaced by arbitrary other R functions, which may depend on other 
## R variables as well. 
##
## The constraint handling approach is a very simple one: death penalty, i.e. if we get an 
## infeasible individual, it is immediately discarded and a new one is drawn from the distribution. 
## (This approach will run into trouble if the current distribution does not allow to reach any  
## feasible solutions.)
## 
library(rCMA)
fitFunc <- function(x) {  sum(x*x); }
isFeasible <- function(x) {  (sum(x) - length(x)) >= 0;  }
n = 2;

cma <- cmaNew(propFile="CMAEvolutionStrategy.properties");
cmaInit(cma,seed=42,dimension=n,initialX=1.5, initialStandardDeviations=0.2);
res1 = cmaOptimDP(cma,fitFunc,iterPrint=30);

cma <- cmaNew(propFile="CMAEvolutionStrategy.properties");
cmaInit(cma,seed=42,dimension=n,initialX=1.5, initialStandardDeviations=0.2);
res2 = cmaOptimDP(cma,fitFunc,isFeasible,iterPrint=30);

fTarget =c(0,n); 
plot(res1$fitnessVec-fTarget[1],type="l",log="y",xlim=c(1,max(res1$nIter,res2$nIter))
    ,xlab="Iteration",ylab="Distance to target fitness");
lines(res2$fitnessVec-fTarget[2],col="red");
legend("topright",legend=c("TR2","sphere"),lwd=rep(1,2),col=c("red","black"))
str(res2);

bestSolution=rCMA::cmaEvalMeanX(cma,fitFunc,isFeasible);
str(bestSolution);

