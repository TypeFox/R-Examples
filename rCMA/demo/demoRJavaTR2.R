################################################################################ 
## Script to show the usage of CMA-ES (Java) from R (based on package rJava)
##
## The code does virtually the same as the TR2-part from demo demoCMA2.r, but 
## it does not use any functions from rCMA. Instead it uses only functions of 
## the underlying package rJava. It may thus serve for the developer as a template
## how to access methods in Java class CMAEvolutionStrategy which are currently
## not in the interface of rCMA. See the documentation page of cmaEvalMeanX.R  
## for further details.
##
## The problem solved is TR2 (tangent problem with dim=2, fitFunc=sphere with a
## linear tangent constraint in function isFeasible, see demoCMA2.R). 
##
require(rJava)
.jinit()
cmaJarFile <- paste(find.package("rCMA"), "java/cma_jAll.jar",sep="/");
if (!file.exists(cmaJarFile)) stop(sprintf("cmaJarFile=%s not found",cmaJarFile));
.jaddClassPath(cmaJarFile);
propFile <- paste(find.package("rCMA"),"CMAEvolutionStrategy.properties",sep="/")


fitFunc <- function(x) {
    sum(x*x);
}
isFeasible <- function(x) {
    (sum(x) - length(x)) > 0;
}

cma <- .jnew("fr/inria/optimization/cmaes/CMAEvolutionStrategy")
props <- .jcall(cma,"Ljava/util/Properties;","readProperties",propFile);
## function readProperties returns an object of class Properties --> the correct
## JNI specification of the return type is "Ljava/util/Properties;" 
## (Note the terminating ";" !). See cmaEvalMeanX.R
## To assess elements of props:
##       initialX = .jcall(props,"S","getProperty","initialX");

.jcall(cma,,"setSeed",.jlong(42));      #  only effective if called before cma.init()
fitness <- .jcall(cma,"[D","init",as.integer(2),.jarray(c(1.5)),.jarray(c(0.2)));
## initialize CMA-ES object:      dimension    , initialX      , initialStdDev  ;
## this is a  short version for:
##    .jcall(cma,,"setDimension",as.integer(2));
##    .jcall(cma,,"setInitialX",.jarray(c(1.5));)
##    .jcall(cma,,"setInitialStandardDeviations",.jarray(c(0.2)));
##    fitness <- .jcall(cma,"[D","init");

cfe = 0;          # number of constraint function evaluations (isFeasible)
fitnessVec=NULL;  # best-ever fitness after each iteration
while(.jcall(cma,"I","getStopConditionNumber")==0) {
    pop = .jcall(cma,"[[D","samplePopulation");   # default population size is 6
    LP = length(pop);       # population size
    popR = sapply(pop,.jevalArray)   
    ## popR is a (2 x LP) R-matrix  with popR[,1] being the first individuum in the population
    cfe = cfe + LP;
    
    resampleFunc <- function(i) {
      while (!isFeasible(popR[,i])) {
        ## death penalty (DP): if an individual is not feasible, throw it away and sample a new one:
        popR[,i] <<- .jcall(cma,"[D","resampleSingle",as.integer(i-1));
        ## it is important to use "<<-" here, otherwise a local copy popR[,i] is 
        ## made and nothing changed in the popR of the while-loop (!)
        cfe <<- cfe+1;
      }
    }
    feasible = sapply(1:LP, resampleFunc);
    fitness = sapply(1:LP, function(i) { x=popR[,i]; fitFunc(x); });
    iter = .jcall(cma,"J","getCountIter")    # note that "J" is the JNI field descriptor for 'long'
    if (iter %% 30 ==0) {
      i = which.min(fitness);    
      dimPrint = min(5,.jcall(cma,"I","getDimension"));
      cat(sprintf("%04d  %e | %s\n", iter, fitness[i], paste(sprintf("%11.4e",popR[1:dimPrint,i]),collapse=", ")));
      flush.console()
    }
    .jcall(cma,,"updateDistribution",fitness);
    fitnessVec[iter] = .jcall(cma,"D","getBestFunctionValue");
}  
ffe = .jcall(cma,"J","getCountEval")    # number of fitness function evaluations (fitFunc)
scMsg = .jcall(cma,"[S","getStopConditionMessages");
cat(sprintf("Terminated due to %s\n",scMsg));
cat(sprintf("cfe,ffe, %%infeasible: %d %d %f\n",cfe,ffe,(cfe-ffe)/ffe));

plot(fitnessVec-2,type="l",log="y",xlab="Iteration",ylab="Distance to target fitness");
