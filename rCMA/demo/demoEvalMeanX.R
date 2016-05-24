## This demo just displays the complete code of rCMA::cmaEvalMeanX
##
## 
cmaEvalMeanX <- function(cma,fitFunc,isFeasible) {
    meanX = rJava::.jcall(cma,"[D","getMeanX");
    rType = "Lfr/inria/optimization/cmaes/CMASolution;";
    if (isFeasible(meanX)) {
      bestSolutionObj = rJava::.jcall(cma,rType,"setFitnessOfMeanX",fitFunc(meanX));
      # setFitnessOfMeanX may update the best-ever solution
    }  else {
      bestSolutionObj = rJava::.jcall(cma,rType,"getBestSolution");
    }
    bestSolution = list(bestX = rJava::.jcall(bestSolutionObj,"[D","getX")
                      , meanX = meanX
                      , bestFitness = rJava::.jcall(bestSolutionObj,"D","getFitness")
                      , bestEvalNum = rJava::.jcall(bestSolutionObj,"J","getEvaluationNumber")
                      , lastEvalNum = rJava::.jcall(cma,"J","getCountEval")
                      );
    if (bestSolution$bestEvalNum == bestSolution$lastEvalNum & isFeasible(meanX))
      bestSolution$bestX = bestSolution$meanX;

    return(bestSolution);
}

## just to show the syntax, without calling cmaOptimDP
fitFunc <- function(x) {  sum(x*x); }
isFeasible <- function(x) {  TRUE;  }
cma <- cmaNew();
cmaInit(cma,dimension=2,initialX=1.5);
bestSolution=cmaEvalMeanX(cma,fitFunc,isFeasible);
str(bestSolution);
