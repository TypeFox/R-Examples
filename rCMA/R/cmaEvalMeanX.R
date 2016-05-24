
######################################################################################
# cmaEvalMeanX
#
#' Evaluate the meanX of the current population.
#'
#' After executing \code{\link{cmaOptimDP}}, there is a current population and a best-ever solution.
#' Evaluate for the mean of the current population whether it is feasible and whether
#' the mean is an even better solution. If so, update the best-ever solution.
#'
#' The code of this function is also instructive as a full example for the extensibility of the
#' \href{http://cran.r-project.org/web/packages/rJava/index.html}{\code{rJava}} 
#' interface to CMA-ES. See the full code in \code{demo/demoEvalMeanX}. Some example \code{rJava}-calls are: 
#' \preformatted{
#'    rJava::.jcall(cma,"[D","getMeanX");                 
#'    bestSolutionObj = 
#'      rJava::.jcall(cma,"Lfr/inria/optimization/cmaes/CMASolution;","setFitnessOfMeanX",fitFunc(meanX));
#'    rJava::.jcall(bestSolutionObj,"J","getEvaluationNumber");
#' }
#' Every direct method of classes in the CMA-ES Java package \code{cmaes} (see [Hansen09] for the complete Javadoc
#' and [Hansen13] for an overview on CMA-ES in total) can be accessed with the \code{.jcall}-mechanism
#' of the \href{http://cran.r-project.org/web/packages/rJava/index.html}{\code{rJava}} R package:
#' \preformatted{
#'      rJava::.jcall(obj,returnType,method,...)
#' }
#' where \code{...} stands for the calling parameter(s) of \code{method}. \cr
#' \code{returnType} is a string following the JNI type convention (see, e.g. [Oracle14])
#'    \tabular{ccc}{
#'      \tab Field Descriptor \tab  Java Language Type  \cr
#'      \tab       Z          \tab  boolean  \cr
#'      \tab       C          \tab  char  \cr
#'      \tab       I          \tab  int  \cr
#'      \tab       J          \tab  long  \cr
#'      \tab       F          \tab  float  \cr
#'      \tab       D          \tab  double  \cr
#'      \tab      [I          \tab  int[]  \cr
#'      \tab     [[D          \tab  double[][]  \cr
#'      \tab Ljava/langString;\tab  java.lang.String  \cr
#'      \tab       S          \tab  java.lang.String  \cr
#'      \tab       T          \tab  short  \cr
#'    }
#' (Note: (a) the terminating \code{";"} in \code{"Ljava/langString;"} (!) and (b) \code{"S"} is a short hand for \code{"Ljava/langString;"} and 
#'  \code{"T"} is the re-mapped code for \code{short}. ) \cr\cr
#' The calling parameters in \code{...} have to be matched exactly. In R, numeric vectors are stored as \code{doubles}, so the calling syntax
#' \preformatted{
#'           bestSolutionObj = .jcall(cma,rType,"setFitnessOfMeanX",fitFunc(meanX));
#' }
#' is just right for the Java method \code{setFitnessOfMeanX(double[]) }. In other cases, the calling R variable \code{x}
#' has to be cast explicitly:
#'    \tabular{ccc}{
#'      \tab     Cast         \tab  Java Language Type  \cr
#'      \tab   .jbyte(x)      \tab  byte  \cr
#'      \tab   .jchar(x)      \tab  char  \cr
#'      \tab  as.integer(x)   \tab  int  \cr
#'      \tab   .jlong(x)      \tab  long  \cr
#'      \tab   .jfloat(x)     \tab  float  \cr
#'    }
#'
#' @param cma CMA-ES Java object, already initialized with \code{\link{cmaInit}}
#' @param fitFunc a function to be minimized. Signature: accepts a vector x, returns a double.
#' @param isFeasible [\code{function(x)\{TRUE\}}] a Boolean function checking the feasibility
#'    of the vector \code{x}. The default is to return always \code{TRUE}.
#'
#' @return \code{bestSolution}, a list with entries:
#'       \item{\code{bestX}}{  a vector of length \code{dimension} containing the best-ever 
#'          solution, including \code{meanX}  }
#'       \item{\code{meanX}}{  a vector of length \code{dimension} containing the mean of 
#'          the current (last) population in \code{cma} }
#'       \item{\code{bestFitness}}{  the best-ever fitness value, including the evaluation of meanX  }
#'       \item{\code{bestEvalNum}}{  the function evaluation count where \code{bestFitness} occured  }
#'       \item{\code{lastEvalNum}}{  the total function evaluation count. If \code{bestEvalNum==lastEvalNum}
#'          then the best-ever fitness occured in the evaluation of \code{meanX}.  }
#'
#' @references
#'  [Hansen09] \url{https://www.lri.fr/~hansen/javadoc} Nikolaus Hansen: Javadoc for CMA-ES Java package fr.inria.optimization.cmaes, 2009.  \cr
#'  [Hansen13] \url{https://www.lri.fr/~hansen/cmaesintro.html} Nikolaus Hansen: The CMA Evolution Strategy, 2013. \cr
#'  [Oracle14] \url{http://docs.oracle.com/javase/7/docs/technotes/guides/jni/spec/jniTOC.html} 
#'    Oracle: The Java Native Interface. Programmer's Guide and Specification. 
#'    Chapter 3 (JNI types), Sec. 'Type Signatures', 2014.
#'
#' @examples
#'    \dontrun{ 
#'    ## just to show the syntax, without calling cmaOptimDP
#'    fitFunc <- function(x) {  sum(x*x); }
#'    isFeasible <- function(x) {  TRUE;  }
#'    cma <- cmaNew(propFile="CMAEvolutionStrategy.properties");
#'    cmaInit(cma,dimension=2,initialX=1.5);
#'    bestSolution=cmaEvalMeanX(cma,fitFunc,isFeasible);
#'    str(bestSolution);
#'    }
#'
#' @seealso   \code{\link{cmaInit}}, \code{\link{cmaOptimDP}}
#' @author Wolfgang Konen, FHK, 2013-2015
#' @export
######################################################################################
cmaEvalMeanX <- function(cma,fitFunc,isFeasible=function(x)TRUE) {
    meanX = rJava::.jcall(cma,"[D","getMeanX");
    rType = "Lfr/inria/optimization/cmaes/CMASolution;";
    if (isFeasible(meanX)) {
      bestSolutionObj = rJava::.jcall(cma,rType,"setFitnessOfMeanX",fitFunc(meanX));
      # setFitnessOfMeanX compares fitFunc(meanX) with the current best-ever solution.
      # If fitFunc(meanX) is better, it updates the best-ever solution.
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
      bestSolution$bestX = bestSolution$meanX
    
    return(bestSolution);
}

