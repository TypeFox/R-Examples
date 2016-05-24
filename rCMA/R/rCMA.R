######################################################################################
######################################################################################
# Package Description for Roxygene:
#' CMA-ES R-to-Java interface
#'
#' \tabular{ll}{
#' Package: \tab rCMA\cr
#' Type: \tab Package\cr
#' Version: \tab 1.1\cr
#' Date: \tab 2015-04-30\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' rCMA is a package to perform CMA-ES optimization, using the *Java* implementation by Niko Hansen [Hansen2009]. \cr
#' CMA-ES [HansOst96, Hansen13] is the Covariance Matrix Adapting Evolutionary Strategy for numeric black box optimization. \cr 
#' The main features of rCMA are: \enumerate{
#'  \item Abiltiy to start the Java CMA-ES optimization with fitness functions defined in R. 
#'  \item Constraint handling: Arbitrary constraints can be incorporated, see function parameter \code{isFeasible} 
#'    in \code{\link{cmaOptimDP}}. 
#'  \item Extensibility: Full access to all methods of the Java class \code{CMAEvolutionStrategy} through package  
#'    \href{http://cran.r-project.org/web/packages/rJava/index.html}{\code{rJava}}. New methods can be added easily.
#'    See the documentation of \code{\link{cmaEvalMeanX}} for further details, explanation of JNI types and a full example. 
#'  \item Test and Debug: The access of Java methods from R allows for easy debugging and test
#'    of programs using \code{CMAEvolutionStrategy} through R scripts without the necessity to
#'    change the underlying JAR file.
#' }
#' The main entry point functions are \code{\link{cmaNew}}, \code{\link{cmaInit}} and \code{\link{cmaOptimDP}}. 
#'
#' Note: To install \code{rJava} properly on some Unix systmes, it might be necessary to issue as
#' root the command \code{R CMD javareconf} once, or, as normal user to issue the command \code{R CMD javareconf -e}
#' prior to installing package \code{rJava} or prior to loading library \code{rJava}. 
#' 
#' @name rCMA-package
#' @aliases rCMA
#' @docType package
#' @title R interface to the Java CMA-ES of Niko Hansen
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de})
#' @references   
#'  [HansOst96] Hansen, N. and Ostermeier, A.: Adapting arbitrary normal mutation distributions in evolution strategies: The covariance matrix adaptation. 
#'    In Proceedings of the 1996 IEEE International Conference on Evolutionary Computation, pp. 312-317, 1996. 
#'    \href{https://www.lri.fr/~hansen/CMAES.pdf}{PDF}. \cr
#'  [Hansen09] \url{https://www.lri.fr/~hansen/javadoc} Nikolaus Hansen: Javadoc for CMA-ES Java package fr.inria.optimization.cmaes, 2009. \cr
#'  [Hansen13] \url{https://www.lri.fr/~hansen/cmaesintro.html} Nikolaus Hansen: The CMA Evolution Strategy Web Page, 2013. \cr
#'  [Urbanek13] \url{http://cran.r-project.org/web/packages/rJava}
#'    Urbanek, S.: rJava: Low-level R to Java interface, 2013. \cr
#'  [Oracle14] \url{http://docs.oracle.com/javase/7/docs/technotes/guides/jni/spec/jniTOC.html} 
#'    Oracle: The Java Native Interface. Programmer's Guide and Specification, 2014.
#' @keywords package CMA Covariance Matrix rJava
#End of Package Description
NA #NULL, ends description without hiding first function
######################################################################################
######################################################################################

## When package rCMA is loaded (e.g. via library(rCMA)), then this function .onLoad is 
## called with libname=<path-to-your-R-library> and pkgname="rCMA". The function 
##    rJava::.jpackage
## will attach all JAR-files found in the java/-subdir of rCMA to the Java class path
## (see .jclassPath() )
#----
#.onLoad <- function(libname, pkgname) {
#  rJava::.jpackage(pkgname, lib.loc=libname)
#}
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("Attaching %s/%s/java/cma_jAll.jar to the Java class path (.jclassPath())",libname,pkgname));
}

######################################################################################
# cmaNew
#
#' Create a new CMA-ES Java object.
#'
#' @param propFile [NULL] filename of a file with property settings. If NULL, read file \cr
#'    \code{CMAEvolutionStrategy.properties} from the package directory (\code{find.package("rCMA")})
#'
#' @return the new CMA-ES Java object of class \code{CMAEvolutionStrategy}, which has as
#'    additional attribute \code{props}, the Java \code{Properties} object as read from \code{propFile}.
#'
#' @note 
#'    The default properties file can be found in \code{CMAEvolutionStrategy.properties}.  
#'    A read-only copy can be inspected by browsing to "Index" (of package rCMA), then "Overview of user guides ...".  \cr
#'    It allows to set more parameter, especially more  \link[=cmaSetStopFitness]{stop conditions}. 
#'
#' @examples
#'    ## show how element initialX can be inferred from attribute props:
#'    ## (see  cmaEvalMeanX-documentation for further details on .jcall and its argument "S")
#'    cma <- cmaNew();
#'    props <- attr(cma,"props");
#'    initialX = rJava::.jcall(props,"S","getProperty","initialX");
#'    print(initialX);
#'
#' @seealso   \code{\link{cmaInit}}
#' @author Wolfgang Konen, FHK, 2013
#' @export
######################################################################################
cmaNew <- function(propFile=NULL) {
  #require(rJava)
  cmaJarFile <- paste(find.package("rCMA"), "java/cma_jAll.jar",sep="/");
  #cmaJarFile = "../inst/java/cma_jAll.jar"     # only for developer version
  if (!file.exists(cmaJarFile)) stop(sprintf("cmaJarFile=%s not found in working directory",cmaJarFile));
  if (is.null(propFile)) propFile <- paste(find.package("rCMA"),"CMAEvolutionStrategy.properties",sep="/")
  ##
  ## This is now necessary again, because the  .onLoad mechanism above does not work with R3.1.3
  rJava::.jinit()          # initialize Java-VM
  rJava::.jaddClassPath(cmaJarFile);
  ##
  ## For unclear reasons, the following two lines gave an error in system("R CMD check rCMA") (and only there!) 
  ## when running any example containing cmaNew(). The error was not observed when executing cmaNew() from the 
  ## R prompt directly. Anyhow, since we do not need them really, we commented the two following lines out:
  #if (length(grep("cma_jAll.jar",.jclassPath()))==0)
  #  warning("cma_jAll.jar is not part of .jclassPath(). Check whether rJava::.jaddClassPath has access to the right directory.");
  cma <- rJava::.jnew("fr/inria/optimization/cmaes/CMAEvolutionStrategy")
  props <- rJava::.jcall(cma,"Ljava/util/Properties;","readProperties",propFile);
  ## To assess elements of props:
  ##       initialX = rJava::.jcall(props,"S","getProperty","initialX");
  attr(cma,"props") <- props
  cma;
}

######################################################################################
# cmaInit
#
#' Initialize a CMA-ES Java object.
#'
#'
#' @param cma CMA-ES Java object, as created by \code{\link{cmaNew}}
#' @param seed [NULL] if not NULL, set the seed to the given value
#' @param dimension [NULL] if not NULL, overwrite the dimension setting from \code{propFile} (\code{\link{cmaNew}})
#' @param initialX [NULL] if not NULL, overwrite the initialX setting from \code{propFile} (\code{\link{cmaNew}}).
#'    initialX can be a double or a double vector of length \code{dimension}.
#' @param initialStandardDeviations [NULL] if not NULL, overwrite the initialStandardDeviations
#'    setting from \code{propFile} \code{\link{cmaNew}}. initialStandardDeviations can be a double or a double 
#'    vector of length \code{dimension}.
#'
#' @return \code{fitness}, a vector of 0's with the length of the intended population.
#'
#'
#' @note
#'  As a side effect, the CMA-ES Java object \code{cma} of class \code{CMAEvolutionStrategy}
#'  is transferred into an augmented state. As a second side effect, the population size is 
#'  set to 
#'        \deqn{\lambda = 4 + 3 floor(ln(n))}
#'  where \eqn{n=}\code{dimension}.
#'
#' @examples
#'    cma <- cmaNew();
#'    cmaInit(cma,seed=42,dimension=2,initialX=1.5);
#'
#' @seealso   \code{\link{cmaNew}}, \code{\link{cmaOptimDP}}
#' @author Wolfgang Konen, FHK, 2013
#' @export
######################################################################################
cmaInit <- function(cma,seed=NULL,dimension=NULL,initialX=NULL, initialStandardDeviations=NULL) {
  if (!is.null(seed)) rJava::.jcall(cma,,"setSeed",rJava::.jlong(seed));
  if (!is.null(dimension)) rJava::.jcall(cma,,"setDimension",as.integer(dimension));
  if (!is.null(initialX)) rJava::.jcall(cma,,"setInitialX",rJava::.jarray(initialX));
  if (!is.null(initialStandardDeviations)) rJava::.jcall(cma,,"setInitialStandardDeviations",rJava::.jarray(initialStandardDeviations));
  fitness <- rJava::.jcall(cma,"[D","init");
}

######################################################################################
# cmaOptimDP
#
#' Perform a CMA-ES optimization with constraints (DP).
#'
#' The optimization uses DP (death penalty) for handling constraint violations: 
#' Each time an infeasible individual is encountered, it is thrown away  
#' and a new individual is resampled from the CMA distribution.
#'
#' This functions loops through iterations (generations) until a \link[=cmaSetStopFitness]{stop condition} is met: 
#' In each iteration, a population is sampled (infeasible individuals are replaced via
#' Java function \code{resampleSingle}), its fitness vector is evaluated and the CMA 
#' distribution is updated according to this fitness vector.
#' 
#' Every \code{iterPrint} generations a one-line diagnostic output of the form
#' \preformatted{
#'      iter  fitness | x1 x2 ... xp
#' } 
#' is printed where \code{fitness} is the current best value of the fitness function to be minimized
#' and \code{x1 x2 ... xp} are the first \code{maxDimPrint} dimensions of the corresponding 
#' best point in input space.
#' 
#' @param cma CMA-ES Java object, already initialized with \code{\link{cmaInit}}
#' @param fitFunc a function to be minimized. Signature: accepts a vector \code{x}, returns a \code{double}.
#' @param isFeasible [\code{function(x)\{TRUE\}}] a Boolean function checking the feasibility
#'    of the vector \code{x}. The default is to return always \code{TRUE}.
#' @param maxDimPrint [5] how many dimensions of vector \code{x} to print in diagnostic output 
#' @param iterPrint [10]  after how many iterations should diagnostic output be printed?
#' @param verbose [2] possible values are 0 (no output), 1, 2 (much output)
#'
#' @return \code{res}, a list with diagnostic output from the optimization run:
#'       \item{\code{sMsg}}{  a string vector with all console output from the optimization run. 
#'            To print it, use: \code{ cat(sMsg) } or \code{ for (x in sMsg) cat(x) }  }
#'       \item{\code{bestX}}{ vector of length \code{dimension} with the best-ever solution X }
#'       \item{\code{bestFitness}}{  the corresponding best-ever fitness  }
#'       \item{\code{bestEvalNum}}{ the fitness function evaluation number which gave this best-ever result }
#'       \item{\code{nIter}}{ number of iterations  }
#'       \item{\code{fitnessVec}}{ vector of length \code{nIter}: the best fitness after each iteration }
#'       \item{\code{xMat}}{  (\code{nIter x dimension})-matrix: \code{xMat[i,]} is the best solution X after iteration \code{i} }
#'       \item{\code{cfe}}{  number of constraint function evaluations (\code{isFeasible}) }
#'       \item{\code{ffe}}{  number of fitness function evaluations (\code{fitFunc}) }
#'
#' @note
#'    If your fitness function depends on other parameters besides \code{x}, then 
#'    encapsulate it in a new function \code{fitFunc} at a place where the other parameters
#'    are accessible and rely on R's mechanism to locate the other parameters 
#'    in the environment surrounding \code{fitFunc}:
#'    \preformatted{
#'      par1 <- someObject;
#'    } \preformatted{
#'      fitFunc <- function(x) {  myFuncWithOtherPars(x,par1); }
#'    }  
#'
#' @example    demo/demoCMA2.R
#'
#' @seealso   \code{\link{cmaNew}}, \code{\link{cmaInit}}
#' @author Wolfgang Konen, FHK, 2013-2015
#' @export
######################################################################################
cmaOptimDP <- function( cma
                      , fitFunc
                      , isFeasible=function(x){TRUE}
                      , maxDimPrint=5, iterPrint=10, verbose=2)
{
  if (rJava::.jcall(cma,"D","getState")<0) stop("CMA object is not yet initialized. Use cmaInit() first.");
  cfe = 0;
  sMsg = NULL;      # holds console output as string vector
  fitnessVec=NULL;  # best-ever fitness after each iteration
  xMat = NULL;      # (nIter x dimension) matrix: row i holds bestX after iteration i
  iter = rJava::.jcall(cma,"J","getCountIter");          # when the while-loop is never entered 
  while(rJava::.jcall(cma,"I","getStopConditionNumber")==0) {
      pop = rJava::.jcall(cma,"[[D","samplePopulation");
      LP = length(pop);       # population size
      popR = sapply(pop,rJava::.jevalArray)
      ## popR is a (2 x LP)-R-matrix  with popR[,1] being the first individuum in the population
      cfe = cfe + LP;

      resampleFunc <- function(i) {
        while (!isFeasible(popR[,i])) {
          ## death penalty (DP): if an individual is not feasible, throw it away and sample a new one:
          popR[,i] <<- rJava::.jcall(cma,"[D","resampleSingle",as.integer(i-1));
          ## it is important to use "<<-" here, otherwise a local copy popR[,i] is
          ## made and nothing changed in the popR of the while-loop (!)
          cfe <<- cfe+1;
        }
      }
      feasible = sapply(1:LP, resampleFunc);
      fitness = sapply(1:LP, function(i) { x=popR[,i]; fitFunc(x); });
      iter = rJava::.jcall(cma,"J","getCountIter")    # note that "J" is the JNI field descriptor for 'long'
      if (iter %% iterPrint == 0) {
        i = which.min(fitness);
        dimPrint = min(maxDimPrint,rJava::.jcall(cma,"I","getDimension"));
        s=sprintf("%04d  %e | %s\n", iter, fitness[i], paste(sprintf("%11.4e",popR[1:dimPrint,i]),collapse=", "));
        sMsg = c(sMsg,s);
        if (verbose>0) {
          cat(s);
          flush.console();
        }
      }
      rJava::.jcall(cma,,"updateDistribution",fitness);
      fitnessVec[iter] = rJava::.jcall(cma,"D","getBestFunctionValue");
      xMat = rbind(xMat,rJava::.jcall(cma,"[D","getBestX"));
  }
  bestSolution = rJava::.jcall(cma,"Lfr/inria/optimization/cmaes/CMASolution;","getBestSolution");
  ffe = rJava::.jcall(cma,"J","getCountEval")    # number of fitness function evaluations (fitFunc)
  s1 = sprintf("Terminated due to %s\n",rJava::.jcall(cma,"[S","getStopConditionMessages"));
  s2 = sprintf("cfe,ffe, %%infeasible: %d %d %f\n",cfe,ffe,(cfe-ffe)/ffe);
  sMsg = c(sMsg,s1,s2);
  if (verbose>0) cat(s1,s2);

  res = list( sMsg=sMsg      # to print it: 'cat(sMsg)' or 'for (x in sMsg) cat(x)'
            , bestX = rJava::.jcall(bestSolution,"[D","getX")
            , bestFitness = rJava::.jcall(bestSolution,"D","getFitness")
            , bestEvalNum = rJava::.jcall(bestSolution,"J","getEvaluationNumber")
            , nIter = iter   # number of iterations
            , fitnessVec = fitnessVec
            , xMat = xMat
            , cfe = cfe      # number of constraint function evaluations (isFeasible)
            , ffe = ffe      # number of fitness function evaluations (fitFunc)
            );
}

######################################################################################
# cmaSamplePopulation
#
#' Sample a population from the current CMA-ES distribution.
#'
#' The population size is given by \code{\link{cmaGetPopulationSize}(cma)}. It can be
#' either set manually with \code{\link{cmaSetPopulationSize}(cma,p)}, prior to 
#' \code{\link{cmaInit}(cma)}, or CMA-ES will use the default population size\cr
#'     \code{popSize = 4 + 3*log(dimension)}.
#'
#' @param cma CMA-ES Java object, already initialized with \code{\link{cmaInit}}
#'
#' @return \code{popR}, a (\code{dimension x popSize}) matrix  with \code{popR[,1]} 
#'    being the first individuum in the population. \cr
#'    \code{dimension = \link{cmaGetDimension}(cma)} \cr
#'    \code{popSize = \link{cmaGetPopulationSize}(cma)} \cr
#'
#'
#' @examples
#'    cma <- cmaNew();
#'    cmaInit(cma,dimension=2,initialX=1.5);
#'    popR <- cmaSamplePopulation(cma);
#'
#' @seealso   \code{\link{cmaUpdateDistribution}}, \code{\link{cmaNew}}
#' @author Wolfgang Konen, FHK, 2013
#' @export
######################################################################################
cmaSamplePopulation <- function(cma) {
      pop = rJava::.jcall(cma,"[[D","samplePopulation");
      popR = sapply(pop,rJava::.jevalArray);
}

######################################################################################
# cmaCalcFitness
#
#' Calculate the fitness of a population.
#'
#' The population is usually obtained by \code{\link{cmaSamplePopulation}}.
#'
#' @param cma     CMA-ES Java object, already initialized with \code{\link{cmaInit}}
#' @param popR    a (\code{dimension x popSize}) matrix  from  \code{\link{cmaSamplePopulation}}
#' @param fitFunc a function to be minimized. Signature: accepts a vector \code{x}, returns a \code{double}.
#'
#' @return \code{fitness}, a vector of length \code{\link{cmaGetPopulationSize}(cma)} with the fitness of each individuum
#'
#' @examples
#'    cma <- cmaNew();
#'    cmaInit(cma,dimension=2,initialX=1.5);
#'    popR <- cmaSamplePopulation(cma);
#'    fitFunc <- function(x) {sum(x*x)};
#'    fitness <- cmaCalcFitness(cma,popR,fitFunc);
#'    cmaUpdateDistribution(cma,fitness);
#'
#' @seealso   \code{\link{cmaSamplePopulation}}, \code{\link{cmaUpdateDistribution}}, \code{\link{cmaNew}}
#' @author Wolfgang Konen, FHK, 2013
#' @export
######################################################################################
cmaCalcFitness <- function(cma,popR,fitFunc) {
      fitness = sapply(1:ncol(popR), function(i) { x=popR[,i]; fitFunc(x); });
}

######################################################################################
# cmaUpdateDistribution
#
#' Update CMA-ES distribution with the fitness vector of the last population.
#'
#'
#' @param cma CMA-ES Java object, already initialized with \code{\link{cmaInit}}
#' @param fitness vector of length \code{\link{cmaGetPopulationSize}(cma)} with the fitness of each individuum
#'
#' @note
#'  As a side effect, the CMA-ES Java object cma of class \code{CMAEvolutionStrategy}
#'  is augmented.
#'
#' @seealso   \code{\link{cmaSamplePopulation}}, \code{\link{cmaNew}}, \code{\link{cmaOptimDP}}
#' @author Wolfgang Konen, FHK, 2013
#' @export
######################################################################################
cmaUpdateDistribution <- function(cma,fitness) {
  rJava::.jcall(cma,,"updateDistribution",fitness);
}


