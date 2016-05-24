# Changes:
# 2013-02-17: Added the seed parts into the subfunction pp

#-----------------------------------------------------------------------------------------------
# Begin Javastuff preparations

.onLoad <- 
  function(libname, pkgname) {
    .jpackage(pkgname, lib.loc = libname)
  }

eppSession<-function(){
	eppSession <- .jnew("epp/EPPLab")
	return(eppSession)
}
# End of Javastuff preparations
#-----------------------------------------------------------------------------------------------


# Internal function that calls the Java Functions

 pp <- function(data,index,alg,nSimulations,iterations=NULL,individuals=NULL,particles=NULL, step_iter, eps){
 
  seed <- runif(1)*10^7

  indexInt <- 0

  rows <- nrow(data)
  cols <- ncol(data)

  index<-match.arg(index,c("KurtosisMax","Friedman","Discriminant","FriedmanTukey","KurtosisMin","Indice4","StahelDonoho"))
  alg<-match.arg(alg,c("GA","PSO","Tribe"))
   
  # Possible options for alg are "GA","Tribes" and "PSO" 
  if(alg=="GA")
  {
    ifelse(is.null(iterations),iterations <- 50, iterations <- iterations)
    ifelse(is.null(individuals),tPara <- 20, tPara <- individuals)
  } else if(alg=="PSO"){
    ifelse(is.null(iterations),iterations <- 200, iterations <- iterations)
    ifelse(is.null(particles),tPara <- 50, tPara <- particles)
  } else if(alg=="Tribe"){
    ifelse(is.null(iterations),iterations <- 200, iterations <- iterations)
    ifelse(is.null(particles),tPara <- 50, tPara <- particles)
  } else {
    stop("Unknown algorithm!\n")
  }

  # Possible Indices as named in the Java code. For easing reasons we switch them to numerical values.
   if(index=="KurtosisMax"){
    indexInt <- 1
  } else if(index=="Friedman"){
    indexInt <- 2
  } else if(index=="Discriminant"){
    indexInt <- 3
  } else if(index=="FriedmanTukey"){
    indexInt <- 4
  } else if(index=="KurtosisMin"){
    indexInt <- 5
 # } else if(index=="Indice4"){
 #   indexInt <- 6
 # } else if(index=="StahelDonoho"){
 #   indexInt <- 7
  } else {
    stop("Unknown (or currently disabled) index!\n")
  }

 eppLabString <- .jcall("epp/EPPLab",returnSig="S",method="eppLabRInterface",indexInt,nSimulations,alg,iterations,tPara,as.double(rows),as.double(cols),as.double(as.vector(t(data))),as.double(seed), as.double(step_iter), as.double(eps))
 eppLabString

}


EpplabOutputConv <- function(x, maxiter)
    {
    xList <- strsplit(x,"\n")
    xList <- xList[[1]]
    indI <-grep("<I>",xList)
    indA <-grep("<A>",xList)
    indNumberIter <- grep("<Number_of_iterations>",xList)
    indHasConverged <- grep("<Has_converged>",xList)    

    Apart <- xList[indA]
    Apart <- gsub("<A> ","",Apart)
    Apart <- gsub(",",".",Apart)
    Apart <- strsplit(Apart," ")
    Apart <- do.call(cbind,Apart)
    Apart <- apply(Apart,2,as.numeric)
    
    Ipart <- xList[indI]
    Ipart <- gsub("<I> ","",Ipart)
    Ipart <- gsub(",",".",Ipart)
    Ipart <- as.numeric(Ipart)

    NumberIterPart <- xList[indNumberIter]
    NumberIterPart <- gsub("<Number_of_iterations> ","",NumberIterPart)
    NumberIterPart <- gsub(" </Number_of_iterations>","",NumberIterPart)
    NumberIterPart <- as.numeric(unlist(strsplit(NumberIterPart," ")))

    HasConvergedPart <- xList[indHasConverged]
    HasConvergedPart <- gsub("<Has_converged> ","",HasConvergedPart)
    HasConvergedPart <- gsub(" </Has_converged>","",HasConvergedPart)
    HasConvergedPart <- unlist(strsplit(HasConvergedPart," "))

    #converged <- NumberIterPart < maxiter
    HasConvergedPart <- ifelse(HasConvergedPart=="true",TRUE,FALSE)
    converged <- HasConvergedPart
    Ipart.c <- Ipart[converged]
    Ipart.nc <- Ipart[!converged] 
    Apart.c <- as.matrix(Apart[,converged])
    Apart.nc <- as.matrix(Apart[,!converged])
    NumberIterPart.c <- NumberIterPart[converged]
    NumberIterPart.nc <- NumberIterPart[!converged] 
    HasConvergedPart.c <- HasConvergedPart[converged]
    HasConvergedPart.nc <- HasConvergedPart[!converged]

    orderI.c <- order(Ipart.c, decreasing=TRUE)
    Ipart.c <- Ipart.c[orderI.c]
    Apart.c <- as.matrix(Apart.c[,orderI.c])
    NumberIterPart.c <- NumberIterPart.c[orderI.c]
    HasConvergedPart.c <- HasConvergedPart.c[orderI.c]

    orderI.nc <- order(Ipart.nc, decreasing=TRUE)
    Ipart.nc <- Ipart.nc[orderI.nc]
    Apart.nc <- as.matrix(Apart.nc[,orderI.nc])
    NumberIterPart.nc <- NumberIterPart.nc[orderI.nc]
    HasConvergedPart.nc <- HasConvergedPart.nc[orderI.nc]

    Ipart <- c(Ipart.c,Ipart.nc)
    Apart <- cbind(Apart.c,Apart.nc)
    NumberIterPart <- c(NumberIterPart.c,NumberIterPart.nc)
    HasConvergedPart <- c(HasConvergedPart.c, HasConvergedPart.nc)

# Uncomment this if the order should be given with respect to value <I> (maybe this could be also an option in the print function?!)
#     orderI <- order(Ipart, decreasing=TRUE)
#     Ipart <- Ipart[orderI]
#     Apart <- Apart[,orderI]
#     NumberIterPart <- NumberIterPart[orderI]

# Order according to convergence
#     orderNumbIter <- order(NumberIterPart, decreasing=FALSE)
#     Ipart <- Ipart[orderNumbIter]
#     Apart <- Apart[,orderNumbIter]
#     NumberIterPart <- NumberIterPart[orderNumbIter]

    if(sum(!converged)>0) {
      if(sum(!converged)==1){ 
        warning("There was ",sum(!converged)," non-converged simulation run!")
      } else {
        warning("There were ",sum(!converged)," non-converged simulation runs!")
      }
    }

    list(PPindex = Ipart, PPdir = Apart, PPiter = NumberIterPart, PPconv= HasConvergedPart)
    }




#' Function for Exploratory Projection Pursuit.
#' 
#' REPPlab optimizes a projection pursuit (PP) index using a Genetic Algorithm
#' (GA) or one of two Particle Swarm Optimisation (PSO) algorithms over several
#' runs, implemented in the Java program EPP-lab.  One of the PSO algorithms is
#' a classic one while the other one is a parameter-free extension called
#' Tribes. The parameters of the algorithms (maxiter and individuals for GA and
#' maxiter and particles for PSO) can be modified by the user.  The PP indices
#' are the well-known Friedman and Friedman-Tukey indices together with the
#' kurtosis and a so-called discriminant index that is devoted to the detection
#' of groups.  At each run, the function finds a local optimum of the PP index
#' and gives the associated projection direction and criterion value.
#' 
#' The function always centers the data using \code{\link{colMeans}} and
#' divides by the standard deviation. Sphering the data is optional. If
#' sphering is requested the function \code{\link{WhitenSVD}} is used, which
#' automatically tries to determine the rank of the data.
#' 
#' Currently the function provides the following projection pursuit indices:
#' \code{KurtosisMax}, \code{Discriminant}, \code{Friedman},
#' \code{FriedmanTukey}, \code{KurtosisMin}.
#' 
#' Three algorithms can be used to find the projection directions. These are a
#' Genetic Algorithm \code{GA} and two Particle Swarm Optimisation algorithms
#' \code{PSO} and \code{Tribe}.
#' 
#' Since the algorithms might find local optima they are run several times. The
#' function sorts the found directions according to the optimization criterion.
#' 
#' The different algorithms have different default settings. It is for GA:
#' \code{maxiter=50} and \code{individuals=20}. For PSO: \code{maxiter=20} and
#' \code{particles=50}. For Tribe: \code{maxiter=20}.
#' 
#' For GA, the size of the generated population is fixed by the user
#' (individuals). The algorithm is based on a tournament section of three
#' participants.  It uses a 2-point crossover with a probability of 0.65 and
#' the mutation operator is applied to all the individuals with a probability
#' of 0.05. The termination criterion corresponds to the number of generations
#' and is also fixed by the user (maxiter).
#' 
#' For PSO, the user can give the number of initial generated particles and
#' also the maximum number of iterations. The other parameters are fixed
#' following Clerc (2006) and using a "cosine" neighborhood adapted to PP for
#' the PSO algorithm. For Tribes, only the maximum number of iterations needs
#' to be fixed. The algorithm proposed by Cooren and Clerc (2009) and adapted
#' to PP using a "cosine neighborhood" is used.
#' 
#' The algorithms stop as soon as one of the two following conditions holds:
#' the maximum number of iterations is reached or the relative difference
#' between the index value of the present iteration i and the value of
#' iteration i-\code{step_iter} is less than \code{eps}. In the last situation,
#' the algorithm is said to converge and \code{EPPlab} will return the number
#' of iterations needed to attain convergence.  If the convergence is not
#' reached but the maximum number of iterations is attained, the function will
#' return some warnings. The default values are 10 for \code{step_iter} and
#' \code{1E-06} for \code{eps}. Note that if few runs have not converged this
#' might not be problem and even non-converged projections might reveal some
#' structure.
#' 
#' @param x Matrix where each row is an observation and each column a
#' dimension.
#' @param PPindex The used index, see details.
#' @param PPalg The used algorithm, see details.
#' @param n.simu Number of simulation runs.
#' @param sphere Logical, sphere the data. Default is \code{FALSE}, in which
#' case the data is only standardized.
#' @param maxiter Maximum number of iterations.
#' @param individuals Size of the generated population in GA.
#' @param particles Number of generated particles in the standard PSO
#' algorithm.
#' @param step_iter Convergence criterium parameter, see details. (Default: 10)
#' @param eps Convergence criterium parameter, see details. (Default: 10^(-6))
#' @return A list with class 'epplab' containing the following components:
#' \item{PPdir}{Matrix containing the PP directions as columns, see details.}
#' \item{PPindexVal}{Vector containing the objective criterion value of each
#' run.} \item{PPindex}{Name of the used projection index.}
#' \item{PPiter}{Vector containing the number of iterations of each run.}
#' \item{PPconv}{Boolean vector. Is TRUE if the run converged and FALSE else.}
#' \item{PPalg}{Name of the used algorithm.} \item{maxiter}{Maximum number of
#' iterations, as given in function call.} \item{x}{Matrix containing the data
#' (centered!).} \item{sphere}{Logical} \item{transform}{The transformation
#' matrix from the whitening or standardization step.} \item{backtransform}{The
#' back-transformation matrix from the whitening or standardization step.}
#' \item{center}{The mean vector of the data}
#' @author Daniel Fischer, Klaus Nordhausen
#' @references \cite{Larabi Marie-Sainte, S., (2011), Biologically inspired
#' algorithms for exploratory projection pursuit, PhD thesis, University of
#' Toulouse.}
#' 
#' \cite{Ruiz-Gazen, A., Larabi Marie-Sainte, S. and Berro, A. (2010),
#' Detecting multivariate outliers using projection pursuit with particle swarm
#' optimization, \emph{COMPSTAT2010}, pp. 89-98.}
#' 
#' \cite{Berro, A., Larabi Marie-Sainte, S. and Ruiz-Gazen, A. (2010). Genetic
#' algorithms and particle swarm optimization for exploratory projection
#' pursuit. Annals of Mathematics and Artifcial Intelligence, 60, 153-178.}
#' 
#' \cite{Larabi Marie-Sainte, S., Berro, A. and Ruiz-Gazen, A. (2010). An
#' effcient optimization method for revealing local optima of projection
#' pursuit indices. \emph{Swarm Intelligence}, pp. 60-71.}
#' 
#' \cite{Clerc, M. (2006). Particle Swarm Optimization. ISTE, Wiley.}
#' 
#' \cite{Cooren, Y., Clerc, M. and Siarry, P. (2009). Performance evaluation of
#' TRIBES, an adaptive particle swarm optimization algorithm. Swarm
#' Intelligence, 3(2), 149-178.}
#'
#' @examples
#' 
#'   library(tourr)
#'   data(olive)
#'   olivePP <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMax",n.simu=5, maxiter=20)
#'   summary(olivePP)
#' 
#'   library(amap)
#'   data(lubisch)
#'   X <- lubisch[1:70,2:7]
#'   rownames(X) <- lubisch[1:70,1]
#'   res <- EPPlab(X,PPalg="PSO",PPindex="FriedmanTukey",n.simu=15, maxiter=20,sphere=TRUE)
#'   print(res)
#'   summary(res)
#'   fitted(res)
#'   plot(res)
#'   pairs(res)
#'   predict(res,data=lubisch[71:74,2:7])
#' 
#' @export EPPlab
EPPlab <- function(x, PPindex="KurtosisMax", PPalg="GA", n.simu=20, sphere=FALSE, maxiter=NULL, individuals=NULL, particles=NULL, step_iter=10, eps=10^(-6))
        {
	# Check if row names are given
        if(is.null(rownames(x))) rownames(x) <- as.character(1:nrow(x))

	# center first the data:
	MEANS <- colMeans(x)
        x.c <- as.matrix(sweep(x,2,MEANS,"-"))
	
	# Whitening or centring the data
	if(sphere==TRUE)
	{		
	  x <- WhitenSVD(x)
	} else {
	  x <- x.c
          SDS <- apply(x,2,sd)
          x <- sweep(x,2,SDS, "/")
	  attr(x,"center") <- MEANS
	  attr(x,"transform") <- diag(1/SDS)
	  attr(x,"backtransform") <- diag(SDS)
	}

        jepplab <- pp(data=x, index=PPindex , alg=PPalg , nSimulations=n.simu , iterations=maxiter,individuals=individuals,particles=particles, step_iter=step_iter, eps=eps)
        AI <- EpplabOutputConv(jepplab, maxiter)
        xNames <- colnames(x)
        
        PPindexVal <- AI$PPindex
	# Give the whitened directions as PPDir
	    PPdir <- attributes(x)$transform %*% AI$PPdir
        rownames(PPdir) <- xNames
        colnames(PPdir) <- paste("Run",1:n.simu,sep="")
        
        RES <- list(PPdir=PPdir, PPindexVal=PPindexVal, PPindex=PPindex, PPiter= AI$PPiter,PPconv= AI$PPconv, PPalg=PPalg, x=as.matrix(x.c), sphered=sphere, whiteMat=attributes(x)$transform, backMat=attributes(x)$backtransform, center=attributes(x)$center, maxiter=maxiter)
        class(RES) <- "epplab"
	RES
     }

