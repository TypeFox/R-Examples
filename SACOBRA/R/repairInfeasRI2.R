######################################################################################
# repairInfeasRI2
#
#'  Repair an infeasible solution with the method RI2
#' 
#'  If the solution x is infeasible, i.e. if \code{any(fReal>0)==TRUE}:
#'  \enumerate{
#'    \item Estimate the gradient of the constraint surrogate function(s) (go a tiny step in each dimension 
#'      in the direction of constraint increase).
#'    \item 
#     OLD=TRUE (not recommended): Move a step into the opposite direction for every violated
#       constraint and sum these steps, such that the new x has all true constraint   
#       responses approximately 0 (or slightly negative, i.e. definitely feasible, if 
#       \code{ri$kappa>1.0}).
#     OLD=FALSE (recommended): 
#'      Take \code{cobra$ri$mmax} random realizations in the 'feasible parallelepiped'
#'      and select among them the best feasible solution, based on the surrogates,
#'    \item Check whether the new solution is for every dimension in the bounds 
#'      \code{[cobra$lower,cobra$upper]} of the search region. 
#'      If not, set the gradient to 0 in these dimensions and re-iterate from step 2. 
#'  }
#'  There is no guarantee but a good chance, that the returned solution \code{z} will be feasible.
#'  
#'  For further details see  [Koch15a] Koch, P.; Bagheri, S.; Konen, W. et al. "A New Repair Method For Constrained 
#'  Optimization". Proc. 17th Genetic and Evolutionary Computation Conference (GECCO), 2015.
#'      
#'  @param x            an infeasible solution (vector of length \code{dimension})
#'  @param fReal        a vector of length nconstraint holding the real constraint values at \code{x}
#'  @param rbf.model    the constraint surrogate models 
#'  @param cobra        parameter list, we need here 
#'     \describe{
#'       \item{\code{lower}}{    lower bounds of search region}
#'       \item{\code{upper}}{    upper bounds of search region}
#'       \item{\code{ri}}{       a list with all parameters for \code{repairInfeasRI2}  }
#'       \item{\code{trueFuncForSurrogate}}{ if TRUE (only for diagnostics), use the true constraint
#'                      functions \code{conFunc} instead of the constraint surrogate models 
#'                      \code{rbf.model} }
#'     }
#'  @param checkIt      [FALSE] if TRUE, perform a check whether the returned solution is really 
#'                      feasible. Needs access to the true constraint function \code{conFunc}
#'  @param conFunc      [NULL] function returning real constraint vector (needed only
#'                      if \code{checkIt==T} or if \code{cobra$trueFuncForSurrogate==T})
#'  @return \code{z},  a vector with a repaired (hopefully feasible) solution
#'  @seealso   \code{\link{repairChootinan}}, \code{\link{cobraPhaseII}}
#'
#'  @author Wolfgang Konen, Cologne Univeristy of Applied Sciences
#'  @export
######################################################################################
repairInfeasRI2 <- function(x,fReal,rbf.model,cobra,checkIt=FALSE,conFunc=NULL) 
{  
  testit::assert("cobra$ri is NULL (not defined)", !is.null(cobra$ri) );
  testit::assert("cobra$ri$q is NULL (not defined)", !is.null(cobra$ri$q) );
  
  ri = cobra$ri
  if (ri$OLD) {
    ri$eps1=0.0
    ri$eps2=0.0
    if (is.null(ri$kappa)) ri$kappa=1.2
    cat("ri$OLD=TRUE: Setting defaults for ri: eps1=eps2=0\n")
  }
  trueFunc = cobra$trueFuncForSurrogate
  lowerP = cobra$lower
  upperP = cobra$upper
  gradEps = cobra$ri$gradEps
  if (is.null(gradEps)) gradEps = min(upperP-lowerP)/1000; 
  # gradEps = stepsize for numerical gradient calculation
  # e.g. 0.001 if the smallest length of search cube is 1.0
  
  #########################################################################
  # helper functions
  #
  
  # isEpsilonFeasible: 
  # return TRUE, if all constraint surrogates have a feasibility 
  # 'better' than eps2 (i.e. they are -eps2 or lower)
  isEpsilonFeasible <- function(y,eps2,rbf.model) {
    if (!trueFunc) {
      fRbf <- interpRBF(y,rbf.model)
    } else {
      fRbf <- conFunc(y)      
    }
    if (any(fRbf+eps2>0)) return(FALSE);
    return(TRUE);
  }
 
  # findBestInfeasible: 
  # Given a set of k points x + deltaMat[k,], k=1,...,nrow(deltaMat). 
  # Find this point x + deltaMat[kBest,] which has
  #   a) the lowest number of eps2-infeasibilities
  #   b) if there is more than one such point, the one with the lowest maximum violation
  #      (if there are more than one point with the same lowest max violation, select the first)
  # Return delta[kBest,]
  findBestInfeasible <- function(deltaMat,x,eps2,rbf.model) {
    numMaxViol <- function(k) {
      if (!trueFunc) {
        fRbf <- interpRBF(x+deltaMat[k,],rbf.model)
      } else {
        fRbf <- conFunc(x+deltaMat[k,])     
      }
      ind <- which(fRbf+eps2>0)
      return(list(numViol=length(ind)
                 ,maxViol=max(fRbf)))            
    }
    res = sapply(1:nrow(deltaMat),numMaxViol)
    numViol = unlist(res[1,])
    maxViol = unlist(res[2,])
    maxViol[numViol>min(numViol)] <- Inf  # invalidate all entries with too high numViol
    kBest=which(maxViol==min(maxViol))[1]
    
    return(deltaMat[kBest,])
  }
  
  # selectBest:
  # among the eps2-feasible points (rows of S) select the one with minimal length  
  selectBest <- function(S,x,ri) {
    L2 = rowSums(S*S)
    ind = which(L2==min(L2))[1]
    return(S[ind,])
  }
  
  # checkSingleConstraints (for debug only):
  # Does the single correction Del[k,] calculated for kth violated constraint bring
  # the solution close to the kth constraint boundary, as it should?
  # If so, abs(fSingle[k]) should be much smaller than fBefore[k]
  checkSingleConstraints <- function(Del,Grad,Viol) {
    cat("\n")
    cat("   Length of Del vectors:\n")
    print(sqrt(rowSums(Del*Del)));
    cat("   Length of Grad vectors:\n")
    print(sqrt(rowSums(Grad*Grad)));
    fBefore = NULL
    fSingle = NULL
    for (k in 1:length(Viol)) {
      v = Viol[k]
      fBefore = c(fBefore,interpRBF(x,rbf.model)[v])
      fSingle = c(fSingle,interpRBF(x+Del[k,],rbf.model)[v])
    }
    cat("   single violations before/after single correction: \n"); 
    print(rbind(fBefore,fSingle));
    cat("   quotient q=fBefore/fSingle (|q| should be larger than approx. 8.0): \n")
    print(fBefore/fSingle);     
  }
  
  checkBestInfeasible <- function(Del) {
    cat("No ri$eps2-feasible point - select the best infeasible point\n")
    # print some debug info
    x_old = x+ri$kappa*colSums(Del);
    fx=interpRBF(x,rbf.model); indx=which(fx+ri$eps2>0);        # infeasible point
    fy=interpRBF(x+Delta,rbf.model); indy=which(fy+ri$eps2>0);  # best infeasible repair 
    fw=interpRBF(x_old,rbf.model); indw=which(fw+ri$eps2>0);    # old repair mechanism
    print(indx); print(indy); print(indw);
  }
  
  #########################################################################
  # start of repairInfeasRI2
  #
  verbosecat(cobra$verbose,important=FALSE,"RI2: repairing the infeasible result ...\n")
  dimension <- length(x)
  nconstraint <- length(fReal) # ncol(rbf.model$coef)
  if (!is.null(rbf.model)) testit::assert("nconstraint", nconstraint==ncol(rbf.model$coef))
  rownames(fReal) <- NULL
  
  nd <- data.frame(outer(rep(1,2*dimension+1),x))
  for (i in 1:dimension) {
    nd[2*i  ,i] <- nd[2*i  ,i] - gradEps
    nd[2*i+1,i] <- nd[2*i+1,i] + gradEps
  }
  
  
  # f is a (2*dimension+1 x nconstraint) matrix containing constraint surrogate 
  # responses at x (row 1) and small '-' and '+' deviations from x for any dimension  
  # in rows 2,...,2*dimension+1
  if (!trueFunc) {
    f <- as.matrix(predict(rbf.model,newdata=nd))
  } else {
    f <- t(apply(nd,1,conFunc))     
  }
  
  testit::assert("Wrong types for f or fReal", is.matrix(f), is.vector(fReal) );
  testit::assert("Columns do not match in f and fReal", ncol(f) == length(fReal))
  
  ix = rep(FALSE,dimension)
  # ix: a boolean vector of length dimension. It indicates which dimensions of the gradient
  # should be set to zero (because a repair in this dimension would cause the repaired solution 
  # to leave the search region [lowerP,upperP]). Initially, all elements of ix are false, i.e. 
  # every dimension is taken.

  while(1) {
    Del <- NULL
    Grad <- NULL
    Viol <- NULL
    for (k in 1:ncol(f)) {
      if (fReal[k]+ri$eps1>0) {   # if the kth constraint is not ri$eps1-feasible
        gradf <- rep(0,dimension)
        for (i in 1:dimension) {
          if (f[2*i,k]>f[2*i+1,k])    # if the penalty increase is larger in direction '-gradEps':
            gradf[i] <- (f[2*i  ,k]-f[1,k])/(-gradEps)
          else                        # else, i.e. '+gradEps':
            gradf[i] <- (f[2*i+1,k]-f[1,k])/gradEps
        }
        gradf[ix] <- 0            # zero all dimensions which led to out-of-search-region
                                  # repairs in previous iterations
        g2 <- sum(gradf*gradf)  
        if (g2==0) {
          warning("Cannot repair infeasible solution w/o moving out of search region")
          # Return the incoming (infeasible) solution instead:
          z = x
          return(z)      
        }
        Del_k =  -gradf * (fReal[k]+ri$eps1)/g2  
        Del = rbind(Del,Del_k)
        # Del is a matrix with as many rows as there are eps1-infeasible constraints.
        # The kth row of Del contains the step suggested for the kth constraint.
        Grad = rbind(Grad,gradf)
        # Similarly, the kth row of Grad contains the gradient for the kth constraint.
        
        Viol = c(Viol,k)  # Viol is the list of violated constraint numbers
      }
    } # for k
    #checkSingleConstraints(Del,Grad,Viol);
    #browser()
    #if (checkIt) print(sqrt(rowSums(Del*Del)));
    
    if (ri$OLD) 
    {
      # this was the old version (before 2014-09-29, as in repairInfeasible-OLD.r)
      Delta = colSums(Del)
      z = x + ri$kappa*Delta      
    } else {
      # this is the new repairInfeasible mechanism (after 2014-09-29, see
      # Notes.d/presentation/present-Wolfgang-2014-09-24-RepairInfeas2): 
      S = NULL
      deltaMat = NULL
      for (m in 1:ri$mmax) {
        alpha = runif(nrow(Del))*ri$q         # random coef. from distribution U[0,a]
        alphaMat = outer(alpha,rep(1,ncol(Del)))
        delta = colSums(alphaMat*Del)
        deltaMat = rbind(deltaMat,delta)
        if (isEpsilonFeasible(x+delta,ri$eps2,rbf.model))
          S = rbind(S,delta)
      }
      if (is.null(S)) {
        # No ri$eps2-feasible point - return the best infeasible solution instead:
        Delta = findBestInfeasible(deltaMat,x,ri$eps2,rbf.model)
        if (checkIt) checkBestInfeasible(Del);
      } else {
        Delta = selectBest(S,x,ri)
      } # else (is.null(S))
      z = x + Delta
    } # else (ri$OLD)
    
    # This should normally not happen:
    testit::assert("New solution z contains NA or NaN!", !any(is.na(z)))
     
    ix2 = (z > upperP) | (z < lowerP) 
    if (! any(ix2)) {
      # we are done: z is inside search region in every dimension
      break # out of while           
    }
    ix = ix2 | ix   # don't forget the dimensions which were 'outside' in previous iterations
  } # while
  
  if (checkIt) {
    fRbf = interpRBF(z,rbf.model)       # fRbf:  constraint surrogate values after repair
    fTrue = conFunc(z)                  # fTrue: true constraint values after repair
    #print(fReal); print(fRbf); print(fTrue)
    violatedConstraints = which(fTrue>0)
    cfcReal <- maxReal <- 0;
    if (any(fTrue>0)) {
      cfcReal = sum(fTrue[violatedConstraints]);
      maxReal = max(fTrue[violatedConstraints]);
    }
    #checkSingleConstraints(Del,Grad,Viol);
    cat("Repaired solution is feasible: ",cfcReal<=0,", cfcReal=",cfcReal,", maxViol=",maxReal, "\n")
    #if (cfcReal>0) {
      cat("   eps1-inf constraints before repair: ",paste(which(fReal+ri$eps1>0) ,collapse=" "),"\n")
      cat("   violated constraints before repair: ",paste(which(fReal>0) ,collapse=" "),"\n")
      cat("   violated constraints  after repair: ",paste(which(fTrue>0),collapse=" "),"\n")
      cat("   violated c-surrogats  after repair: ",paste(which(fRbf>0),collapse=" "),"\n")
    #}
    
  }
  
  #browser()
  
  return (z)
}

######################################################################################
# repairInfeasibleW
#
#'  Wrapper for \code{\link{repairInfeasRI2}} (needed by RBFsearch.R).
# 
#'  @param resNM        a list as returned from optimizer (i.e. nmkb) with an infeasible 
#'                      solution in \code{x = resNM$par}
#'  @param fReal        a (1 x nconstraint) matrix holding the real constraint values at \code{x}
#'  @param rbf.model    the constraint surrogate models 
#'  @param cobra        parameter list, we need here 
#'         lower        lower bounds of search region
#'         upper        upper bounds of search region
#'         ri           a list with all parameters for repairInfeasRI2
#'  @param checkIt      [FALSE] if TRUE, perform a check whether the returned solution is really feasible.
#'                      Needs access to the true constraint function \code{conFunc}
#'  @param conFunc      [NULL] function returning real constraint vector      (needed if checkIt==T)
#'   
#'   @return \code{resRI}  a list containing:
#'      \item{\code{par}}{ the repaired (feasible) solution }
#'      \item{\code{value}}{ copied from \code{resNM$value} }
#'      \item{\code{feval}}{ 0 (to indicate that this solution comes from \code{repairInfeasRI2}) }
#'      \item{\code{convergence}}{ 0 }
#' @seealso   \code{\link{repairInfeasRI2}}, \code{RBFsearch}
#' @keywords internal
######################################################################################
repairInfeasibleW <- function(resNM,fReal,rbf.model,cobra,checkIt=FALSE,conFunc=NULL) 
{  
  x <- resNM$par
  
  z <- repairInfeasRI2(x,fReal,rbf.model,cobra,checkIt,conFunc) 
  
  resRI <- resNM
  resRI$par <- z
  resRI$convergence <- 0
  resRI$feval <- 0
  return(resRI)  
}


