#active set methods for different solver


activeSet <- function(isomat, mySolver = "LS", x0 = NULL, ups = 1e-12, check = TRUE, maxiter = 100, ...) 
{
  a <- isomat                      #matrix with order restrictions
  if (ncol(isomat) != 2) stop("isomat must have 2 columns!")
  
  extra <- list(...)
  if (length(x0) == 0) x0 <- rep(0,max(isomat))
   
  n <- length(x0)
  xold <- x0                       #starting values
  a <- isomat[,c(2,1)]
  ax <- aTx(a, xold)               #difference between order restrictions

  if (any(ax < 0)) stop("Starting solution not feasible. Ax must be >= 0!")

  ia <- is.active(ax, ups = ups)   #which constraints are active 
  iter <- 0

  #--------------- solver specification --------------------
  mySolverDB <- list()
  mySolverDB$chebyshev <- mSolver
  mySolverDB$LS <- lsSolver
  mySolverDB$L1 <- dSolver
  mySolverDB$quantile <- pSolver
  mySolverDB$GLS <- lfSolver
  mySolverDB$poisson <- sSolver
  mySolverDB$Lp <- oSolver
  mySolverDB$asyLS <- aSolver
  mySolverDB$L1eps <- eSolver
  mySolverDB$huber <- hSolver
  mySolverDB$SILF <- iSolver

  if(is.character(mySolver)) {
      pos <- pmatch(tolower(mySolver),
                    tolower(names(mySolverDB)))
      if(is.na(pos))
          stop(gettextf("Invalid skmeans method '%s'.", mySolver))
      mySolver <- mySolverDB[[pos]]
   }
  #------------ end solver specification -----------------


  #-------------- start active set iterations ------------------------
  repeat {
    iter <- iter + 1
    if (length(ia)==0) {           #no set active (typically 1st iteration)
      aia <- NULL                  
    } else {                       #active set
      aia <- a[ia,]                #active constraint
    }
    yl <- mySolver(xold, aia, extra)  #call solver
    
    y <- yl$x                      #fitted values
    lbd <- yl$lbd                  #Lagrange multiplier (KKT vector lambda)
    fy <- yl$f                     #value target function
    gy <- yl$gx                    #gradient
    ay <- aTx(a,y)                 #compute Ax
    iy <- which.min(ay)            #restriction with the largest violation of Ax >= 0
    my <- ay[iy]                   #value of this restriction (worst one)
    
    if (length(lbd)==0) {          #no lambda
      ml <- Inf
    } else {
      il <- which.min(lbd)         #index minimum Lagrange 
      ml <- lbd[il]                #value minimum Lagrange (worst one)
    }
    
    if (is.pos(my, ups)) {         #no violation of Ax >= 0 (feasible)
      if (is.pos(ml, ups)) break() #convergence reached, all lambda >= 0
      xnew <- y                    #Ax >= 0 ok; lamba >= 0 not
      ax <- ay
      ia <- ia[-il]                #add constraint to active set with worst lambda
    } else {                       #still constraint violations (infeasible)
      k <- which((ax>-ups) & (ay<ups))  #index where we have violations
      #k <- which((ax>0) & (ay<0))  #index where we have violations
      rat <- -ay[k]/(ax[k]-ay[k])  #line search starts
      ir <- which.max(rat)
      alw <- rat[ir]               #alpha
      xnew <- y+alw*(xold-y)       #update function values
      ax <- aTx(a,xnew)            #update constraints Ax
      ia <- sort(c(ia,k[ir]))      #collect active sets
    }
    xold <- xnew                   #end iteration, start new one
    if (iter == maxiter) {
      warning("Maximum number of iterations reached!")
      break()
    }
    
  }
  #---------------------- end active set iterations --------------------
  options(warn = 0)

  lup <- rep(0, length(ay)) 
  lup[ia] <- lbd                   #final vector of lambdas (0 where there was no active set)
  hl <- taTx(a, lup, n)            #A'lambda (should be equal to gradient)
  if (check) {                     #check KKT
    ck <- checkSol(y, gy, a, ay, hl, lup, ups)   #checks feasibility of 4 KKT conditions
    ck <- as.list(ck)
    names(ck) <- c("stationarity","primal.feasibility","dual.feasibility","complementary.slackness")
  } else {
    ck <- NULL
  }

  result <- list(x = y, y = extra$y, lambda = lup, fval = fy, constr.val = ay, Alambda = hl, 
  gradient = gy, isocheck = ck, niter = iter, call = match.call())
  class(result) <- "activeset"
  result
}
