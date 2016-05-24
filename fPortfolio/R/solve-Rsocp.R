
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  solveRsocp                Portfolio interface to solver Rsocp
#  .rsocpArguments           Returns arguments for solver
#  .rsocp                    Wrapper to solver function
#  .rsocpControl             Returns default controls for solver
################################################################################


solveRsocp <-
  function(data, spec, constraints)
  {
    # Description:
    #   Portfolio interface to solver Rsocp
    
    # Example:
    #   ans = solveRquadprog(.lppData, .mvSpec, "LongOnly")[-3]
    #   .mvSpec2 = .mvSpec; setTargetRisk(.mvSpec2) = ans$targetRisk
    #   solveRsocp(.lppData, .mvSpec2, "LongOnly")[-3]; ans
    
    #   efficientortfolio
    
    # FUNCTION:
    
    # Transform Data and Constraints:
    Data = portfolioData(data, spec)
    
    # Trace:
    trace <- getTrace(spec)
    if(trace) cat("\nPortfolio Optimiziation:\n Using Rsocp ...\n\n")
    
    # Get Specifications:
    nAssets = getNAssets(Data)
    
    # Create '.rsocp' conform arguments:
    args <- .rsocpArguments(data, spec, constraints)
    
    # Optimize:
    ans <- .rsocp(
      f = args$f, 
      A = args$A, 
      b = args$b, 
      C = args$C, 
      d = args$d, 
      N = args$N, 
      targetRisk = args$targetRisk, 
      mu = args$mu,
      Scale = args$Scale)
    
    # Return Value:
    class(ans) = c("solver", "list")
    ans
  }


################################################################################


.rsocpArguments <-
  function(data, spec, constraints)
  {
    # Description:
    #   Returns socp conform arguments for the solver
    
    # Example:
    #   .rsocpArguments(data, spec, constraints)
    
    # FUNCTION:
    
    # Settings:
    Data = portfolioData(data, spec)
    nAssets = getNAssets(Data)
    Scale = 1.0e6 * sd(as.vector(data))
    mu = getMu(Data) / Scale
    Sigma = getSigma(Data) / Scale^2
    targetRisk = getTargetRisk(spec) / Scale
    
    # Objective Function:
    f <- -mu
    
    # Constraints:
    eqsumW = eqsumWConstraints(data, spec, constraints)
    
    # C - Cone Constraints:
    C1 <- rep(0, nAssets)                                     # xCx
    C2 <- eqsumW[2, -1]                                       # sum(x)
    C3 <- rbind(diag(nAssets), -diag(nAssets) )               # x[i]>0
    
    # d - Cone Constraints:
    d1 <- targetRisk                                          # xCx = risk
    d2 <- eqsumW[2, 1]                                        # sum(x) <= 1
    d3 <- c(rep(0, nAssets), rep(-1, nAssets))                # x[i] > 0
    
    # A - Cone Constraints:
    A1 <- Rsocp::.SqrtMatrix(Sigma)
    A2 <- matrix(0, ncol = nAssets)
    A3 <- matrix(0, nrow = nrow(C3), ncol = nAssets)
    
    # b - Cone Constraints:
    b1 <- rep(0, nAssets)                                     # xCx
    b2 <- 0                                                   # sum(x)
    b3 <- rep(0, nrow(C3))                                    # x[i]>0
    
    # N - Cone Constraints:
    N1 <- nAssets                                             # dim(C)
    N2 <- 1                                                   # Full Investment
    N3 <- rep(1, nrow(C3))                                    # Long
    
    # Combine Constraints for SOCP:
    A <- rbind(A1, A2, A3)
    b <- c(b1, b2, b3)
    C <- rbind(C1, C2, C3)  
    d <- c(d1, -d2, -d3)  
    N <- c(N1, N2, N3)      
    
    # Return Value:
    list(f = f, A = A, b = b, C = C, d = d, N = N, 
         targetRisk = targetRisk * Scale, mu = mu * Scale, Scale = Scale)
  }


################################################################################


.rsocp  <-
  function(f, A, b, C, d, N, x = NULL, z = NULL, w = NULL, 
           targetRisk, mu = mu, Scale = Scale, control = .rsocpControl())
  {
    # Description:
    #   SOCP solver function for portfolios
    
    # Details:
    #   Package: fPortfolio
    #   Title: An R extenstion library to use SOCP from R.
    #   Version: 0.1 as of 2008-31-01
    #   Author: Yohan Chalabi and Diethelm Wuertz
    #   Description: Second-order cone programming solver 
    #       written by M. Lobo, L. Vandenberghe, and S. Boyd. 
    #       R.socp is a wrapper library to use it from R.
    #   License: R.socp - GPL
    
    # FUNCTION
    
    # Solve Portfolio:
    optim <- Rsocp::socp(f, A, b, C, d, N, x, z, w, control)
    
    # Extract Weights:
    weights = .checkWeights(optim$x)
    attr(weights, "invest") = sum(weights)
    
    # Prepare Output List: 
    ans <- list(
      type = "MV",
      solver = "solveRsocp",
      optim = optim,
      weights = weights,
      targetReturn = (weights %*% mu)[[1]],
      targetRisk = targetRisk,
      objective = (weights %*% mu)[[1]],
      status = as.integer(!optim$convergence),
      message = optim$message)
    
    # Return Value: 
    ans
  }


################################################################################


.rsocpControl <- 
  function(abs.tol = 1.0e-18, rel.tol = 1.0e-16, target = 0, 
           max.iter = 1000, Nu = 10, out.mode = 0, BigM.K = 2, BigM.iter = 5)
  { 
    # Description:
    #   Control list for portfolio SOCP optimization
    
    # FUNCTION:
    
    # Return Value:
    list(
      abs.tol = abs.tol, 
      rel.tol = rel.tol,
      target = target,
      max.iter = max.iter,
      Nu = Nu,
      out.mode = out.mode,
      BigM.K = BigM.K,
      BigM.iter = BigM.iter)
  }


################################################################################

