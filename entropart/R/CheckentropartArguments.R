CheckentropartArguments <-
function() {

  # Get the list of arguments of the parent function
  ParentFunction <- sys.call(-1)[[1]]
  # If apply() or similar was used, the function name is not in ParentFunction: sys.call(-1)[[1]] returns "FUN"
  if (ParentFunction == "FUN") {
    warning("Function arguments cannot be checked, probably because you used apply(). Add CheckArguments=FALSE to suppress this warning.")
    return (TRUE)
  }
  
  ErrorFunction <- paste("Error in ", ParentFunction, ":")
  Args <- formals(match.fun(ParentFunction))

  ErrorMessage <- function(Message, Argument) {
    cat(deparse(substitute(Argument)), "cannot be:\n")
    print(utils::head(Argument))
    cat(paste(ErrorFunction, Message))
    stop("Check the function arguments.", call. = FALSE)
  }

  # Correction 
  if (!is.na(names(Args["Correction"]))) {
    Correction <- eval(expression(Correction), parent.frame())
    if (!is.character(Correction))
      ErrorMessage("Correction must be a string.", Correction)
  }
  # RCorrection 
  if (!is.na(names(Args["RCorrection"]))) {
    RCorrection <- eval(expression(RCorrection), parent.frame())
    if (!is.character(RCorrection))
      ErrorMessage("RCorrection must be a string.", RCorrection)
  }
  
  # MC 
  if (!is.na(names(Args["MC"]))) {
    MC <- eval(expression(MC), parent.frame())
    if (!is.MetaCommunity(MC))
      ErrorMessage("MC must be a MetaCommunity object.", MC)
  }
  
  # MClist
  if (!is.na(names(Args["MClist"]))) {
    MClist <- eval(expression(MClist), parent.frame())
    if (!is.list(MClist))
      ErrorMessage("MClist must be a list.", MClist)
    if (any(!unlist(lapply(MClist, function(x) is.MetaCommunity(x)))))
      ErrorMessage("All elements of MClist must be of class MetaCommunity.", MClist)
  }

  # q 
  if (!is.na(names(Args["q"]))) {
    q <- eval(expression(q), parent.frame())
    if (!is.numeric(q) | length(q)!=1)
      ErrorMessage("q must be a number.", q)
  }
  # q.seq 
  if (!is.na(names(Args["q.seq"]))) {
    q.seq <- eval(expression(q.seq), parent.frame())
    if (!is.vector(q.seq))
      ErrorMessage("q.seq must be a numeric vector.", q.seq)
  }

  # alpha
  if (!is.na(names(Args["alpha"]))) {
    alpha <- eval(expression(alpha), parent.frame())
    if (!is.numeric(alpha) | length(alpha)!=1)
      ErrorMessage("alpha must be a number.", alpha)
    if (any(alpha < 0))
      ErrorMessage("alpha must be positive.", alpha)
  }
  
  # Alpha 
  if (!is.na(names(Args["Alpha"]))) {
    Alpha <- eval(expression(Alpha), parent.frame())
    if (!is.numeric(Alpha))
      ErrorMessage("Alpha must be a number.", Alpha)    
    if (Alpha <= 0 | Alpha >= 1)
      ErrorMessage("Alpha must be strictly between 0 and 1.", Alpha)    
  }

  # BootstrapMethod 
  if (!is.na(names(Args["BootstrapMethod"]))) {
    BootstrapMethod <- eval(expression(BootstrapMethod), parent.frame())
    if (!is.character(BootstrapMethod))
      ErrorMessage("BootstrapMethod must be a string.", BootstrapMethod)
  }

  # Estimator 
  if (!is.na(names(Args["Estimator"]))) {
    Estimator <- eval(expression(Estimator), parent.frame())
    if (!is.character(Estimator))
      ErrorMessage("Estimator must be a string.", Estimator)
  }
  # CEstimator 
  if (!is.na(names(Args["CEstimator"]))) {
    CEstimator <- eval(expression(CEstimator), parent.frame())
    if (!is.character(CEstimator))
      ErrorMessage("CEstimator must be a string.", CEstimator)
  }
  
  # k
  if (!is.na(names(Args["k"]))) {
    k <- eval(expression(k), parent.frame())
    if (!is.numeric(k) | length(k)!=1)
      ErrorMessage("k must be a number.", k)
    if (k < 1)
      ErrorMessage("k must be at least 1.", k)
    if (as.integer(k) != k)
      ErrorMessage("k must be an integer.", k)
  }

  # mean
  if (!is.na(names(Args["mean"]))) {
    mean <- eval(expression(mean), parent.frame())
    if (!is.numeric(mean) | length(mean)!=1)
      ErrorMessage("mean must be a number.", mean)
  }

  # n
  if (!is.na(names(Args["n"]))) {
    n <- eval(expression(n), parent.frame())
    if (!is.numeric(n) | length(n)!=1)
      ErrorMessage("n must be a number.", n)
    if (any(n < 1))
      ErrorMessage("n must be at least 1.", n)
    if (as.integer(n) != n)
      ErrorMessage("n must be an integer.", n)
  }
  
  # Normalize 
  if (!is.na(names(Args["Normalize"]))) {
    Normalize <- eval(expression(Normalize), parent.frame())
    if (!is.logical(Normalize))
      ErrorMessage("Normalize must be TRUE or FALSE.", Normalize)
  }

  # nPoints
  if (!is.na(names(Args["nPoints"]))) {
    nPoints <- eval(expression(nPoints), parent.frame())
    if (!is.numeric(nPoints) | length(nPoints)!=1)
      ErrorMessage("nPoints must be a number.", nPoints)
    if (any(nPoints < 1))
      ErrorMessage("nPoints must be at least 1.", nPoints)
    if (as.integer(nPoints) != nPoints)
      ErrorMessage("nPoints must be an integer.", nPoints)
  }
  
  # Ns 
  if (!is.na(names(Args["Ns"]))) {
    Ns <- eval(expression(Ns), parent.frame())
    if (!is.null(Ns)) { 
      if (!is.numeric(Ns))
        ErrorMessage("Ns must be numeric.", Ns)
      if (any(Ns < 0))
        ErrorMessage("All abundance values must be positive.", Ns)
    }
  }
  # Nexp 
  if (!is.na(names(Args["Nexp"]))) {
    Nexp <- eval(expression(Nexp), parent.frame())
    if (!is.null(Nexp)) {   
      if (!is.numeric(Nexp))
        ErrorMessage("Nexp must be numeric.", Nexp)
      if (any(Nexp < 0))
        ErrorMessage("All abundance values must be positive.", Nexp)
    }
  } 
  
  # NumberOfSimulations 
  if (!is.na(names(Args["NumberOfSimulations"]))) {
    NumberOfSimulations <- eval(expression(NumberOfSimulations), parent.frame())
    if (!is.numeric(NumberOfSimulations))
      ErrorMessage("NumberOfSimulations must be a number.", NumberOfSimulations)
     if (NumberOfSimulations < 0)
       ErrorMessage("NumberOfSimulations must be positive.", NumberOfSimulations)
  }
  
  # Ps 
  if (!is.na(names(Args["Ps"]))) {
    Ps <- eval(expression(Ps), parent.frame())
    if (!is.null(Ps)) {
      if (!is.numeric(Ps))
        ErrorMessage("Ps must be numeric.", Ps)    
      # Probabilities must sum to 1
      if (!isTRUE(all.equal(sum(Ps), 1)))
        ErrorMessage("Probabilities must sum to 1.", Ps)
      if (any(Ps < 0))
        ErrorMessage("All probabilities must be positive.", Ps)
    }
  }
  # Pexp 
  if (!is.na(names(Args["Pexp"]))) {
    Pexp <- eval(expression(Pexp), parent.frame())
    if (!is.null(Pexp)) {   
      if (!is.numeric(Pexp))
        ErrorMessage("Pexp must be numeric.", Pexp) 
      # Probabilities must sum to 1
      if (!isTRUE(all.equal(sum(Pexp), 1)))
        ErrorMessage("Probabilities must sum to 1.", Pexp) 
      if (any(Pexp < 0))
        ErrorMessage("All probabilities must be positive.", Pexp)
    }
  }

  # NorP 
  if (!is.na(names(Args["NorP"]))) {
    NorP <- eval(expression(NorP), parent.frame())
    if (!is.numeric(NorP))
      ErrorMessage("NorP must be numeric.", NorP)
    if (any(NorP < 0))
      ErrorMessage("All NorP values must be positive.", NorP)   
    if (!is.vector(NorP) & !is.SpeciesDistribution(NorP)) {
      # NorP may be a true vector or a SpeciesDistribution. Then dim(NorP) is NULL, and nothing more has to be checked
      # or a "named vector" whose attributes are not "names". Then dim() returns the vector's length.
      if (length(dim(NorP)) != 1) {
        # or a 2D numeric object
        if (length(dim(NorP)) == 2) {
          if (dim(NorP)[2] > 2) {
            # then it must have 1 or 2 columns
            ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)                      
          }
        } else {
          ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)          
        }
      } 
    }
  }
  # NorPexp 
  if (!is.na(names(Args["NorPexp"]))) {
    NorPexp <- eval(expression(NorPexp), parent.frame())
    if (!is.numeric(NorPexp))
      ErrorMessage("NorPexp must be numeric.", NorPexp)
    if (any(NorPexp < 0))
      ErrorMessage("All NorPexp values must be positive.", NorPexp)   
    if (!is.vector(NorPexp) & !is.SpeciesDistribution(NorPexp)) {
      # NorPexp may be a true vector or a SpeciesDistribution. Then dim(NorPexp) is NULL, and nothing more has to be checked
      # or a "named vector" whose attributes are not "names". Then dim() returns the vector's length.
      if (length(dim(NorPexp)) != 1) {
        # or a 2D numeric object
        if (length(dim(NorPexp)) == 2) {
          if (dim(NorPexp)[2] > 2) {
            # then it must have 1 or 2 columns
            ErrorMessage("NorPexp may be a vector or a two-column matrix.", NorPexp)                      
          }
        } else {
          ErrorMessage("NorPexp may be a vector or a two-column matrix.", NorPexp)          
        }
      } 
    }
  }
  
  # r
  if (!is.na(names(Args["r"]))) {
    r <- eval(expression(r), parent.frame())
    if (!is.numeric(r) | length(r)!=1)
      ErrorMessage("r must be a number.", r)
    if (r < 1)
      ErrorMessage("r must be at least 1.", r)
    if (as.integer(r) != r)
      ErrorMessage("r must be an integer.", r)
  }

    # sd
  if (!is.na(names(Args["sd"]))) {
    sd <- eval(expression(sd), parent.frame())
    if (!is.numeric(sd) | length(sd)!=1)
      ErrorMessage("sd must be a number.", sd)
    if (any(sd < 0))
      ErrorMessage("sd must be positive.", sd)
  }

  # Simulations 
  if (!is.na(names(Args["Simulations"]))) {
    Simulations <- eval(expression(Simulations), parent.frame())
    if (!is.numeric(Simulations))
      ErrorMessage("Simulations must be numeric.", Simulations)
    if (!is.vector(Simulations) | length(Simulations) > 1)
      ErrorMessage("Simulations must be a single number.", Simulations)
    if (any(Simulations < 1))
      ErrorMessage("Simulations must be at least 1.", Simulations)
    if (as.integer(Simulations) != Simulations)
      ErrorMessage("Simulations must be an integer.", Simulations)
  }

  # Tree 
  if (!is.na(names(Args["Tree"]))) {
    Tree <- eval(expression(Tree), parent.frame())
    if (!is.null(Tree)) {
      if (!inherits(Tree, "phylo") & !inherits(Tree, "phylog") & !inherits(Tree, "hclust") & !inherits(Tree, "PPtree"))
        ErrorMessage("Tree may be NULL or an object of class hclust or phylo or phylog or PPtree.", Tree)
      if (inherits(Tree, "phylog")) {
        if (is.null(Tree$Wdist))
          ErrorMessage("phylog Tree must contain a distance matrix (use add.tools=TRUE when creating it).", Tree)
      }
    }
  }
  # PhyloTree
  if (!is.na(names(Args["PhyloTree"]))) {
    PhyloTree <- eval(expression(PhyloTree), parent.frame())
    if (!is.null(PhyloTree)) {
      if (!inherits(PhyloTree, "phylo") & !inherits(PhyloTree, "phylog") & !inherits(PhyloTree, "hclust") & !inherits(PhyloTree, "PPtree"))
        ErrorMessage("PhyloTree may be NULL or an object of class hclust or phylo or phylog or PPtree", PhyloTree)
      if (inherits(PhyloTree, "phylog")) {
        if (is.null(PhyloTree$Wdist))
          ErrorMessage("phylog PhyloTree must contain a distance matrix (use add.tools=TRUE when creating it).", PhyloTree)
      }
    }
  }
  
  # prob
  if (!is.na(names(Args["prob"]))) {
    prob <- eval(expression(prob), parent.frame())
    if (!is.numeric(prob) | length(prob)!=1)
      ErrorMessage("prob must be a number.", prob)
    if (any(prob < 0) | any(prob > 1))
      ErrorMessage("prob must be between 0 and 1.", prob)
  }

  # RealValue
  if (!is.na(names(Args["RealValue"]))) {
    RealValue <- eval(expression(RealValue), parent.frame())
    if (!is.numeric(RealValue) | length(RealValue)!=1)
      ErrorMessage("RealValue must be a number.", RealValue)
  }

  # S
  if (!is.na(names(Args["S"]))) {
    S <- eval(expression(S), parent.frame())
    if (!is.numeric(S) | length(S)!=1)
      ErrorMessage("S must be a number.", S)
    if (any(S < 1))
      ErrorMessage("S must be at least 1.", S)
    if (as.integer(S) != S)
      ErrorMessage("S must be an integer.", S)
  }
  
  # SimulatedValues
  if (!is.na(names(Args["SimulatedValues"]))) {
    SimulatedValues <- eval(expression(SimulatedValues), parent.frame())
    if (!is.vector(SimulatedValues))
      ErrorMessage("SimulatedValues must be a numeric vector.", SimulatedValues)
  }

  # size
  if (!is.na(names(Args["size"]))) {
    size <- eval(expression(size), parent.frame())
    if (!is.numeric(size) | length(size)!=1)
      ErrorMessage("size must be a number.", size)
    if (any(size < 1))
      ErrorMessage("size must be at least 1.", size)
    if (as.integer(size) != size)
      ErrorMessage("size must be an integer.", size)
  }
  
  # Weights
  if (!is.na(names(Args["Weights"]))) {
    Weights <- eval(expression(Weights), parent.frame())
    if (!is.vector(Weights))
      ErrorMessage("Weights must be a numeric vector.", Weights)
    if (any(Weights < 0))
      ErrorMessage("All weights must be positive.", Weights)
  }

  # xy 
  if (!is.na(names(Args["xy"]))) {
    xy <- eval(expression(xy), parent.frame())
    if (!is.numeric(xy) | length(xy)!=2)
      ErrorMessage("xy must be a numeric vector of length 2.", xy)
  }

  # Z
  if (!is.na(names(Args["Z"]))) {
    Z <- eval(expression(Z), parent.frame())
    if (!is.null(Z)) {
      if (!is.matrix(Z)) {
        ErrorMessage("Z must be a square matrix.", Z)
      } else {
        if (dim(Z)[1] != dim(Z)[2])
          ErrorMessage("Z must be a square matrix.", Z)
        if (!is.null(colnames(Z)) | !is.null(rownames(Z))) {
          # If the matrix is named, rows and columns must have the same names
          if (!identical(colnames(Z), rownames(Z)))
            ErrorMessage("Z row and column names must be identical.", Z)
        }
        # Must be a relatedness matrix
        if (any(Z<0))
          ErrorMessage("All terms of the relatedness matrix Z must be positive.", Z)
        if (any(diag(Z)<0))
          ErrorMessage("All terms of the relatedness matrix Z diagonal must be strictly positive.", Z)
      }
    }  
  }
  
  return (TRUE)
}
