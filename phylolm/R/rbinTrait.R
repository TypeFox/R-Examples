rbinTrait <- function(n=1, phy, beta, alpha, X = NULL, model = c("LogReg"))
{
  model = match.arg(model)
  if (is.null(n)) stop("n needs to be an integer (number of replicates)")
  if (length(n)>1) stop("n needs to be an integer (number of replicates)")
  n = as.numeric(n)
  if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo".')
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")		
  phy = reorder(phy,"pruningwise")
  ntip <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- ntip + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  d = length(beta)	

  if ((d==1)&(!is.null(X)))
    stop("The design matrix is not needed when the coefficient is a scalar.")
  if (d>1) {
    if (is.null(X)) stop("there is no independent variables.")
    X = as.matrix(X)
    if (length(beta)!=ncol(X))
      stop("the number of columns in the design matrix does not
match the length of the vector of coefficients.")
    if (nrow(X)!=ntip)
      stop("the number of rows in the design matrix does not
match the number of tips in the tree.")
    if (is.null(rownames(X))) {
      warning("independent variables have no tip labels,
order assumed to be the same as in the tree.\n")
      data.names = phy$tip.label 
    } else data.names = rownames(X)
    order = match(data.names, phy$tip.label)
    if (sum(is.na(order))>0)
      stop("data names do not match with the tip labels.\n")
    g = X%*%beta 
    mu = as.vector(exp(g)/(1+exp(g)))
    p = mean(mu)
  } else {
    g = as.numeric(beta)
    p = exp(g)/(1+exp(g))		
  }
  q = 1-p
	
  y <- matrix(0, ntip + phy$Nnode, n) # simulated values at all nodes, n reps
  el = phy$edge.length	
  y[ROOT,] = as.numeric(runif(n)<p)
  for (i in N:1){ # Markov process along tree
    j0 = (y[anc[i],] == 0)
    y[des[i], j0] = as.numeric(runif(sum( j0))<(p-p*exp(-el[i]*alpha)))
    y[des[i],!j0] = as.numeric(runif(sum(!j0))<(p+q*exp(-el[i]*alpha)))
  }
  y <- y[1:ntip,]

  if (d>1) {     # influence of predictors at tips
    et = matrix(pmin(mu/p,(1-mu)/(1-p)), ntip,n) 
    b  = matrix(as.numeric(mu>p),        ntip,n)
    ij0 = (y==0)
    # next: flip some 0's to 1's if mu>p, and some 1's to 0's if mu<p
    y[ ij0] = as.numeric(runif(sum( ij0)) < (b[ ij0]-b[ ij0]*et[ ij0]))
    y[!ij0] = as.numeric(runif(sum(!ij0)) < (b[!ij0]+(1-b[!ij0])*et[!ij0]))
  }

  if (n==1) names(y) <- phy$tip.label
  else   rownames(y) <- phy$tip.label
  return(y)
}


