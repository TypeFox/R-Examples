rTrait <- function(n=1, phy, model=c("BM","OU","lambda","kappa","delta","EB","trend"),
	parameters=NULL,plot.tree=FALSE)
{
  ## initialize
  if (is.null(n)) stop("n needs to be an integer (number of replicates)")
  if (length(n)>1) stop("n needs to be an integer (number of replicates)")
  n = as.numeric(n)
  if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo".')
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
  phy = reorder(phy,"pruningwise")
  model = match.arg(model)
  if ((model=="trend")&(is.ultrametric(phy)))
   stop("the trend is unidentifiable for ultrametric trees.")
  ntip <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- ntip + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]

  ## Default parameters
  parameters.default = c(0,1,0,0,1,1,1,0,0)
  names(parameters.default) = c("ancestral.state", "sigma2", 
   "optimal.value", "alpha", "lambda", "kappa", "delta", "rate", "trend")

  ## User defined parameters
  if (is.null(parameters)) { parameters = parameters.default } else { 
   if (class(parameters)!= "list") {
    stop("please specify parameters as a list().")
   } else {
    specified <- !c(is.null(parameters$ancestral.state), 
                    is.null(parameters$sigma2), 
                    is.null(parameters$optimal.value),
                    is.null(parameters$alpha),
                    is.null(parameters$lambda),
                    is.null(parameters$kappa),
                    is.null(parameters$delta),
                    is.null(parameters$rate),
                    is.null(parameters$trend))
    parameters.user <- c(parameters$ancestral.state,
                    parameters$sigma2,
                    parameters$optimal.value,
                    parameters$alpha,
                    parameters$lambda,
                    parameters$kappa,
                    parameters$delta,
                    parameters$rate,
                    parameters$trend)
    names(parameters.default) = c("ancestral.state","sigma2","optimal.value",
                                  "alpha", "lambda", "kappa", "delta", "rate", "trend")
    parameters <- parameters.default
    parameters[specified] <- parameters.user 
   }				
  }
  p = list(ancestral.state = parameters[1],
           sigma2 = parameters[2],
           optimal.value = parameters[3],
           alpha = parameters[4],
           lambda = parameters[5],
           kappa = parameters[6],
           delta = parameters[7],
           rate = parameters[8],
           trend = parameters[9])

  if (model == "OU" & p$alpha == 0)
    model = "BM"
  x <- matrix(0, ntip + phy$Nnode, n) # simulated values at all nodes, n reps
  if (model == "OU"){ # forward simulation along each edge, OU model
    gam = p$sigma2/(2*p$alpha)
    el = phy$edge.length
    eps <- matrix(rnorm(N*n, 0, sqrt(rep(gam,each=N)*(1-exp(-2*el*rep(p$alpha,each=N))))),
                  N,n) # increments along N edges, n replicates
    x[ROOT,] = p$ancestral.state # the root edge is ignored
    for (i in N:1)
      x[des[i],] = x[anc[i],]*exp(-p$alpha*el[i]) +
        p$optimal.value*(1-exp(-p$alpha*el[i])) + eps[i,]
  }
  else { # BM model on tree with transformed branch lengths + possible trend
    tree = transf.branch.lengths(phy,model,parameters=p,check.pruningwise=F)$tree
    if (plot.tree) { plot(tree,no.margin=T); add.scale.bar()}
    if (model != "trend") p$trend=0;
    eps <- matrix(rnorm(N*n, tree$edge.length*rep(p$trend,each=N),
                        sqrt(tree$edge.length*rep(p$sigma2,each=N))),
                  N,n) # increments along N edges, n replicates
    x[ROOT,] <- rnorm(n, p$ancestral.state + tree$root.edge*p$trend,
                      sqrt(tree$root.edge*p$sigma2))
    for (i in N:1) x[des[i],] = x[anc[i],] + eps[i,]
  }
  sim <- x[1:ntip,]
  if (n==1) names(sim) <- phy$tip.label
  else   rownames(sim) <- phy$tip.label
  return(sim)
}
