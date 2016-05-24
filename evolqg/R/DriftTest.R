#'Test drift hypothesis
#'
#'Given a set of covariance matrices and means for terminals, test the hypothesis
#'that obseved divergency is larger/smaller than expected by drift alone using a regression of
#' the between-group varicances on the within-group eigenvalues.
#'
#'@param means list or array of species means being compared. array must have means in the rows.
#'@param cov.matrix ancestral covariance matrix for all populations
#'@param show.plot boolean. If TRUE, plot of eigenvalues of ancetral matrix by between group variance is showed.
#'@return list of results containing:
#'@return regresion: the linear regression between the log of the eigenvalues of the ancestral matrix and the log of the between group variance (projected on the eigenvectors of the ancenstral matrix)
#'@return coefficient_CI_95: confidence intervals for the regression coefficients
#'@return log.between_group_variance: log of the between group variance (projected on the ancestral matrix eigenvectors)
#'@return log.W_eVals: log of the ancestral matrix eigenvalues
#'@return plot: plot of the regression using ggplot2
#'@note If the regression coefficient is significantly different to one, the null hypothesis of drift is rejected.
#'@references Marroig, G., and Cheverud, J. M. (2004). Did natural selection or genetic drift 
#'produce the cranial diversification of neotropical monkeys? The American Naturalist, 163(3), 417-428. doi:10.1086/381693
#'@references Proa, M., O'Higgins, P. and Monteiro, L. R. (2013), Type I error rates for testing genetic drift with phenotypic covariance matrices: A simulation study. Evolution, 67: 185-195. doi: 10.1111/j.1558-5646.2012.01746.x
#'@author Ana Paula Assis, Diogo Melo
#'@export
#'@import plyr
#'@importFrom ggplot2 ggplot geom_text geom_smooth labs theme_bw aes_string
#'@importFrom stats na.omit
#'@examples
#'
#' #Input can be an array with means in each row or a list of mean vectors
#'means = array(rnorm(40*10), c(10, 40)) 
#'cov.matrix = RandomMatrix(40, 1, 1, 10)
#'DriftTest(means, cov.matrix)
DriftTest <- function(means, cov.matrix, show.plot=TRUE)
{
  if(is.data.frame(means) | (!is.array(means) & !is.list(means)))
    stop("means must be in a list or an array.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(is.list(means)){
    mean.array <- laply(means, identity)
  }  else{ 
    mean.array <- means
  }
  W.pc <- eigen(cov.matrix)
  projection.Wpc <- as.matrix(mean.array) %*% W.pc$vectors #projecting the means in 
                                                           #the principal components of W
  log.B_variance <- log(apply(projection.Wpc, 2, var)) #variance between groups
  log.W_eVals <- log(W.pc$values) 
  regression <- lm(log.B_variance~log.W_eVals)
  reg.plot <- ggplot(data.frame(log.B_variance, 
                                log.W_eVals, 
                                names = 1:(dim(mean.array)[2])), 
                     aes_string('log.W_eVals', 'log.B_variance')) + 
              geom_text(aes_string(label = 'names')) + 
              geom_smooth(method = "lm", color = 'black') + 
              labs(x = "log(W Eigenvalues)", y = "log(B variances)") + theme_bw()
  if(show.plot) print(reg.plot)
  containsOne <- function(x) ifelse(x[1] < 1 & x[2] > 1, TRUE, FALSE)
  test <- !containsOne(confint(regression)[2,])
  names(test) <- '5 %'
  objeto <- list("regression" = regression,
                 "coefficient_CI_95" = confint(regression),
                 "log.between_group_variance" = log.B_variance, 
                 "log.W_eVals" = log.W_eVals,
                 "drift_rejected" = test,
                 "plot" = reg.plot)
  return(objeto)
}

#'Drift test along phylogeny
#'
#'Performs a regression drift test along a phylogeny using DriftTest function.
#'
#'@param tree phylogenetic tree
#'@param mean.list list of tip node means. Names must match tip node labels.
#'@param cov.matrix.list list of tip node covariance matrices. Names must match tip node labels.
#'@param sample.sizes vector of tip nodes sample sizes
#'@return A list of regression drift tests performed in nodes with over 4 descendant tips.
#'@export
#'@seealso DriftTest PlotTreeDriftTest
#'@author Diogo Melo
#'@examples
#'library(ape)
#'data(bird.orders)
#'
#'tree <- bird.orders
#'mean.list <- llply(tree$tip.label, function(x) rnorm(5))
#'names(mean.list) <- tree$tip.label
#'cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
#'names(cov.matrix.list) <- tree$tip.label
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'
#'test.list <- TreeDriftTest(tree, mean.list, cov.matrix.list, sample.sizes)
#'
#'#Ancestral node plot:
#'test.list[[length(test.list)]]$plot
TreeDriftTest <- function(tree, mean.list, cov.matrix.list, sample.sizes = NULL){
  if(!all(tree$tip.label %in% names(mean.list))) stop("All tip labels must be in names(mean.list).")
  if(!all(tree$tip.label %in% names(cov.matrix.list))) stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  if(!any(node.mask)) stop("At least one node must have more than 4 descendents in mean.list")
  test.list <- llply(nodes[node.mask], function(node) DriftTest(getMeans(mean.list, tree, node), 
                                                                cov.matrices[[node]], FALSE))
  names(test.list) <- nodes[node.mask]
  return(test.list)
}

#'Plot results from TreeDriftTest
#'
#'Plot which labels reject drift hypothesis.
#'
#'@param test.list Output from TreeDriftTest
#'@param tree phylogenetic tree
#'@param ... adition arguments to plot
#'@seealso DriftTest TreeDriftTest
#'@importFrom ape nodelabels
#'@author Diogo Melo
#'@export
#'@examples
#'library(ape)
#'data(bird.orders)
#'
#'tree <- bird.orders
#'mean.list <- llply(tree$tip.label, function(x) rnorm(5))
#'names(mean.list) <- tree$tip.label
#'cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
#'names(cov.matrix.list) <- tree$tip.label
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'
#'test.list <- TreeDriftTest(tree, mean.list, cov.matrix.list, sample.sizes)
#'PlotTreeDriftTest(test.list, tree)
PlotTreeDriftTest <- function(test.list, tree, ...){
  tested.nodes <- as.numeric(names(test.list))
  non.drift.nodes <- laply(test.list, function(x) x$drift_rejected)
  plot(tree, ...)
  nodelabels(node = tested.nodes, thermo = as.numeric(non.drift.nodes))  
}

#'@importFrom phytools getDescendants
getMeans <- function(mean.list, tree, node){
  means <- mean.list[na.omit(tree$tip.label[getDescendants(tree, node)])]
}
