#' Random Skewers projection
#' 
#' Not tested!
#' Uses MCMC Bayeisian posterior samples of a set of covariance matrices to identify 
#' directions of the morphospace in which these matrices differ in their amount of genetic variance.
#' 
#' @param cov.matrix.array Array with dimentions traits x traits x populations x MCMCsamples
#' @param p significance treashhold for comparison of variation in each random direction
#' @param num.vectors number of random vectors
#' @return projection of all matrices in all random vectors
#' @return set of random vectors and confidence intervals for the projections
#' @return eigen decomposition of the random vectors in directions with significant differences of variations
#' #@export
#' @importFrom coda HPDinterval as.mcmc
#' @references Aguirre, J. D., E. Hine, K. McGuigan, and M. W. Blows. "Comparing G: multivariate analysis of genetic variation in multiple populations." Heredity 112, no. 1 (2014): 21-29.
#' @examples
#' #random set of covariance matrices 
#' cov.matrices = aperm(aaply(1:15, 1, function(x) 
#'                      laply(RandomMatrix(10, 100, 
#'                            variance = runif(10, 1, 10)), 
#'                            identity)), 
#'                      c(3, 4, 1, 2))
#' #rs_proj = evolqg:::RSProjection(cov.matrices, p = 0.8)  
#' #plot(rs.proj, cov.matrices)
RSProjection <- function(cov.matrix.array, p = 0.95, num.vectors = 1000){
  if (dim(cov.matrix.array)[[1]] != dim(cov.matrix.array)[[2]]){
    stop("Covariance matrix array must be of order n x n x m x MCMCsamplele")
  }
  if (is.na(dim(cov.matrix.array)[4])) {
    stop("There are no MCMCsampleles")
  }
  n <- dim(cov.matrix.array)[[1]]
  m <- dim(cov.matrix.array)[[3]]
  MCMCsample <- dim(cov.matrix.array)[[4]]
  
  rand.vec = aaply(matrix(rnorm(n*num.vectors), n, num.vectors), 2, Normalize)

  # Projections of all matrices on the random vectors
  cov.projection <- apply(cov.matrix.array, 3:4, function(mat) diag(rand.vec %*% mat %*% t(rand.vec)))
  cov.projection <- aperm(cov.projection, c(3, 2, 1))
  colnames(cov.projection) <- dimnames(cov.matrix.array)[[3]]

  #setting up an index for HPD comparisons
  prs <- cbind(rep(1:m, each = m), 1:m)
  prs.comp <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  
  proj.score <-matrix(NA, num.vectors, ((m^2 - m)/2))
  
  #for a given random vector, examine if the HPDintervals of any pair of  G matrices overlap
  HPD.int <- alply(cov.projection, 3, function(proj) HPDinterval(as.mcmc(proj), prob = p))
  proj.score = t(sapply(HPD.int, function(HPD) ifelse(HPD[prs.comp[,1], 1] > HPD[prs.comp[,2],2] | 
                                                        HPD[prs.comp[,2], 1] > HPD[prs.comp[,1],2],
                                                      1, 0)))
  
  #collate the random vectors and the outcome of their projection on the G matrices
  vec.score <- cbind(rand.vec, proj.score)
  colnames(vec.score) <- c(1:n, 
                           paste0(dimnames(cov.matrix.array)[[3]][prs.comp[, 1]], 
                                  ".vs.", 
                                  dimnames(cov.matrix.array)[[3]][prs.comp[, 2]]))
  
  sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0)

  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger num.vectors or lower p"); 
    eig.R <- NA}
  else{
    #eigen decomposition of the matrix of significant random vectors
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(cov.matrix.array)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }
  
  out = list(cov.projection = cov.projection, vec.score = vec.score, eig.R = eig.R)
  class(out) = "rs_proj"
  return(out)
}
# 
# plot.rs_proj <- function (rs_proj, cov.matrix.array){
#   
#   n <- length(rs_proj$eig.R$values)
# 
#   HPD.int <- aaply(rs_proj$cov.projection, 3, function(proj) HPDinterval(as.mcmc(proj), prob = 0.95))
#   HPD.int <- aperm(HPD.int, c(2, 3, 1))
#   dimnames(HPD.int) = list(colnames(rs_proj$cov.projection), NULL, NULL)
#   dat = adply(HPD.int, 1:3)
#   
#   names(dat) = c('Population', 'interval', 'trait', 'value')
#   dat = dcast(dat, Population+trait~interval)
#   names(dat) = c('Population', 'trait', 'lower', 'upper')
#   dat$trait = paste0("PC", 1:n)
#   order.list = paste0("PC", 1:n) 
#   dat$trait = factor(dat$trait, levels = order.list)
#   dat$mean = rowMeans(cbind(dat$upper, dat$lower))
#   plot = ggplot(dat, aes(Population, mean)) + geom_point() +
#     geom_errorbar(aes( ymin=lower, ymax=upper)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ylab('Genetic Variance') + facet_wrap(~ trait, scales = "free_y")
#   return(plot)
# }