#'PC Score Correlation Test
#'
#'Given a set of covariance matrices and means for terminals, test the hypothesis
#'that obseved divergency is larger/smaller than expected by drift alone using the correlation on
#'principal component scores.
#'
#'@param means list or array of species means being compared. array must have means in the rows.
#'@param cov.matrix ancestral covariance matrix for all populations
#'@param taxons names of taxons being compared. Must be in the same order of the means.
#'@param show.plots boolean. If TRUE, plot of eigenvalues of ancetral matrix by between group variance is showed.
#'@return list of results containing:
#'@return correlation matrix of principal component scores and p.values for each correlation. Lower triangle of outputput are correlations, and upper triangle are p.values.
#'@return if show.plots is TRUE, also returns a list of plots of all projections of the nth PCs, where n is the number of taxons.
#'@references Marroig, G., and Cheverud, J. M. (2004). Did natural selection or genetic drift 
#'produce the cranial diversification of neotropical monkeys? The American Naturalist, 163(3), 417-428. doi:10.1086/381693
#'@author Ana Paula Assis, Diogo Melo
#'@export
#'@import plyr
#'@importFrom ggplot2 ggplot geom_text geom_smooth labs theme_bw aes_string
#'@examples
#' #Input can be an array with means in each row or a list of mean vectors
#'means = array(rnorm(40*10), c(10, 40)) 
#'cov.matrix = RandomMatrix(40, 1, 1, 10)
#'taxons = LETTERS[1:10]
#'PCScoreCorrelation(means, cov.matrix, taxons)
#'
#'##Plots list can be displayed using grid.arrange()
#'#library(gridExtra)
#'#pc.score.output <- PCScoreCorrelation(means, cov.matrix, taxons, TRUE)
#'#do.call(grid.arrange, c(pc.score.output$plots,list(nrow=4,ncol=6)))
#'##Or we can print to file:
#'#ggsave("multipage.pdf", do.call(marrangeGrob, c(pc.score.output$plots, list(nrow=2, ncol=2))))
PCScoreCorrelation <- function(means, cov.matrix, taxons = names(means), show.plots = FALSE){
  if(is.data.frame(means) | (!is.array(means) & !is.list(means)))
    stop("means must be in a list or an array.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(is.list(means)){
    mean.array <- laply(means, identity)
  }  else{ 
    mean.array <- means
    if(is.null(taxons)) taxons <- rownames(means)
  }
  if(is.null(taxons)) taxons <- paste0("taxon", 1:nrow(mean.array))
  n.taxon <- length(taxons)-1
  if(n.taxon > dim(cov.matrix)[1]) n.taxon <- dim(cov.matrix)[1]
  
  W.pc <- eigen(cov.matrix)
  proj.med <- as.matrix(mean.array) %*% W.pc$vectors
  colnames(proj.med) <- paste0("PC", 1:ncol(proj.med))
  
  mat.pcs.deriva <- matrix(1 , ncol = n.taxon, nrow = n.taxon)
  if(show.plots){ 
    plots <- vector("list", (n.taxon^2 - n.taxon)/2)
    current.plot <- 1
  }
  for(i in 1:n.taxon){
    for (j in 1:i){
      if(j != i){
        test = cor.test(proj.med[,i], proj.med[,j])
        mat.pcs.deriva[i,j] <- test$estimate
        mat.pcs.deriva [j,i] <- test$p.value
        if(show.plots){
          plots[[current.plot]] <- ggplot(data.frame(x = proj.med[,i], 
                                                     y = proj.med[,j], 
                                                     taxons = taxons),
                                          aes_string('x', 'y')) +
                                   geom_text(aes_string(label = 'taxons')) + 
                                   geom_smooth(method = "lm", color = 'black') + 
                                   labs(x = paste0("PC", i), y = paste0("PC", j)) + 
                                   theme_bw()
          current.plot <- current.plot + 1
        }
      }
    }
  }
  if(show.plots) return(list("correlation_p.value" = mat.pcs.deriva, "plots" = plots))
  else return(mat.pcs.deriva) ## a diagonal inferior sao as correlacoes e a superior o valor de p
}
