#' Remove collinearity among variables of a a raster stack
#' 
#' This functions analyses the correlation among variables of the provideded
#' stack of environmental variables (using Pearson's R), and can return a 
#' vector containing names of variables that are not intercorrelated, or a list
#' containing grouping variables according to their degree of collinearity.
#' 
#' @param raster.stack a RasterStack object, in which each layer represent an environmental 
#' variable.
#' @param multicollinearity.cutoff a numeric value corresponding to the cutoff
#' of correlation above which to group variables.
#' @param select.variables \code{TRUE} or \code{FALSE}. If \code{TRUE}, then the
#' function will choose one variable among each group to return a vector of
#' non correlated variables (see details). If \code{FALSE}, the function will return a list
#' containing the groups of correlated variables.
#' @param sample.points \code{TRUE} or \code{FALSE}. If you have a large
#' raster file then use this parameter to sample a number of points equal to
#' \code{nb.points}.
#' @param nb.points a numeric value. Only useful if \code{sample.points = TRUE}.
#' The number of sampled points from the raster, to perform the PCA. A too small
#' value may not be representative of the environmental conditions in your raster.
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the hierarchical
#' ascendant classification used to group variables will be plotted.
#' @return
#' a vector of non correlated variables, or a list where each element is a
#' group of non correlated variables.
#' @details
#' This function uses the Pearson's correlation coefficient to analyse 
#' correlation among variables. This coefficient is then used to compute a
#' distance matrix, which in turn is used it compute an ascendant hierarchical
#' classification, with the '\emph{complete}' method (see 
#' \code{\link[stats]{hclust}}). If at least one correlation above the \code{
#' multicollinearity.cutoff} is detected, then the variables will be grouped
#' according to their degree of correlation. 
#' 
#' If \code{select.variables = TRUE}, then the function will return a vector
#' containing variables that are not intercorrelated.
#' The variables not correlated to any other variables are automatically included
#' in this vector. For each group of intercorrelated variables, one variable will
#' be randomly chosen and included in this vector.
#' 
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @examples
#' # Create an example stack with six environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100),
#'              raster(a * logisticFun(1:100, alpha = 10, beta = 70)),
#'              raster(t(a)),
#'              raster(exp(a)),
#'              raster(log(a)))
#' names(env) <- paste("Var", 1:6, sep = "")   
#'    
#' # Defaults settings: cutoff at 0.7
#' removeCollinearity(env, plot = TRUE)
#' 
#' # Changing cutoff to 0.5
#' removeCollinearity(env, plot = TRUE, multicollinearity.cutoff = 0.5)
#' 
#' # Automatic selection of variables not intercorrelated
#' removeCollinearity(env, plot = TRUE, select.variables = TRUE)
#' 
#' # Assuming a very large raster file: selecting a subset of points
#' removeCollinearity(env, plot = TRUE, select.variables = TRUE,
#'                    sample.points = TRUE, nb.points = 5000)
#' 
#' 

removeCollinearity <- function(raster.stack, multicollinearity.cutoff = .7,
                               select.variables = FALSE, sample.points = FALSE, 
                               nb.points = 10000, plot = FALSE)
{
  if(sample.points)
  {
    if(!is.numeric(nb.points))
    {stop("nb.points must be a numeric value corresponding to the number of pixels to sample from raster.stack")}
    env.df <- sampleRandom(raster.stack, size = nb.points, na.rm = TRUE)
  } else
  {
    env.df <- getValues(raster.stack)
    if(any(is.na(env.df)))
    {
      env.df <- env.df[-unique(which(is.na(env.df), arr.ind = T)[, 1]), ] # Removing NAs 
    }
  }
  
  if(!is.numeric(multicollinearity.cutoff))
  {stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")} else
    if(multicollinearity.cutoff > 1 | multicollinearity.cutoff < 0)
    {stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")} 
  
  # Correlation matrix creation
  cor.matrix <- matrix(data = 0,
                       nrow = nlayers(raster.stack),
                       ncol = nlayers(raster.stack),
                       dimnames = list(names(raster.stack), names(raster.stack)))
  
  # Correlation based on Pearson
  cor.matrix <- 1 - abs(cor(env.df, method = "pearson"))
  
  # Transforming the correlation matrix into an ascendent hierarchical classification
  dist.matrix <- as.dist(cor.matrix)
  ahc <- stats::hclust(dist.matrix, method = "complete")
  groups <- stats::cutree(ahc, h = 1 - multicollinearity.cutoff)
  if(length(groups) == max(groups))
  {
    message(paste("  - No multicollinearity detected in your data at threshold ", multicollinearity.cutoff, "\n", sep = ""))
    mc <- FALSE
  } else
  { mc <- TRUE }
  
  if(plot)
  {
    par(mar = c(5.1, 5.1, 4.1, 3.1))
    plot(ahc, hang = -1, xlab = "", ylab = "Distance (1 - Pearson's r)", 
         main = "", las = 1,
         sub = "", axes = F)
    axis(2, at = seq(0, 1, length = 6), las = 1)
    if(mc)
    {
      title(paste('Groups of intercorrelated variables at cutoff', 
                  multicollinearity.cutoff))
      par(xpd = T)
      rect.hclust(ahc, h = 1 - multicollinearity.cutoff)
    } else
    {
      title(paste('No intercorrelation among variables at cutoff', 
                  multicollinearity.cutoff))
    }

  }
  
  # Random selection of variables
  if(select.variables)
  {
    sel.vars <- NULL
    for (i in 1:max(groups))
    {
      sel.vars <- c(sel.vars, sample(names(groups[groups == i]), 1))
    }
  } else
  {
    if(mc)
    {
      sel.vars <- list()
      for (i in groups)
      {
        sel.vars[[i]] <- names(groups)[groups == i]
      }
    } else
    {
      sel.vars <- names(raster.stack)
    }
  }
  return(sel.vars)
}
