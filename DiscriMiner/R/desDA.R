#' Descriptive Discriminant Analysis
#' 
#' Performs a Descriptive Discriminant Analysis (a.k.a. Factorial Discriminant
#' Analysis from the french \emph{Analyse Factorielle Discriminante})
#' 
#' When \code{covar="within"} the estimated pooled within-class covariance
#' matrix is used in the calculations. \cr When \code{covar="total"} the total
#' covariance matrix is used in the calculations. \cr The difference between
#' \code{covar="within"} and \code{covar="total"} is in the obtained
#' eigenvalues.
#' 
#' The estiamted pooled within-class covariance matrix is actually the
#' within-class covariance matrix divided by the number of observations minus
#' the number of classes (see \code{\link{getWithin}})
#' 
#' @param variables matrix or data frame with explanatory variables
#' @param group vector or factor with group memberships
#' @param covar character string indicating the covariance matrix to be used.
#' Options are \code{"within"} and \code{"total"}
#' @return An object of class \code{"desda"}, basically a list with the
#' following elements
#' @return \item{power}{table with discriminant power of the 
#' explanatory variables}
#' @return \item{values}{table of eigenvalues}
#' @return \item{discrivar}{table of discriminant variables, 
#' i.e. the coefficients of the linear discriminant functions}
#' @return \item{discor}{table of correlations between the variables and the
#' discriminant axes}
#' @return \item{scores}{table of discriminant scores for each observation}
#' @author Gaston Sanchez
#' @seealso \code{\link{discPower}}
#' @references Lebart L., Piron M., Morineau A. (2006) \emph{Statistique
#' Exploratoire Multidimensionnelle}. Dunod, Paris.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load bordeaux wines dataset
#'   data(bordeaux)
#' 
#'   # descriptive discriminant analysis with within covariance matrix
#'   my_dda1 = desDA(bordeaux[,2:5], bordeaux$quality)
#'   my_dda1
#' 
#'   # descriptive discriminant analysis with total covariance matrix
#'   my_dda2 = desDA(bordeaux[,2:5], bordeaux$quality, covar="total")
#'   my_dda2
#'   
#'   # plot factor coordinates with ggplot
#'   library(ggplot2)
#'   bordeaux$f1 = my_dda1$scores[,1]
#'   bordeaux$f2 = my_dda1$scores[,2]
#'   ggplot(data=bordeaux, aes(x=f1, y=f2, colour=quality)) + 
#'   geom_hline(yintercept=0, colour="gray70") +
#'   geom_vline(xintercept=0, colour="gray70") +
#'   geom_text(aes(label=year), size=4) + 
#'   opts(title="Discriminant Map - Bordeaux Wines (years)")
#'  }
#' 
desDA <-
function(variables, group, covar="within")
{
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  check_cov <- covar %in% c("within", "total")
  if (!check_cov)
    warning("\nInvalid covar value; 'covar = within' is used")
  # get ingredients
  X = verify_Xy$X
  y = verify_Xy$y
  
  # how many obs and variables
  n = nrow(X)
  p = ncol(X)
  # group levels and number of levels
  glevs = levels(y)
  ng = nlevels(y)
  # how many obs in each group
  nk = as.vector(table(y))
  # number of factors
  nf = min(p, ng-1)
  # global mean
  gm = colMeans(X)
  
  # total covariance
  V = var(X)
  # between-class covariance matrix
  B = my_betweenCov(X, y)
  # estimated within-class covariance matrix
  W = my_withinCov(X, y)
  # within-class covariance matrix for disc-power
  Ww = matrix(0, p, p)
  for (k in 1:ng)
  {
    tmp <- y == glevs[k]
    nk = sum(tmp)
    Wk = ((nk-1)/(n-1)) * var(X[tmp,])
    Ww = Ww + Wk
  }
    
  ## Discriminant importance of explanatory variables
  # F-statistics
  F_statistic = ((n-ng)/(ng-1)) * (diag(B) / diag(Ww))
  p_values = 1 - pf(F_statistic, ng-1, n-ng)
  # Wilk's lambdas
  wilks_lamb = diag(Ww / V)
  # correlation ratio
  cor_ratio = diag(B) / diag(V)
  # table of disc power
  tab1 = cbind(cor_ratio, wilks_lamb, F_statistic, p_values)
  
  ## Discriminant axes
  # select covariance matrix
  if (covar == "within") Cov = W else Cov = ((n-1)/n) * var(X)
  # group means matrix
  GM = groupMeans(X, y)
  # decomposing between-class matrix:  B = CC'
  GM_gm = sweep(GM, 1, gm, FUN="-")
  C = sweep(GM_gm, 2, sqrt(nk/n), FUN="*")
  # eigen-decomposition
  EIG = eigen(t(C) %*% solve(Cov) %*% C)
  # eigenvalues
  lam = EIG$values[1:nf]
  # eigenvectors
  U = solve(V) %*% C %*% EIG$vectors[,1:nf]
  # normalizing eigenvectors
  aux = sqrt(diag(t(U) %*% Cov %*% U))
  Unorm = sweep(U, 2, aux, FUN="/")
  # table of eigenvalues
  tab2 = cbind(lam, 100*lam/sum(lam), 100*cumsum(lam)/sum(lam))
  colnames(tab2) = c("value", "proportion", "accumulated")
  rownames(tab2) = paste("DF", 1:nf, sep="")
  
  # Linear Discriminant Functions
  alphas = (-1) * gm %*% Unorm
  tab3 = rbind(alphas, Unorm)
  rownames(tab3) = c("constant", colnames(X))
  colnames(tab3) = paste("DF", 1:nf, sep="")
  
  # factor coordinates and correlations with expl variables
  Fs = X %*% Unorm + matrix(rep(alphas,each=n), n, nf)
  colnames(Fs) = paste("z", 1:nf, sep="")
  tab4 = cor(X, Fs)
  colnames(tab4) = paste("DF", 1:nf, sep="")
  
  # results
  structure(list(power = tab1, 
                 values = tab2, 
                 discrivar = tab3, 
                 discor = tab4, 
                 scores = Fs),
            class = "desda")
}
