#'@title NIPALS: Non-linear Iterative Partial Least Squares
#'
#'@description
#'Principal Components Analysis with NIPALS algorithm
#'
#'@details
#'The function \code{nipals} performs Principal Components Analysis of a data
#'matrix that may contain missing values.
#'
#'@param Data A numeric matrix or data frame (which may contain
#'missing values).
#'@param comps Number of components to be calculated (by default 2)
#'@param scaled A logical value indicating whether to scale the data
#'(\code{TRUE} by default).
#'@return An object of class \code{"nipals"}, basically a list with the
#'following elements:
#'
#'When the analyzed data contain missing values, the help interpretation tools
#'(e.g. \code{cor.xt, disto, contrib, cos, dmod}) may not be meaningful, that
#'is to say, some of the results may not be coherent.
#'@return \item{values}{The pseudo eigenvalues}
#'@return \item{scores}{The extracted scores (i.e. components)}
#'@return \item{loadings}{The loadings}
#'@return \item{cor.xt}{Correlations between the variables and the scores}
#'@return \item{disto}{Squared distance of the observations to the origin}
#'@return \item{contrib}{Contributions of the observations (rows)}
#'@return \item{cos}{Squared cosinus}
#'@return \item{dmod}{Distance to the Model}
#'@author Gaston Sanchez
#'@seealso \code{\link{plot.nipals}}, \code{\link{plsreg1}}
#'@references Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et
#'Pratique.} Paris: Editions TECHNIP.
#'
#'Tenenhaus, M. (2007) \emph{Statistique. Methodes pour decrire, expliquer et
#'prevoir}. Paris: Dunod.
#'@export
#'@examples
#'
#'  \dontrun{
#'  # load datasets carscomplete and carsmissing
#'  data(carscomplete) # complete data
#'  data(carsmissing)  # missing values
#'
#'  # apply nipals
#'  my_nipals1 = nipals(carscomplete)
#'  my_nipals2 = nipals(carsmissing)
#'
#'  # plot variables (circle of correlations)
#'  plot(my_nipals1, what="variables", main="Complete data")
#'  plot(my_nipals2, what="variables", main="Missing data")
#'
#'  # plot observations with labels
#'  plot(my_nipals1, what="observations", show.names=TRUE, main="Complete data")
#'  plot(my_nipals2, what="observations", show.names=TRUE, main="Missing data")
#'  
#'  # compare results between my_nipals1 and my_nipals2
#'  plot(my_nipals1$scores[,1], my_nipals2$scores[,1], type="n")
#'  title("Scores comparison: my_nipals1  -vs-  my_nipals2", cex.main=0.9)
#'  abline(a=0, b=1, col="gray85", lwd=2)
#'  points(my_nipals1$scores[,1], my_nipals2$scores[,1], pch=21, 
#'         col="#5592e3", bg = "#5b9cf277", lwd=1.5)
#'  }
#'
nipals <-
function(Data, comps = 2, scaled = TRUE)
{
  # =======================================================
  # checking arguments
  # =======================================================
  X = as.matrix(Data)
  n = nrow(X)
  p = ncol(X)
  if (!n || !p) stop("dimension 0 in Data")
  if (p==1) stop("\nData must be a numeric matrix or data frame")
  if (is.null(colnames(X)))
    colnames(X) = paste(rep("X",p), 1:p, sep="")
  if (is.null(rownames(X)))
    rownames(X) = rep(1:n)
  if (!is.logical(scaled)) scaled = TRUE
  if (scaled) X <- scale(X) else X <- scale(X, scale=FALSE)    
  nc = comps
  if (mode(nc)!="numeric" || length(nc)!=1 || 
    nc<=1 || (nc%%1)!=0 || nc>min(n,p))
    nc <- 2
  if (nc==n) nc = n - 1
  if (any(is.na(X))) na.miss = TRUE else na.miss = FALSE
  
  # =======================================================
  # NIPALS algorithm
  # =======================================================
  # Prepare ingredients
  X.old = X
  Th = matrix(NA, n, nc) # scores
  Ph = matrix(NA, p, nc) # loadings
  eigvals = rep(NA, nc) # eigenvalues
  # iterative process
  for (h in 1:nc)
  {
    th.new = X.old[,1]
    ph.old = rep(1, p)
    ph.new = rep(1, p)
    iter = 1
    repeat
    {
      # missing data
      if (na.miss)
      {
        for (j in 1:p)
        {
          i.exist = intersect(which(complete.cases(X[,j])), which(complete.cases(th.new)))
          ph.new[j] = sum(X.old[i.exist,j] * th.new[i.exist]) / sum(th.new[i.exist]^2)
        }
        ph.new = ph.new / sqrt(sum(ph.new^2))       
        for (i in 1:n)
        {
          j.exist = which(complete.cases(X[i,]))
          th.new[i] = sum(X.old[i,j.exist] * ph.new[j.exist]) / sum(ph.new[j.exist]^2)
        }
      } else {
        # no missing data
        ph.new = t(X.old) %*% th.new / sum(th.new^2)
        ph.new = ph.new / sqrt(sum(ph.new^2))
        th.new = X.old %*% ph.new
      }
      ph.aux = sum((ph.new - ph.old)^2)
      ph.old = ph.new
      # check convergence
      if (ph.aux < 1e-06 || iter==100) break
      iter = iter + 1
    }
    Th[,h] = th.new
    Ph[,h] = ph.new
    X.new = X.old - th.new %*% t(ph.new)
    X.old = X.new
    eigvals[h] = sum(th.new^2) / (n-1)
  }
  dimnames(Th) = list(rownames(X), paste(rep("t",nc), 1:nc, sep=""))
  dimnames(Ph) = list(colnames(X), paste(rep("p",nc), 1:nc, sep=""))
  
  # =======================================================
  # eigenvalues
  # =======================================================
  # eigenvalues
  eig.perc = 100 * eigvals / p
  eigs = data.frame(values=eigvals, percentage=eig.perc, cumulative=cumsum(eig.perc))
  rownames(eigs) = paste(rep("v",nc), 1:nc, sep="")

  # =======================================================
  # interpretation tools
  # =======================================================
  # correlation between components and variables
  if (na.miss) {
    cor.sco = matrix(NA, p, nc)
    for (j in 1:p) {
      i.exist = which(complete.cases(X[,j]))
      cor.sco[j,] = round(cor(X[i.exist,j], Th[i.exist,]), 4)
    }
    dimnames(cor.sco) = list(colnames(X), colnames(Th))
  } else
    cor.sco = cor(X, Th)
  # individuals contribution
  ConInd = (100/n) * Th^2 %*% diag(1/eigvals)
  dimnames(ConInd) = list(rownames(X), paste(rep("ctr",nc),1:nc,sep=".") )
  # square distance of projections
  D.proj = Th^2
  D.orig = rowSums(X^2, na.rm=TRUE) # square distance of centered data
  # square cosinus
  Cos.inds = matrix(0, n, nc)
  for (i in 1:n)
    Cos.inds[i,] = D.proj[i,] / D.orig[i] 
  dimnames(Cos.inds) = list(rownames(X), paste(rep("cos2",nc),1:nc,sep=".") )
  # Distance to the model: DModX 
  DMX2 = matrix(NA, n, nc)
  DMX2[,1] = (D.orig * (1-Cos.inds[,1])) / (p-eigvals[1])
  for (j in 2:nc)
    DMX2[,j] = (D.orig * (1 - rowSums(Cos.inds[,1:j]))) / (p-sum(eigvals[1:j]))
  # Distance to the Model X (normalized)
  DMX = DMX2 
  dimnames(DMX) = list(rownames(X), colnames(Th))

  # results 
  res = list(values = eigs, 
              scores = Th, 
              loadings = Ph, 
              cor.xt = cor.sco, 
              disto = D.orig, 
              contrib = ConInd, 
              cos = Cos.inds, 
              dmod = DMX)
  class(res) = "nipals"
  return(res)
}

