#'@title PLS-CA: Partial Least Squares Canonical Analysis
#'
#'@description
#'Performs partial least squares canonical analysis for two blocks of data.
#'Compared to PLSR2, the blocks of variables in PLS-CA play a symmetric role
#'(i.e. there is neither predictors nor responses)
#'
#'@param X A numeric matrix or data frame (X-block) with more than one variable. 
#'No missing data are allowed
#'@param Y A numeric matrix or data frame (Y-block) with more than one variable. 
#'No missing data are allowed
#'@param comps The number of extracted PLS components (\code{NULL} by default)
#'When \code{comps=NULL} the number of components is determined by taking the
#'minimum between the number of columns from \code{X} and \code{Y}.
#'@param scaled A logical value indicating whether scaling data should be
#'performed (\code{TRUE} by default).
#'#'When \code{scaled=TRUE} the data is scaled to standardized values (mean=0,
#'variance=1). Otherwise the data will only be centered (mean=0).
#'@return An object of class \code{"plsca"}, basically a list with the
#'following elements:
#'@return \item{x.scores}{scores of the X-block (also known as T
#'components)}
#'@return \item{x.wgs}{weights of the X-block}
#'@return \item{x.loads}{loadings of the X-block}
#'@return \item{y.scores}{scores of the Y-block (also known as U
#'components)}
#'@return \item{y.wgs}{weights of the Y-block}
#'@return \item{y.loads}{loadings of the Y-block}
#'@return \item{cor.xt}{correlations between X and T}
#'@return \item{cor.yu}{correlations between Y and U}
#'@return \item{cor.tu}{correlations between T and U}
#'@return \item{cor.xu}{correlations between X and U}
#'@return \item{cor.yt}{correlations between Y and T}
#'@return \item{R2X}{explained variance of X by T}
#'@return \item{R2Y}{explained variance of Y by U}
#'@return \item{com.xu}{communality of X with U}
#'@return \item{com.yt}{communality of Y with T}
#'@author Gaston Sanchez
#'@seealso \code{\link{plot.plsca}}
#'@references Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et
#'Pratique}. Editions TECHNIP, Paris.
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## example of PLSCA with the vehicles dataset
#'  data(vehicles)
#'  
#'  # apply plsca
#'  my_plsca = plsca(vehicles[,1:12], vehicles[,13:16])
#'  my_plsca
#'  
#'  # plot variables
#'  plot(my_plsca)
#'  }
#'
plsca <-
function(X, Y, comps=NULL, scaled=TRUE)
{
  # ============ checking arguments ============
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(X) != nrow(Y))
    stop("\nDifferent number of rows in X and Y")
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y) 
  if (p == 1) stop("\nX must contain more than one variable")
  if (q == 1) stop("\nY must contain more than one variable")
  if (any(is.na(X)) | any(is.na(Y))) stop("\nSorry, No missing data are allowed")
  if (is.null(colnames(X)))
    colnames(X) = paste(rep("X",p), 1:p, sep="")
  if (is.null(rownames(X)))
    rownames(X) = rep(1:n)
  if (is.null(colnames(Y)))
    colnames(Y) = paste(rep("Y",q),1:q,sep="")
  if (is.null(rownames(Y)))
    rownames(Y) = rep(1:n)        
  if (!is.logical(scaled)) scaled<-TRUE
  if (!is.null(comps)) 
  {
    nc = comps
    if (mode(nc)!="numeric" || length(nc)!=1 || 
      nc<=1 || (nc%%1)!=0 || nc>min(n,p,q))
      nc = min(n, p, q)   
  } else nc = min(n, p, q)
  
  # ============ setting inputs ==============    
  if (scaled) X <- scale(X) else X <- scale(X, scale=F)
  if (scaled) Y <- scale(Y) else Y <- scale(Y, scale=F)
  X.old = scale(X)
  Y.old = scale(Y)
  Xt = matrix(NA, n, nc)
  Yu = matrix(NA, n, nc)
  Ah = matrix(NA, p, nc)
  Bh = matrix(NA, q, nc)
  Ch = matrix(NA, p, nc)
  Eh = matrix(NA, q, nc)
  
  # ============ pls CA algorithm ==============
  for (h in 1:nc)
  {
    th = X.old[,1]
    uh = Y.old[,1]
    ah.old = rep(1,p)
    iter = 1
    repeat
    {
      ah = t(X.old)%*%uh / sum(uh^2)
      ah = ah / sqrt(sum(ah^2))
      th = X.old %*% ah 
      bh = t(Y.old)%*%th / sum(th^2)
      bh = bh / sqrt(sum(bh^2))
      uh = Y.old %*% bh 
      ah.dif = sum((ah - ah.old)^2)
      if (ah.dif < 1e-06 || iter == 100) break
      iter = iter + 1
      ah.old = ah
    }        
    ch = t(X.old) %*% th / sum(th^2)
    eh = t(Y.old) %*% uh / sum(uh^2)
    X.old = X.old - th%*%t(ch)
    Y.old = Y.old - uh%*%t(eh)
    Xt[,h] = th
    Yu[,h] = uh
    Ah[,h] = ah
    Bh[,h] = bh
    Ch[,h] = ch
    Eh[,h] = eh
  }
  dimnames(Ch) = list(colnames(X), paste(rep("c",nc),1:nc,sep=""))
  dimnames(Eh) = list(colnames(Y), paste(rep("e",nc),1:nc,sep=""))
  dimnames(Xt) = list(rownames(X), paste(rep("t",nc),1:nc,sep=""))
  dimnames(Yu) = list(rownames(Y), paste(rep("u",nc),1:nc,sep=""))
  As = Ah %*% solve(t(Ch)%*%Ah)
  Bs = Bh %*% solve(t(Eh)%*%Bh)
  dimnames(As) = list(colnames(X), paste(rep("t",nc),1:nc,sep=""))
  dimnames(Bs) = list(colnames(Y), paste(rep("u",nc),1:nc,sep="")) 
  XT = scale(X) %*% As
  YU = scale(Y) %*% Bs
  dimnames(XT) = list(rownames(X), paste(rep("t",nc),1:nc,sep=""))
  dimnames(YU) = list(rownames(Y), paste(rep("u",nc),1:nc,sep=""))
  cor.xt = cor(X, Xt)
  cor.yu = cor(Y, Yu)
  cor.tu = cor(Xt, Yu)
  cor.xu = cor(X, Yu)
  cor.yt = cor(Y, Xt)
  R2X = t(apply(cor.xt^2, 1, cumsum))
  R2Y = t(apply(cor.yu^2, 1, cumsum))
  Com.XU = t(apply(cor.xu^2, 1, cumsum))
  Com.YT = t(apply(cor.yt^2, 1, cumsum))
  
  # ============ results ==============
  res = list(x.scores = XT, 
             x.wgs = As, 
             x.loads = Ch, 
             y.scores = YU, 
             y.wgs = Bs, 
             y.loads = Eh,
             cor.xt = cor.xt, 
             cor.yu = cor.yu, 
             cor.tu = cor.tu, 
             cor.xu = cor.xu, 
             cor.yt = cor.yt,
             R2X = R2X, 
             R2Y = R2Y, 
             com.xu = Com.XU, 
             com.yt = Com.YT)
  class(res) = "plsca"
  return(res)
}
