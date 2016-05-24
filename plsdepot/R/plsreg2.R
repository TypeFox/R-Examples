#'@title PLS-R2: Partial Least Squares Regression 2
#'
#'@description
#'The function plsreg2 performs partial least squares regression for the multivariate case (i.e.
#'more than one response variable)
#'
#'@details
#'The minimum number of PLS components \code{comps} to be extracted is 2.
#'
#'The data is scaled to standardized values (mean=0, variance=1).
#'
#'The argument \code{crosval} gives the option to perform cross-validation.
#'This parameter takes into account how \code{comps} is specified.
#'When \code{comps=NULL}, the number of components is obtained by cross-validation.
#'When a number of components is specified, cross-validation results are calculated
#'for each component.
#'
#'@param predictors A numeric matrix or data frame containing the predictor variables.
#'@param responses A numeric matrix or data frame containing the response variables.
#'@param comps The number of extracted PLS components (2 by default)
#'@param crosval Logical indicating whether cross-validation should be performed
#'(\code{TRUE} by default). No cross-validation is done if there is missing data or
#'if there are less than 10 observations.
#'@return An object of class \code{"plsreg2"}, basically a list with the
#'following elements:
#'@return \item{x.scores}{components of the predictor variables (also known as T-components)}
#'@return \item{x.loads}{loadings of the predictor variables}
#'@return \item{y.scores}{components of the response variables (also known as U-components)}
#'@return \item{y.loads}{loadings of the response variables}
#'@return \item{cor.xt}{correlations between X and T}
#'@return \item{cor.yt}{correlations between Y and T}
#'@return \item{cor.xu}{correlations between X and U}
#'@return \item{cor.yu}{correlations between Y and U}
#'@return \item{cor.tu}{correlations between T and U}
#'@return \item{raw.wgs}{weights to calculate the PLS scores with the deflated
#'matrices of predictor variables}
#'@return \item{mod.wgs}{modified weights to calculate the PLS scores with the
#'matrix of predictor variables}
#'@return \item{std.coefs}{Vector of standardized regression coefficients (used with 
#'scaled data)}
#'@return \item{reg.coefs}{Vector of regression coefficients (used with the original
#'data)}
#'@return \item{y.pred}{Vector of predicted values}
#'@return \item{resid}{Vector of residuals}
#'@return \item{expvar}{table with R-squared coefficients}
#'@return \item{VIP}{Variable Importance for Projection}
#'@return \item{Q2}{table of Q2 indexes (i.e. leave-one-out cross validation)}
#'@return \item{Q2cum}{table of cummulated Q2 indexes}
#'@author Gaston Sanchez
#'@seealso \code{\link{plot.plsreg2}},
#'\code{\link{plsreg1}}.
#'@references Geladi, P., and Kowlaski, B. (1986) Partial Least Squares
#'Regression: A Tutorial. \emph{Analytica Chimica Acta}, \bold{185}, pp. 1-17.
#'
#'Hoskuldsson, A. (1988) PLS Regression Methods. \emph{Journal of
#'Chemometrics}, \bold{2}, pp. 211-228.
#'
#'Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et Pratique}. Editions
#'TECHNIP, Paris.
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## example of PLSR2 with the vehicles dataset
#'  data(vehicles)
#'  
#'  # apply plsreg2 extracting 2 components (no cross-validation)
#'  pls2_one = plsreg2(vehicles[,1:12], vehicles[,13:16], comps=2, crosval=FALSE)
#'  
#'  # apply plsreg2 with selection of components by cross-validation
#'  pls2_two = plsreg2(vehicles[,1:12], vehicles[,13:16], comps=NULL, crosval=TRUE)
#'  
#'  # apply plsreg2 extracting 5 components with cross-validation
#'  pls2_three = plsreg2(vehicles[,1:12], vehicles[,13:16], comps=5, crosval=TRUE)
#'  
#'  # plot variables
#'  plot(pls2_one)
#'  }
#'
plsreg2 <-
function(predictors, responses, comps = 2, crosval = TRUE)
{  
  # =======================================================
  # checking arguments
  # =======================================================
  # check predictors
  X = as.matrix(predictors)
  if (any(is.na(X))) stop("\nNo missing data are allowed")
  n = nrow(X)
  p = ncol(X)
  if (is.null(colnames(X))) colnames(X) = paste(rep("X",p), 1:p, sep="")
  if (is.null(rownames(X))) rownames(X) = 1:n
  # check responses  
  Y = as.matrix(responses)
  if (any(is.na(Y))) stop("\nNo missing data are allowed")
  if (nrow(X) != nrow(Y)) stop("\ndifferent number of rows in predictors and responses")
  q = ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) = paste(rep("Y",q), 1:q, sep="")
  if (is.null(rownames(Y))) rownames(Y) = 1:n
  if (p<2 || q<2) stop("predictors and responses must have more than one column")
  # number of components
  if (!is.null(comps))
  {
    nc = comps
    if (mode(nc) != "numeric" || length(nc) != 1 || 
      nc <= 1 || (nc%%1) != 0 || nc > min(n,p))
      nc = min(n, p)   
    if (nc == n) nc = n - 1    
  } else {
    if (n >= 10) {
      crosval = TRUE
      nc = min(n, p)
    } else {
      crosval = FALSE
      nc = 2
      message("\nSorry, no cross-validation with less than 10 observations")
    }
  }
  if (!is.logical(crosval)) crosval = FALSE
  
  # =======================================================
  # prepare ingredients
  # =======================================================
  # scaling 
  X.old = scale(X)
  Y.old = scale(Y)
  # initialize
  Wh = matrix(0, p, nc)
  Uh = matrix(0, n, nc)
  Th = matrix(0, n, nc)
  Ch = matrix(0, q, nc)
  Ph = matrix(0, p, nc)
  bh = rep(0, nc)
  # crosval
  if (crosval)
  {
    RSS = rbind(rep(n-1,q), matrix(NA, nc, q))
    PRESS = matrix(NA, nc, q)
    Q2 = matrix(NA, nc, q)
    # getting segments for CV
    sets_size = c(rep(n %/% 10,9), n-9*(n %/% 10))
    # randomize observations
    obs = sample(1:n, size=n)
    # get segments
    segments = vector("list", length=10)
    ini = cumsum(sets_size) - sets_size + 1
    fin = cumsum(sets_size)
    for (k in 1:10)
      segments[[k]] = obs[ini[k]:fin[k]]
  }
  
  # =======================================================
  # PLSR1 algorithm
  # =======================================================
  h = 1
  repeat
  {
    u.new = Y.old[,1] # first column of Y.old
    w.old = rep(1, p)
    iter = 1
    repeat
    {
      w.new = t(X.old) %*% u.new / sum(u.new^2)
      w.new = w.new / sqrt(sum(w.new^2)) # normalize w.new
      t.new = X.old %*% w.new
      c.new = t(Y.old) %*% t.new / sum(t.new^2)
      u.new = Y.old %*% c.new / sum(c.new^2)
      w.dif = w.new - w.old
      w.old = w.new
      if (sum(w.dif^2)<1e-06 || iter==100) break
      iter = iter + 1
    } 
    p.new = t(X.old) %*% t.new / sum(t.new^2)
    
    # in case of cross-validation
    if (crosval)
    {
      # cross validation "leave-one-out"
      RSS[h+1,] =  colSums((Y.old - t.new%*%t(c.new))^2)
      press = matrix(0, 10, q)
      for (i in 1:10)
      {
        aux = segments[[i]]
        uh.si = Y.old[-aux,1]
        wh.siold = rep(1, p)
        itcv = 1
        repeat
        {
          wh.si = t(X.old[-aux,]) %*% uh.si / sum(uh.si^2)
          wh.si = wh.si / sqrt(sum(wh.si^2))
          th.si = X.old[-aux,] %*% wh.si
          ch.si = t(Y.old[-aux,]) %*% th.si / sum(th.si^2)
          uh.si = Y.old[-aux,] %*% ch.si / sum(ch.si^2)
          wsi.dif = wh.si - wh.siold
          wh.siold = wh.si
          if (sum(wsi.dif^2)<1e-06 || itcv==100) break
          itcv = itcv + 1
        }
        Yhat.si = (X.old[aux,] %*% wh.si) %*% t(ch.si) 
        press[i,] = colSums((Y.old[aux,] - Yhat.si)^2)
      }
      PRESS[h,] = colSums(press)
      Q2[h,] = 1 - (PRESS[h,] / RSS[h,])
    } # end crosval
    
    # deflate
    X.old = X.old - (t.new %*% t(p.new))
    Y.old = Y.old - (t.new %*% t(c.new))
    # populate
    Wh[,h] = w.new
    Uh[,h] = u.new
    Th[,h] = t.new
    Ch[,h] = c.new
    Ph[,h] = p.new
    bh[h] = t(u.new) %*% t.new
    # do we need to stop?
    if (is.null(comps) && crosval) {
      if (sum(Q2[h,]<0.0975) == q || h == nc) break
    } else {
      if (h == nc) break
    }
    h = h + 1
  } # end repeat

  # =======================================================
  # PLSR2 results
  # =======================================================
  Th = Th[,1:h]
  Ph = Ph[,1:h]
  Wh = Wh[,1:h]
  Uh = Uh[,1:h]
  Ch = Ch[,1:h]
  Ph = Ph[,1:h]
  # modified weights
  Ws = Wh %*% solve(t(Ph) %*% Wh)
  # std beta coeffs
  Bs = Ws %*% t(Ch)
  # regular coefficients
  Br = diag(1/apply(X,2,sd)) %*% Bs %*% diag(apply(Y,2,sd))
  cte = as.vector(apply(Y,2,mean) - apply(X,2,mean)%*%Br)
  # Y predicted
  Y.hat = X %*% Br + matrix(rep(cte,each=n),n,q)
  resids = Y - Y.hat  # residuals
  # correlations
  cor.xt = cor(X, Th)
  cor.yt = cor(Y, Th)
  cor.tu = cor(Th, Uh)
  cor.xu = cor(X, Uh)
  cor.yu = cor(Y, Uh)
  # explained variance
  R2x = cor(X, Th)^2  # R2 coefficients
  R2y = cor(Y, Th)^2  # R2 coefficients
  Rdx = colMeans(R2x) # redundancy
  Rdy = colMeans(R2y) # redundancy
  EV = cbind(Rdx, cumsum(Rdx), Rdy, cumsum(Rdy))
  Rd.mat = matrix(0, h, h)
  for (j in 1:h)
    Rd.mat[1:j,j] = Rdy[1:j]
  VIP = sqrt((Wh^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), h, h))
  # if crosval
  if (crosval)
  {
    PRESS = PRESS[1:h,]
    RSS = RSS[1:(h+1),]
    Q2 = Q2[1:h,]
    Q2G = 1 - (rowSums(PRESS) / rowSums(RSS[1:h,]))
    Q2cum = Q2
    Q2cum[1,] = 1 - (PRESS[1,] / RSS[1,])
    for (i in 2:h)
      Q2cum[i,] = 1 - apply(PRESS[1:i,]/RSS[1:i,], 2, prod)
    Q2Gcum = Q2G
    for (i in 1:h)
      Q2Gcum[i] = 1 - prod((rowSums(PRESS)/rowSums(RSS[-h,]))[1:i])
    Q2T = cbind(Q2, Q2G)
    Q2TC = cbind(Q2cum,Q2Gcum)
    # add names
    q2 = c(paste(rep("Q2",q),colnames(Y),sep="."),"Q2")
    q2c = c(paste(rep("Q2cum",q),colnames(Y),sep="."),"Q2cum")
    dimnames(Q2T) = list(paste(rep("t",h), 1:h, sep=""), q2)
    dimnames(Q2TC) = list(paste(rep("t",h), 1:h, sep=""), q2c)
  } else {
    Q2T = NULL
    Q2TC = NULL
  }
  # add names
  dimnames(Wh) = list(colnames(X), paste(rep("w",h),1:h,sep=""))
  dimnames(Ws) = list(colnames(X), paste(rep("w*",h),1:h,sep=""))
  dimnames(Uh) = list(rownames(Y), paste(rep("u",h),1:h,sep=""))
  dimnames(Th) = list(rownames(X), paste(rep("t",h),1:h,sep=""))
  dimnames(Ch) = list(colnames(Y), paste(rep("c",h),1:h,sep=""))
  dimnames(Ph) = list(colnames(X), paste(rep("p",h),1:h,sep=""))
  dimnames(Bs) = list(colnames(X), colnames(Y))
  dimnames(Br) = list(colnames(X), colnames(Y))
  dimnames(cor.xt) = list(colnames(X), paste(rep("t",h),1:h,sep=""))
  dimnames(cor.yt) = list(colnames(Y), paste(rep("t",h),1:h,sep=""))
  dimnames(cor.xu) = list(colnames(X), paste(rep("u",h),1:h,sep=""))
  dimnames(cor.yu) = list(colnames(Y), paste(rep("u",h),1:h,sep=""))
  dimnames(cor.tu) = list(paste(rep("t",h),1:h,sep=""), paste(rep("u",h),1:h,sep=""))
  dimnames(EV) = list(paste(rep("t",h),1:h,sep=""), c("R2X","R2Xcum","R2Y","R2Ycum"))
  dimnames(Y.hat) = list(rownames(Y), colnames(Y))
  dimnames(resids) = list(rownames(Y), colnames(Y))
  dimnames(VIP) = list(colnames(X), paste(rep("t",h),1:h,sep=""))
  coeffs = rbind(Br, INTERCEPT=cte)
  
  # results
  structure(list(x.scores = Th, 
                 x.loads = Ph, 
                 y.scores = Uh, 
                 y.loads = Ch,
                 cor.xt = cor.xt, 
                 cor.yt = cor.yt,
                 cor.xu = cor.xu, 
                 cor.yu = cor.yu,
                 cor.tu = cor.tu,
                 raw.wgs = Wh, 
                 mod.wgs = Ws, 
                 std.coefs = Bs, 
                 reg.coefs = coeffs, 
                 y.pred = Y.hat, 
                 resid = resids,
                 expvar = EV, 
                 VIP = VIP,
                 Q2 = Q2T, 
                 Q2cum = Q2TC),
            class = "plsreg2")
}
