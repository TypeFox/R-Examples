#'@title PLS-R1: Partial Least Squares Regression 1
#'
#'@description
#'The function plsreg1 performs Partial Least Squares Regression for the univariate case (i.e. one
#'response variable)
#'
#'@details
#'The minimum number of PLS components (\code{comps}) to be extracted is 2.
#'
#'The data is scaled to standardized values (mean=0, variance=1).
#'
#'The argument \code{crosval} gives the option to perform cross-validation.
#'This parameter takes into account how \code{comps} is specified.
#'When \code{comps=NULL}, the number of components is obtained by cross-validation.
#'When a number of components is specified, cross-validation results are calculated
#'for each component.
#'
#'@param predictors A numeric matrix or data frame with the predictor variables (which
#'may contain missing data).
#'@param response A numeric vector for the reponse variable. 
#'No missing data allowed.
#'@param comps The number of extracted PLS components (2 by default).
#'@param crosval Logical indicating whether cross-validation should be performed
#'(\code{TRUE} by default). No cross-validation is done if there is missing data or
#'if there are less than 10 observations.
#'@return An object of class \code{"plsreg1"}, basically a list with the
#'following elements:
#'@return \item{x.scores}{PLS components (also known as T-components)}
#'@return \item{x.loads}{loadings of the predictor variables}
#'@return \item{y.scores}{scores of the response variable (also known as U-components)}
#'@return \item{y.loads}{loadings of the response variable}
#'@return \item{cor.xyt}{Correlations between the variables and the PLS
#'components}
#'@return \item{raw.wgs}{weights to calculate the PLS scores with the deflated
#'matrices of predictor variables}
#'@return \item{mod.wgs}{modified weights to calculate the PLS scores with the
#'matrix of predictor variables}
#'@return \item{std.coefs}{Vector of standardized regression coefficients}
#'@return \item{reg.coefs}{Vector of regression coefficients (used with the original
#'data scale)}
#'@return \item{R2}{Vector of PLS R-squared}
#'@return \item{R2Xy}{explained variance of variables by PLS-components}
#'@return \item{y.pred}{Vector of predicted values}
#'@return \item{resid}{Vector of residuals}
#'@return \item{T2}{Table of Hotelling T2 values (used to detect atypical
#'observations)}
#'@return \item{Q2}{Table with the cross validation results. Includes: PRESS, RSS,
#'Q2, and cummulated Q2. Only available when \code{crosval=TRUE}}
#'@author Gaston Sanchez
#'@seealso \code{\link{plot.plsreg1}}, \code{\link{plsreg2}}.
#'@references Geladi, P., and Kowalski, B. (1986) Partial Least Squares
#'Regression: A Tutorial. \emph{Analytica Chimica Acta}, \bold{185}, pp. 1-17.
#'
#'Tenenhaus, M. (1998) \emph{La Regression PLS. Theorie et Pratique}. Editions
#'TECHNIP, Paris.
#'
#'Tenenhaus, M., Gauchi, J.-P., and Menardo, C. (1995) Regression PLS et
#'applications. \emph{Revue de statistique appliquee}, \bold{43}, pp. 7-63.
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## example of PLSR1 with the vehicles dataset
#'  # predictand variable: price of vehicles
#'  data(vehicles)
#'  
#'  # apply plsreg1 extracting 2 components (no cross-validation)
#'  pls1_one = plsreg1(vehicles[,1:12], vehicles[,13,drop=FALSE], comps=2, crosval=FALSE)
#'  
#'  # apply plsreg1 with selection of components by cross-validation
#'  pls1_two = plsreg1(vehicles[,1:12], vehicles[,13,drop=FALSE], comps=NULL, crosval=TRUE)
#'  
#'  # apply plsreg1 extracting 5 components with cross-validation
#'  pls1_three = plsreg1(vehicles[,1:12], vehicles[,13,drop=FALSE], comps=5, crosval=TRUE)
#'  
#'  # plot variables
#'  plot(pls1_one)
#'  }
#'
plsreg1 <-
  function(predictors, response, comps = 2, crosval = TRUE)
  {
    # =======================================================
    # checking arguments
    # =======================================================
    # check predictors
    X = as.matrix(predictors)
    n = nrow(X)
    p = ncol(X)
    if (p < 2) stop("\npredictors must contain more than one column")
    if (is.null(colnames(X)))
      colnames(X) = paste(rep("X",p), 1:p, sep="")
    if (is.null(rownames(X)))
      rownames(X) = 1:n
    # check response
    Y = as.matrix(response)
    if (ncol(Y) != 1)
      stop("\nresponse must be a single variable")
    if (any(is.na(response))) 
      stop("\nresponse must not contain missing values")    
    if (nrow(X) != nrow(Y)) 
      stop("\npredictors and response have different number of rows")
    if (is.null(colnames(Y)))
      colnames(Y) = "Y"
    if (is.null(rownames(Y)))
      rownames(Y) = 1:n
    # number of components
    if (any(is.na(X))) na.miss = TRUE else na.miss = FALSE
    if (!is.null(comps))
    {
      nc = comps
      if (mode(nc) != "numeric" || length(nc) != 1 || 
        nc <= 1 || (nc%%1) != 0 || nc > min(n,p))
        nc = min(n, p)   
      if (nc == n) nc = n - 1    
    } else {
      if (na.miss) { # comps=NULL && NA's 
        crosval = FALSE 
        nc = 2
      } else { # comps=NULL && complete data
        if (n >= 10) crosval = TRUE else crosval = FALSE
        nc = min(n, p)
      }
    }
    if (!is.logical(crosval)) crosval = FALSE
    
    # =======================================================
    # prepare ingredients
    # =======================================================
    # scaling 
    Xx = scale(X) 
    Yy = scale(Y)
    # initialize
    X.old = Xx
    Y.old = Yy
    Th = matrix(NA, n, nc)  # matrix of X-scores
    Ph = matrix(NA, p, nc)  # matrix of X-loadings
    Wh = matrix(NA, p, nc)  # matrix of raw-weights
    Uh = matrix(NA, n, nc)  # matrix of Y-scores
    ch = rep(NA, nc)  # vector of y-loadings
    Hot = matrix(NA, n, nc)  # matrix of T2 Hotelling  
    hlim = rep(NA, nc)  # vector of T2 thresholds
    # crosval
    if (crosval)
    {
      RSS = c(n - 1, rep(NA, nc))
      PRESS = rep(NA, nc)
      Q2 = rep(NA, nc)
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
    w.old = rep(1, p)
    t.new = rep(1, n)
    p.new = rep(NA, p)
    h = 1
    repeat
    {
      if (na.miss) 
      {
        ## missing data
        for (j in 1:p) {
          i.exist = which(complete.cases(X[,j]))
          w.old[j] = sum(X.old[i.exist,j] * Y.old[i.exist])
        }
        w.new = w.old / sqrt(sum(w.old^2))       
        for (i in 1:n) {
          j.exist = which(complete.cases(X[i,]))
          t.new[i] = sum(X.old[i,j.exist] * w.new[j.exist])
        }
        for (j in 1:p) {
          i.exist = intersect(which(complete.cases(X[,j])), which(complete.cases(t.new)))
          p.new[j] = sum(X.old[i.exist,j] * t.new[i.exist]) / sum(t.new[i.exist]^2)
        }
        c.new = t(Y.old) %*% t.new / sum(t.new^2)
        u.new = Y.old / as.vector(c.new)
      }
      if (!na.miss) # no missing data
      {
        ## no missing data
        w.old = t(X.old) %*% Y.old / sum(Y.old^2)
        w.new = w.old / sqrt(sum(w.old^2)) # normalization
        t.new = X.old %*% w.new
        p.new = t(X.old) %*% t.new / sum(t.new^2)
        c.new = t(Y.old) %*% t.new / sum(t.new^2)
        u.new = Y.old / as.vector(c.new)
        # in case of cross-validation
        if (crosval)
        {
          # cross validation
          RSS[h+1] =  sum((Y.old - t.new%*%c.new)^2)
          press = rep(0, 10)
          for (i in 1:10)
          {
            aux = segments[[i]]
            Xy.aux = t(X.old[-aux,]) %*% Y.old[-aux]
            wh.si = Xy.aux %*% sqrt(solve(t(Xy.aux) %*% Xy.aux))
            th.si = X.old[-aux,] %*% wh.si
            ch.si = t(Y.old[-aux]) %*% th.si %*% solve(t(th.si) %*% th.si)
            ch.si = as.vector(ch.si)
            Yhat.si = ch.si * X.old[aux,] %*% wh.si
            press[i] = sum((Y.old[aux] - Yhat.si)^2)
          }
          PRESS[h] = sum(press)
          Q2[h] = 1 - PRESS[h]/RSS[h]
        } # end crosval  
      }
      # deflate
      Y.old = Y.old - (t.new %*% c.new)
      X.old = X.old - (t.new %*% t(p.new))
      # populate
      Th[,h] = t.new
      Ph[,h] = p.new
      Wh[,h] = w.new
      Uh[,h] = u.new
      ch[h] = c.new   
      Hot[,h] = (n/(n-1)) * t.new^2 / (sum(t.new^2)/(n-1))
      hlim[h] = qf(.95, h, n-h) * (h*(n^2-1))/(n*(n-h))
      # do we need to stop?
      if (is.null(comps) && crosval) {
        if (Q2[h] < 0.0975 || h == nc) break
      } else {
        if (h == nc) break
      }
      h = h + 1
    } # end repeat
    
    # =======================================================
    # PLSR1 results
    # =======================================================
    if (crosval) {
      q2cum = rep(NA, h)
      for (k in 1:h)
        q2cum[k] = prod(PRESS[1:k]) / prod(RSS[1:k])
      Q2cum = 1 - q2cum
      Q2cv = cbind(PRESS[1:h], RSS[1:h], Q2[1:h], rep(0.0975,h), Q2cum)
      dimnames(Q2cv) = list(1:h, c("PRESS","RSS","Q2","LimQ2","Q2cum"))
      if (is.null(comps))
        h = h - 1
    }
    if (!crosval) Q2cv = NULL
    Th = Th[,1:h]
    Ph = Ph[,1:h]
    Wh = Wh[,1:h]
    Uh = Uh[,1:h]
    ch = ch[1:h]
    # modified weights
    Ws = Wh %*% solve(t(Ph) %*% Wh) 
    # std beta coeffs
    Bs = as.vector(Ws %*% ch)
    if (!na.miss) 
    {
      # beta coeffs
      Br = Bs * (rep(apply(Y,2,sd),p) / apply(X, 2, sd))
      # intercept
      cte = as.vector(colMeans(Y) - Br %*% apply(X, 2, mean))
      # y predicted
      y.hat = as.vector(X%*%Br + cte)
      # correlations
      cor.xyt = cor(cbind(Xx, y=Yy), Th)
    } else {
      mu.x <- attributes(Xx)$'scaled:center'
      sd.x <- attributes(Xx)$'scaled:scale'
      X.hat = Th%*%t(Ph) %*% diag(sd.x,p,p) + matrix(rep(mu.x,each=n),n,p)
      # beta coeffs (unstandardized)
      Br = Bs * (rep(apply(Y,2,sd),p) / sd.x)
      # intercept
      cte = as.vector(colMeans(response) - Br%*%mu.x)
      y.hat = as.vector(X.hat%*%Br + cte)
      cor.xyt = matrix(NA, p+1, h)
      for (j in 1:p) {
        i.exist <- which(complete.cases(X[,j]))
        cor.xyt[j,] = cor(Xx[i.exist,j], Th[i.exist,])
      }
      cor.xyt[p+1,] = cor(Yy, Th)
    }
    resid = as.vector(Y - y.hat) # residuals
    R2 = as.vector(cor(Th, Yy))^2  # R2 coefficients
    # explained variance
    R2Xy = t(apply(cor.xyt^2, 1, cumsum))
    # Hotelling T2
    T2hot = rbind(hlim[1:h], t(apply(Hot[,1:h],1,cumsum)))
    # add names
    dimnames(Wh) = list(colnames(X), paste(rep("w",h), 1:h, sep=""))
    dimnames(Ws) = list(colnames(X), paste(rep("w*",h), 1:h, sep=""))    
    dimnames(Th) = list(rownames(X), paste(rep("t",h), 1:h, sep=""))
    dimnames(Ph) = list(colnames(X), paste(rep("p",h), 1:h, sep=""))
    dimnames(Uh) = list(rownames(Y), paste(rep("u",h), 1:h, sep=""))
    names(ch) = paste(rep("c",h), 1:h, sep="")
    dimnames(T2hot) = list(c("T2",rownames(X)), paste(rep("H",h), 1:h, sep=""))
    names(Bs) = colnames(X)
    names(Br) = colnames(X)
    names(resid) = rownames(Y)
    names(y.hat) = rownames(Y)
    names(R2) = paste(rep("t",h), 1:h, sep="")
    colnames(R2Xy) = paste(rep("t",h), 1:h, sep="")
    dimnames(cor.xyt) = list(c(colnames(X), colnames(Y)), colnames(Th))
    
    # results
    res = list(x.scores = Th, 
               x.loads = Ph,
               y.scores = Uh, 
               y.loads = ch, 
               cor.xyt = cor.xyt, 
               raw.wgs = Wh,
               mod.wgs = Ws, 
               std.coefs = Bs, 
               reg.coefs = c(Intercept=cte, Br), 
               R2 = R2,
               R2Xy = R2Xy,
               y.pred = y.hat, 
               resid = resid,
               T2 = T2hot,
               Q2 = Q2cv,
               y = response)    
    class(res) = "plsreg1"
    return(res)
  }


