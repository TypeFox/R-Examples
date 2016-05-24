riv_withDummies <- function(Y, Xend, Xex = NULL,
                            Zinst, dummies,
                            method = c('S-est', 'SD-est')) {
  
  # Defining parameters and variables
  n <- length(Y)
  Zinst <- as.matrix(Zinst)
  if (is.null(colnames(Zinst)))
    colnames(Zinst) <- paste('Zinst',
                             seq(ncol(Zinst)),
                             sep='')
  
  Xend <- as.matrix(Xend)
  if (is.null(colnames(Xend)))
    colnames(Xend) <- paste('Xend',
                            seq(ncol(Xend)),
                            sep='')
  
  kend <- ncol(Xend)
  k <- ncol(Zinst)
  
  if (kend > k)
    stop('The number of instruments must be equal or bigger to the number of endogenous variables')
  
  # Dummies (assumed exogenous)
  dummies <- as.matrix(dummies)
  if (is.null(colnames(dummies)))
    colnames(dummies) <- paste('D',
                               seq(ncol(dummies)),
                               sep='')
  
  # Create a full data matrix Z = cbind(Xend, Zinst, Xex, D, Y)
  
  if (is.null(Xex)) {
    X <- Xend
  } else {
    Xex <- as.matrix(Xex)
    if (is.null(colnames(Xex)))
      colnames(Xex) <- paste('Xex',
                             seq(ncol(Xex)),
                             sep='')
    X <- cbind(Xend, Xex)
  }
  
  # Size parameters
  p <- ncol(X) + ncol(dummies)          # number of covariates
  r <- k + p + 1                        # number of variables
  
  # Missing values
  if (any(is.na(c(Y, Xend, Xex, dummies, Zinst))))
    stop('missing values are not allowed')
    
  # Adapting the algorithm/estimator proposed by Yohai and Maronna, 2000
  
  # Step 1: Transform the model to remove the dummies
  
  # The intercept is considered a dummy variable and it is included in X1
  t.tilda <- coef(rq(Y ~ dummies, tau = 0.5))
  Tx.tilda <- apply(X, 2, function(x) {
    coef(rq(x ~ dummies))
  })
  Tw.tilda <- apply(Zinst, 2, function(x) {
    coef(rq(x ~ dummies))
  })
  
  # Transforming the continuous variables
  Y.tilda <- Y - cbind(1, dummies) %*% t.tilda
  X2.tilda <- X - cbind(1, dummies) %*% Tx.tilda
  Zinst.tilda <- Zinst - cbind(1, dummies) %*% Tw.tilda
  
  # RIV using the S-est
  # NOTE: The intercept is not estimated in this step! the model is centered anyway
  if (is.null(Xex)) {
    sum.RIV_S <- riv_noDummies(Y.tilda, X2.tilda,
                               Xex = NULL, Zinst.tilda,
                               intercept = FALSE,
                               method = method)
  } else {
    X2.end <- X2.tilda[, 1:kend]
    X2.ex <- X2.tilda[, -(1:kend)]
    
    sum.RIV_S <- riv_noDummies(Y.tilda, X2.end,
                               X2.ex, Zinst.tilda,
                               intercept = FALSE,
                               method = method)
  }
  beta.2 <- sum.RIV_S$Summ[, 'Coef']
  
  beta.1 <- coef(rq(Y.tilda - X2.tilda %*% beta.2 ~ dummies, tau = 0.5))
  
  
  # re-iterate to see if we can improve beta.1 without updating Zinst.tilda
  # The problem here is that we are not using the dummy variables as instruments.
  # The 2SLS will use all of them
  
  for (i in 1:10) {
    Y.tilda <- Y - cbind(1, dummies) %*% beta.1
    if (is.null(Xex)) {
      sum.RIV_S <- riv_noDummies(Y.tilda, Xend,
                                 Xex = NULL, Zinst,
                                 intercept = FALSE,
                                 method = method)
      
    } else {
      sum.RIV_S <- riv_noDummies(Y.tilda, Xend,
                                 Xex, Zinst,
                                 intercept = FALSE,
                                 method = method)
    }
    
    beta.2 <- sum.RIV_S$Summ[, 'Coef']
    beta.1 <- coef(rq(Y - X %*% beta.2 ~ dummies, tau = 0.5))
  }
  
  
  tabRIV <- matrix(c(beta.1, beta.2), ncol = 1,
                   dimnames=list(c('Intercept',
                                   colnames(dummies),
                                   colnames(X)),
                                 'Coef'))
  
  list(Summary.Table = tabRIV)
}
