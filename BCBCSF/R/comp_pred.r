
probs_attrue <- function (probs_pred, y)
{
  tp <- rep(0, nrow(probs_pred))
  names (tp) <- rownames (probs_pred)
  for(i in 1:nrow(probs_pred)) tp[i] <- probs_pred[i,y[i]]

  tp
}

comp_amlp <- function(probs_pred, y)
{
  mean (-log (probs_attrue (probs_pred, y)))
}

## Mloss -- a matrix specifying losses, with row for true values, and
## column for predicted values.
comp_loss <- function(probs_pred, y, Mloss = NULL)
{
     G <- ncol (probs_pred)

     if (is.null (Mloss))
     {
        Mloss <- matrix(1,G,G)
        diag(Mloss) <- 0
     }

     loss_pred <- probs_pred %*% Mloss
     y_pred <- apply(loss_pred,1,which.min)

     loss <- 0
     for(i in 1:nrow(probs_pred)) {
         loss <- loss + Mloss[y[i],y_pred[i]]
     }

     loss / length (y)
}

comp_eer <- function (probs_pred)
{
 mean (1 - apply (probs_pred, 1, max))
}

eval_pred <- function (out_pred, y_ts, Mloss = NULL)
{
  array_probs_pred <- out_pred$array_probs_pred
  nos_fsel <- out_pred$nos_fsel

  amlp <- er <- loss <- NULL

  amlp <- apply (array_probs_pred, 3, comp_amlp, y = y_ts)
  er <- apply (array_probs_pred, 3, comp_loss, y = y_ts)
  probs_at_truelabels <- apply (array_probs_pred, 3, probs_attrue, y = y_ts)

  summary <- 
  data.frame (No.Features = nos_fsel, Error.Rate = er, AMLP = amlp)
   
  if (!is.null (Mloss))
  {
      loss <- apply (array_probs_pred, 3, comp_loss, y = y_ts, Mloss = Mloss)
      result <- cbind (summary, Loss = loss)
  }
  
  list (probs_at_truelabels = probs_at_truelabels, summary = summary)
}


## partition all cases into nfold subsets
## This function partitions a set of observations into subsets of almost
## equal size. The result is used in crossvalidation
mk_folds <- function(y, nfold = 10, random = FALSE)
{
    n <- length (y)
    nos_g <- table (y)
    G <- length (nos_g)
    nfold <- min (nfold, n)
    
    reduced.rep <- TRUE
    while (reduced.rep){
       if (!random)        folds <- rep (1:nfold, length = n)
       else folds <- sample(rep (1:nfold, length=n))
       ## check any fold has reduced class representation
       reduced.rep <- FALSE
       for (i in 1:nfold)
       {
           G_ifold <- length(unique (y[folds!=i]))
           if (G_ifold < G){
               reduced.rep <- TRUE
               random <- TRUE
               break
           } 
       }
    }
    
    
    ## create fold list 
    foldlist <- rep (list (""),nfold)
    for (i in 1:nfold)
    {
        foldlist [[i]] <- which (folds == i)
    }
    
    foldlist
 }


#################### a generic crossvalidation function ####################
## X --- features with rows for cases
## y --- a vector of response values
## nfold --- number of folds in cross validation
##  fitpred_func --- function for training and prediction:
## 	the arguments of fitpred_func must include X_tr, y_tr, X_ts
## 	the outputs of fitpred_func must include probs_pred
## ... --- other arguments needed by fitpred_func other than X_tr, y_tr, X_ts
cross_vld <- function (
     X, y, nfold = 10, folds = NULL, fitpred_func = bcbcsf_fitpred,  ...)
{
  if (!is.matrix(X)) stop ("'X' must be a matrix with rows for cases")

  n <- nrow(X)


  if (is.null (folds))
  {
    folds <- mk_folds (y, nfold, random = FALSE)
  }

  nfold <- length (folds)

  array_probs_pred <- NULL
  vector_ts <- NULL

  for (i_test in 1:nfold)
  {
    cat(sprintf ("Fold%2d: ", i_test) )
    ts <- folds [[i_test]]
    vector_ts <- c (vector_ts, ts)
    tr <- (1:n)[- (ts)]

    onetrpr <- fitpred_func (
          X_tr = X[tr,, drop = FALSE], y_tr = y[tr],
          X_ts = X[ts,, drop = FALSE], ...)
    one_array_probs_pred <- onetrpr$array_probs_pred
    array_probs_pred <- abind ( array_probs_pred, one_array_probs_pred,
      along = 1)
    onetrpr <- onetrpr[names(onetrpr) != "array_probs_pred"]
  }
  cat ("\n")
  
  ## make the order of cases in array_probs_pred is the same as X
  array_probs_pred <- array_probs_pred[order (vector_ts),,, drop = FALSE]
  dims <- dim (array_probs_pred)
  dimnames (array_probs_pred)  <- list(paste("Case", 1:dims[1], sep=""),
                                       paste("Class", 1:dims[2], sep=""),
                                       paste("fsel", 1:dims[3], sep=""))
  
  
  #dimnames (array_probs_pred) [[1]] <- paste("Case", 1:n, sep="")

  c (onetrpr, list (folds = folds, array_probs_pred = array_probs_pred) )
}



