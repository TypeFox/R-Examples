cv.nng <- function (y.train, x.train, intercept = TRUE, y.vali = NULL, x.vali = NULL, lambda.nng, family, penalty, standardize = TRUE, n.fold, cv.folds, loss.func = aic.loss, control = lqa.control (), ...)
{

  call <- match.call ()

  if ((var (x.train[,1]) > control$var.eps) & (intercept == TRUE))    # adjusting x.train if column of ones is missing (in the 'intercept = TRUE' case)
    x.train <- cbind (1, x.train)

  if (!is.null (x.vali))
  {
    if ((var (x.vali[,1]) > control$var.eps) & (intercept == TRUE))   # adjusting x.vali if column of ones is missing (in the 'intercept = TRUE' case)
      x.vali <- cbind (1, x.vali)
  }


  if (missing (loss.func))
  {
    loss.func <- "aic.loss"
    warning ("aic.loss is used for loss function")
  }

  if (missing (family))
    stop ("'family' must be specified")


  if (!is.character (loss.func))
    stop ("loss.func must be given as character \n")

  lambda.candidates <- lambda.nng

  x.train <- as.matrix (x.train)
  nobs.train <- length (y.train)
  digits <- 5 #control$digits

  if (loss.func == "gcv.loss")
  {
    x.vali <- x.train
    y.vali <- y.train
  }

  if (is.null (x.vali) || is.null (y.vali))
  {
    exist.vali <- FALSE

    if (missing (n.fold))
      stop ("Either the number of cv splits or a validation data set must be given! \n")

    if (missing (cv.folds))
      cv.folds <- split (sample (seq (nobs.train)), rep (1 : n.fold, length = nobs.train))
  }
  else
  {
    exist.vali <- TRUE
    x.vali <- as.matrix (x.vali)
    n.fold <- 1
    cv.folds <- list (1 : nobs.train)
  }

##### Construction the loss.array (critical issue is getting the number of necessary dimensions!): ##########

  if (is.list (lambda.candidates))
    {
      ll1 <- length (lambda.candidates)    # identifies the dimension of tuning parameter 'lambda' 
      seq1 <- 1 : ll1
      no <- sapply (seq1, function (seq1) {length (lambda.candidates[[seq1]])})  # identifies the number of TP candidates per dimension

      if (ll1 < 3)
      {
        no <- c (no, rep (1, 3-ll1))
        lambda.candidates[[3]] <- 1
      }
      if (ll1 < 2)
        lambda.candidates[[2]] <- 1    # just alibi names in order to set the dimnames...
   
    } 
    else
    {
      ll1 <- 1
      seq1 <- 1
      no <- c (length (lambda.candidates), 1, 1)
    }

  loss.array <- tr.array <- array (0, dim = c (no[1], n.fold, no[-1]))
  dimnames (loss.array) <- list (paste ("lambda1 = ", round (lambda.candidates[[1]], digits = digits), sep = ""), paste ("fold ", 1 : n.fold, sep = ""), lambda2 = round (lambda.candidates[[2]], digits = digits), lambda3 = round (lambda.candidates[[3]], digits = digits))
  current.lambda <- rep (0, ll1)     # for initialization...


### Fill in the 'loss.array': ##############################################################################

  for (i in 1 : n.fold)
  {
    if (exist.vali)
    {
      cy.train <- y.train
      cx.train <- x.train
      cy.vali <- y.vali
      cx.vali <- x.vali
    }
    else
    {
      cy.train <- y.train[-cv.folds[[i]]] 
      cx.train <- x.train[-cv.folds[[i]],]
      cy.vali <- y.train[cv.folds[[i]]]
      cx.vali <- x.train[cv.folds[[i]],]
    }

    for (i1 in 1 : no[1])
      for (i2 in 1 : no[2])
        for (i3 in 1 : no[3])
        {
           current.no <- c (i1, i2, i3)
           current.lambda <- sapply (seq1, function (seq1) {current.lambda[seq1] <- lambda.candidates[[seq1]][current.no[seq1]]})

           train.obj <- if (ncol (cx.train) < nrow (cx.train))
                            lqa.default (x = cx.train, y = cy.train, lambda.nng = current.lambda, family = family, penalty = penalty, intercept = intercept, standardize = standardize, control = control, method = "nng.update", ...)
                         else
                            lqa.default (x = cx.train, y = cy.train, lambda.nng = current.lambda, family = family, penalty = penalty, intercept = intercept, standardize = standardize, control = control, method = "nng.update", ...)

           pred.obj <- predict.lqa (train.obj, new.x = cx.vali, new.y = cy.vali, ...)

           loss.array[i1,i,i2,i3] <- do.call (loss.func, list (pred.obj))
           tr.array[i1,i,i2,i3] <- train.obj$fit.obj$tr.H
        }
  }

### Find the optimal tuning parameters: #####################################################################

  mean.array <- array (0, dim = no)
  dimnames (mean.array) <- list (paste ("lambda1 = ", round (lambda.candidates[[1]], digits = digits), sep = ""), paste ("lambda2 = ", round (lambda.candidates[[2]], digits = digits), sep = ""), lambda3 = round (lambda.candidates[[3]], digits = digits))
    for (i1 in 1 : no[1])
      for (i2 in 1 : no[2])
        for (i3 in 1 : no[3])
           mean.array[i1,i2,i3] <- mean (loss.array[i1,,i2,i3])

  if (dim (loss.array)[2] > 1)
    loss.array <- drop (loss.array)
 
  best.pos <- which (mean.array == min (mean.array), arr.ind = TRUE)[1,]  # bei mehrfachen Minima nehmen wir einfach die erste Zeile!
  lambda.opt <- sapply (seq1, function (seq1) {current.lambda[seq1] <- lambda.candidates[[seq1]][best.pos[seq1]]})
  names (lambda.opt) <- "lambda"

  best.obj <- if (ncol (cx.train) < nrow (cx.train))
                 lqa.default (x = cx.train, y = cy.train, lambda.nng = lambda.opt, family = family, penalty = penalty, intercept = intercept, standardize = standardize, control = control, method = "nng.update", ...)
              else
                 lqa.default (x = cx.train, y = cy.train, lambda.nng = lambda.opt, family = family, penalty = penalty, intercept = intercept, standardize = standardize, control = control, method = "nng.update", ...)

  cv.obj <- list (call = call, lambda.opt = lambda.opt, beta.opt = best.obj$coef, best.pos = best.pos, loss.mat = loss.array, best.obj = best.obj, loss.func = loss.func, exist.vali = exist.vali, cv.folds = cv.folds, n.fold = n.fold, mean.array = drop (mean.array), lambda.candidates = lambda.candidates, tr.array = tr.array)
  class (cv.obj) <- c ("cv.lqa")
  cv.obj
}
