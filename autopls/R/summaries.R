# Summary and print functions for autopls objects

summary.autopls <- function (object, ...)
{

  # Get parameters
  iter <- object$metapls$current.iter
  lv <- get.lv (object)
  N <- length (object$metapls$Y)
  pred <- sum (object$predictors) 
  val <- object$metapls$val
  prep <- object$metapls$preprocessing
  test <- ifelse (is.null (object$metapls$X.testset), FALSE, TRUE)
  
  if (!test) 
  {
    r2.all <- unlist (R2 (object, c('train', 'CV'), 
      nc = lv, ic = FALSE))
    rmse.all <- unlist (RMSEP (object, c('train', 'CV'), 
      nc = lv, ic = FALSE))
  }
  else 
  {
    sXt <- object$metapls$X.testset [,object$predictors]
    Yt <- object$metapls$Y.testset
    if (prep != 'none') sXt <- prepro (sXt, method = 'bn')
    sYX <- data.frame (Y = Yt, X = I (sXt))
    r2.all <- unlist (R2 (object, c('train', 'CV', 'test'), 
      nc = lv, ic = FALSE, newdata = sYX))  
    rmse.all <- unlist (RMSEP (object, c('train', 'CV', 'test'), 
      nc = lv, ic = FALSE, newdata = sYX))  
    rmse.test <- rmse.all$val3
    r2.test <- r2.all$val3
  }
  
  r2.cal <- r2.all$val1
  r2.val <- r2.all$val2
  if (test) r2.test <- r2.all$val3
  rmse.cal <- rmse.all$val1
  rmse.val <- rmse.all$val2
  if (test) rmse.test <- rmse.all$val3
  scaling <- object$metapls$scaling
  preprocessing <- object$metapls$preprocessing

  # Output object
  output <- list (predictors          = pred,
                  observations        = N,
                  lv                  = lv,
                  run                 = iter,
                  val                 = val,
                  rmse.cal            = rmse.cal,
                  rmse.val            = rmse.val,
                  rmse.test           = NA,
                  r2.cal              = r2.cal,
                  r2.val              = r2.val,
                  r2.test             = NA,
                  scaling             = scaling,
                  preprocessing       = preprocessing)

                  if (test) 
                  {
                    output$rmse.test  <- rmse.test
                    output$r2.test    <- r2.test
                  }                 
  return (output)
}

print.autopls <- function (x, ...)
{

  obj <- summary (x)
  test <- ifelse (is.na (obj$rmse.test), FALSE, TRUE)

  # Screen output
  cat ('\n')
  cat (paste ('Predictors:', 
    obj$predictors, '  '))
  cat (paste ('Observations:', 
    obj$observations, '  '))
  cat (paste ('Latent vectors:', 
    obj$lv, '  '))
  cat (paste ('Run:', 
    obj$run, '\n'))
  cat (paste ('RMSE(CAL):', 
    format (obj$rmse.cal, digits = 3)), '  ')
  cat (paste ('RMSE(', obj$val, '): ', 
    format (obj$rmse.val, digits = 3), '   ', sep = ''))
  if (test) cat (paste ('RMSE(test):', 
    format (obj$rmse.test, digits = 3)))
  cat ('\n')
  cat (paste ('R2(CAL):', 
    format (obj$r2.cal, digits = 3), '   '))
  cat (paste ('R2(', obj$val, '): ', 
    format (obj$r2.val, digits = 3), '   ', sep = ''))
  if (test) cat (paste ('R2(test):', 
    format (obj$r2.test, digits = 3)))
  cat ('\n')

  # Output
  invisible (obj)
}

summary.slim <- function (object, ...)
{

  # Get parameters
  iter <- object$metapls$current.iter
  lv <- object$metapls$current.lv
  N <- object$slimobj$N
  pred <- sum (object$predictors)
  val <- object$metapls$val
  r2.cal <- object$slimobj$r2$val [1]
  r2.val <- object$slimobj$r2$val [2]
  rmse.cal <- object$slimobj$rmse$val [1]
  rmse.val <- object$slimobj$rmse$val [2]

  # Output
  output <- list (N          = N,
                  predictors = pred,
                  iter       = iter,
                  lv         = lv,
                  val        = val,
                  rmse.cal   = rmse.cal,
                  rmse.val   = rmse.val,
                  r2.cal     = r2.cal,
                  r2.val     = r2.val)
  
  return (output)
}

print.slim <- function (x, ...)
{

  obj <- summary (x)
  
  # Screen output
  cat ('\n')
  cat (paste ('Predictors:', obj$predictors, '  '))
  cat (paste ('Observations:', obj$N, '  '))
  cat (paste ('Latent vectors:', obj$lv, '  '))
  cat (paste ('Run:', obj$iter, '\n'))
  cat (paste ('RMSE(CAL):', 
    format (obj$rmse.cal, digits = 3), '   '))
  cat (paste ('RMSE(', obj$val, '): ', 
    format (obj$rmse.val, digits = 3), '\n', sep = ''))
  cat (paste ('R2(CAL):', 
    format (obj$r2.cal, digits = 3), '   '))
  cat (paste ('R2(', obj$val, '): ',
    format (obj$r2.val, digits = 3), '\n\n', sep = ''))

  # Output
  invisible (obj)
}
