print.cv.lqa <-
function (x, ...)
{
  cat ("\n")
  print (x$call)

  cat ("\nloss function: ", x$loss.func, "\n")
  cat ("validation data set used: ", x$exist.vali, "\n")
  cat ("number of folds = ", x$n.fold, "\n\n")
  
  cat ("loss matrix: \n")
  loss.mat <- x$loss.mat
  mean.array <- x$mean.array

  if (length (dim (loss.mat)) == 2)
  {
    loss.mat.extended <- cbind (loss.mat, mean.array)
    colnames (loss.mat.extended)[x$n.fold + 1] <- "mean"
    print (loss.mat.extended)
  }
  else
  {
    print (loss.mat)
    cat ("\nmean array: \n")
    print (mean.array)

  }
  cat ("\nlambda.opt = ", x$lambda.opt, "\n\n")

  print (x$best.obj)
}

