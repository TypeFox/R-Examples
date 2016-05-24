predict.wsrf <- function(object,
                         newdata,
                         type=c("response", "class", "prob", "vote", "aprob", "waprob"), ...)
{
  if (!inherits(object, "wsrf")) 
    stop("Not a legitimate wsrf object")

  if (missing(type)) type <- "class"
  if (type=="response") type <- "class"

  type <- match.arg(type)

  # The C code for predict does not handle missing values. So handle
  # them here by removing them from the dataset and then add in, in
  # the correct places, NA as the results from predict.

  complete <- complete.cases(newdata)
  rnames   <- rownames(newdata)
  newdata  <- newdata[complete,]

  hasmissing <- !all(complete)
  nobs       <- length(complete)
  
  # function "predict()" in C returns "votes" by default, 
  # and can also directly returns "aprob" or "waprob" correspondingly in terms of <type>
  # but "class" and "prob" will be treated as "votes",
  # so "class" and "prob" still need to be calculated in R below

  # .Call("predict") returns a factor vector of labels or 
  # a numeric matrix of nclass * nobservations with labels as rownames.

  res <- .Call("predict", object, newdata, type, PACKAGE="wsrf")
  if (type != "class") res <- t(res)

  # Deal with observations with missing values.

  if (hasmissing)
  {
    if (type == "class")
    {
      cl <- factor(rep(NA, nobs), levels=levels(res))
      cl[complete] <- res
      names(cl) <- rnames
      return(cl)
    }
    else
    {
      fin <- matrix(NA_real_, nrow=nobs, ncol=ncol(res))
      fin[complete, ] <- res
      rownames(fin) <- rnames
      colnames(fin) <- colnames(res)
      return(fin)
    }
  }

  if (type == "class")
    names(res) <- rnames
  else
    rownames(res) <- rnames

  return(res)
}
