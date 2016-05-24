`lambda.stem.ci` <-
function(tb, r, eps=0)
{
  ci<- list()
  betaF <- (exp(r*tb)-1)/(exp(r*tb)-eps)
  ci$upper <- 1 + log(.025)/log(betaF)
  ci$lower <-  1 + log(.975)/log(betaF)
  ci <-as.data.frame(ci)
  ci
}

