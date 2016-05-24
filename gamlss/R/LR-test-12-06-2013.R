#---------------------------------------------------------------------------------------
# this is to test two nested GAMLSS models using a chi-square test 
#even thought no check is done  whether the models are nested
LR.test <- function(null, alternative, print=TRUE)
  {
  if (!("gamlss"%in%class(null))) stop("the null model is not a gamlss model") # ammended 12-06-2013
  if (!("gamlss"%in%class(alternative))) stop("the alternative model is not a gamlss model")
     D0 <- deviance(null)
     D1 <- deviance(alternative)
  if (D1>D0) stop("The null model has smaller deviance than the alternative \n",
                  "The models should be nested" )
    df0 <- null$df.fit
    df1 <- alternative$df.fit
    chi <- D0-D1
     df <- df1-df0
     if (df<0) stop("The difference in df's is negative. \n",
                    "Are the models nested? \n" )
  p.val <- 1-pchisq(chi,df)
  if (print)
  {
  cat(" Likelihood Ratio Test for nested GAMLSS models. \n",
  "(No check whether the models are nested is performed). \n \n",
      "      Null model: deviance=", D0, "with ", df0, "deg. of freedom \n",
      "Altenative model: deviance=", D1, "with ", df1, "deg. of freedom \n \n",
      "LRT =", chi, "with",  df,   "deg. of freedom and p-value=", p.val, "\n")
  }
  else
  {
   return(list(chi=chi, df=df, p.val=p.val))
  }
  }
#----------------------------------------------------------------------------------------
# this does not seems appropriate sinse the demominator can be negative
# Here I am taking abs() to avoid this
#F.Test <- function(null, alternative)
#  {
#  if (!(class(null)[1]%in%"gamlss")) stop("the null model is not a gamlss model")
#  if (!(class(alternative)[1]%in%"gamlss")) stop("the alternative model is not a gamlss model")
#  D0 <- deviance(null)
#  D1 <- deviance(alternative)
#  if (D1>D0) stop("the null model has smaller deviance than the alternative")
#  df0 <- null$df.fit
#  df1 <- alternative$df.fit
#  nom <- D0-D1
#  dfnom  <- df1-df0
#  dfden  <-   alternative$df.fit
#  F <- (nom/dfnom)/(abs(D1)/dfden)
#  p.val <- 1-pf(nom,dfnom, dfden)
#  list(F=F, df1=dfnom, df2=dfden, p.val=p.val)
#  }
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
