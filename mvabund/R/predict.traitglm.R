predict.traitglm = function(object, newR=NULL, newQ=NULL, newL=NULL, type="response", ...)
{
  # predict function. takes as arguments:
  # object - fitted "trait" object
  # newR   - for new env predictions
  # newQ   - for new traits
  # newL   -  for new spp data (only relevant for type="ll")
  # type   - "response" (default) for predicted response, "link" for linear predictor, or
  #          "ll" for mean log-likelihood
  # which.lambda - "best" (default) for the best lambda as chosen in the trait.mod call,
  #                "all" for predictions at all values of lambda considered in trait.mod.
  
  # get new polynomial values for R and Q, if required
  if(is.null(newR))
    R.des.test = object$R.des
  else
  {
    if(is.null(object$formula))
      R.des.test = get.polys(newR, object$R.des)
    else
      R.des.test = list(X=newR)
  }

  if(is.null(newQ))
    Q.des.test = object$Q.des
  else
  {
    if(is.null(object$formula))
      Q.des.test = get.polys(newQ, object$Q.des)
    else
      Q.des.test = list(X=newQ)
  }

  if(is.null(newL))
    newL = object$L
  
  n.sites = dim(R.des.test$X)[1]
  n.spp   = dim(Q.des.test$X)[1]

  if( "composition" %in% names(object$call) == FALSE )
    object$call$composition = FALSE
  if( "col.intercepts" %in% names(object$call) == FALSE )
    object$call$col.intercepts = TRUE
    
  # get new design matrix values for L
  X.des.test = get.design( R.des=R.des.test, Q.des=Q.des.test, L.names=rownames(Q.des.test$X), formula = object$formula, marg.penalty=TRUE, composition =  object$call$composition, col.intercepts = object$call$col.intercepts, any.penalty=object$any.penalty, scaling=object$scaling )
    
  #    recover()
  # get predicted eta and store in out
  out = X.des.test$X %*% coef(object)
  
  # get predicted mu (if required) and overwrite prev value of out.
  if(type=="response" | type=="logL")
  {
    # First get the family function sorted out
    if(is.character(object$family))
    {
      if (object$family == "negbinomial" || object$family == "negative.binomial")
        object$family = negative.binomial(theta=1/object$phi)
      else if (object$family == "binomial(link=logit)")
        object$family = binomial()
      else if (object$family == "binomial(link=cloglog)")
        object$family = binomial("cloglog")
      else if (object$family == "poisson")
        object$family = poisson()
      else if (object$family == "gaussian")
        object$family = gaussian()
      else
        family = get(family, mode = "function", envir = parent.frame())
    }
    
    # now use this family argument to get the inverse link function
    out = as.matrix(out)  
    out = object$family$linkinv(out)
  }
  # get mean predicted ll (if required) and overwrite prev value of out.
  if(type=="logL")
  {
    Lm   = as.matrix(newL)
    newl = as.vector(Lm)
    newl = as.matrix(newl)
    n.vars = dim(Lm)[2]
    rm(Lm)
    out = -0.5 * object$family$aic( newl, 1, out, 1, 1 )
  }
  out = matrix(as.vector(out), n.sites, n.spp)
  dimnames(out)[[1]]=dimnames(R.des.test$X)[[1]]
  dimnames(out)[[2]]=dimnames(Q.des.test$X)[[1]]
  return(out)
}
