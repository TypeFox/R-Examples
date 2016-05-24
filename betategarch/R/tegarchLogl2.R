tegarchLogl2 <-
function(y, pars, lower=-Inf, upper=Inf,
  lambda.initial=NULL, logl.penalty=-10e+100, c.code=TRUE,
  aux=NULL)
{
if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
  chk.conds <- FALSE
}else{
  chk.conds <- TRUE
}
if(!aux$skew){ pars <- c(pars,1) }

if(chk.conds){
  lambda <- tegarchRecursion2(y, omega=pars[1], phi1=pars[2], phi2=pars[3],
    kappa1=pars[4], kappa2=pars[5], kappastar=pars[6], df=pars[7], skew=pars[8],
    lambda.initial=lambda.initial, c.code=c.code, aux=aux)

  #iN <- length(y)
  term1 <- aux$iN*( log(2)-log(pars[8]+1/pars[8])+lgamma((pars[7]+1)/2)-lgamma(pars[7]/2)-log(pi*pars[7])/2 )

  term2 <- sum(lambda)

  yterm <- y + STmean(pars[7], skew=pars[8])*exp(lambda)
  num.term <- yterm^2
  denom.term <- pars[8]^(2*sign(yterm))*pars[7]*exp(2*lambda)
  term3 <- sum( (pars[7]+1)*log(1 + num.term/denom.term)/2 )

  logl <- term1 - term2 - term3
  if(is.nan(logl) || is.na(logl) || abs(logl) == Inf) logl <- logl.penalty
}else{ logl <- logl.penalty }

return(logl)
}
