tegarchLogl <-
function(y, pars, lower=-Inf, upper=Inf,
  lambda.initial=NULL, logl.penalty=-1e+100, c.code=TRUE, aux=NULL)
{
if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
  chk.conds <- FALSE
}else{
  chk.conds <- TRUE
}
if(!aux$skew){ pars <- c(pars,1) }
if(!aux$asym){ pars <- c(pars[1:3],0,pars[4:5]) }

if(chk.conds){
  lambda <- tegarchRecursion(y, omega=pars[1], phi1=pars[2],
    kappa1=pars[3], kappastar=pars[4], df=pars[5],
    skew=pars[6], lambda.initial=lambda.initial, c.code=c.code,
    verbose=FALSE, aux=aux)

  term1 <- aux$iN*( log(2)-log(pars[6]+1/pars[6])+lgamma((pars[5]+1)/2)-lgamma(pars[5]/2)-log(pi*pars[5])/2 )

  term2 <- sum(lambda)

  yterm <- y + STmean(pars[5], skew=pars[6])*exp(lambda)
  num.term <- yterm^2
  denom.term <- pars[6]^(2*sign(yterm))*pars[5]*exp(2*lambda)
  term3 <- sum( (pars[5]+1)*log(1 + num.term/denom.term)/2 )

  logl <- term1 - term2 - term3

  if(is.nan(logl) || is.na(logl) || abs(logl) == Inf) logl <- logl.penalty
}else{ logl <- logl.penalty }

return(logl)
}
