irates.ratio <-
function ( 	irates, 
							covar.code, 
							ci.fun = NULL, 
							ci.level = NULL )
{
 
  object = irates
  event.code = object$event.code
  event.lab = object$event.lab

  covar.code = as.character(covar.code)

  given.covar.code = c(object$covar.code, object$full.sample.code)
  given.covar.lab = c(object$covar.lab, object$full.sample.lab)  
  
  if(length(covar.code) != 2){
  	stop("Only two different covariate values can be compared")
  	}
 
  if(any(!(covar.code %in% given.covar.code))){
  	stop(paste("covar.code", paste(covar.code[which(!(covar.code %in% given.covar.code))], collapse = ", "), "is not contained in irates object"))
  	}
  	
  covar.lab = rep(NA, length(covar.code))
  for(l in 1:length(given.covar.code)){
  	covar.lab[which(covar.code == given.covar.code[l])] = given.covar.lab[l]
  	}
  
  if ( is.null(ci.fun) ){ci.fun = object$ci.fun}
  else if(!(ci.fun %in% c("lin", "log"))) stop(paste("ci.fun", ci.fun, "is not implemented in irates.ratio"))

  if ( is.null(ci.level) )
    ci.level = object$ci.level
  
  quant = 1-(1-ci.level)/2

  empty.vector = rep(NA,length(event.code))
  names(empty.vector) = event.code
  
  irr <- var <- conf.lower <- conf.upper <- empty.vector

  i <- 1
  for ( h in 1:length(event.code) )
  {
    irr[i] <- object$ir[covar.code[1], event.code[h]]/object$ir[covar.code[2],event.code[h]]
    var[i] <- object$var[covar.code[1],event.code[h]]/(object$ir[covar.code[2],event.code[h]]^2) + object$var[covar.code[2],event.code[h]]*(object$ir[covar.code[1],event.code[h]]^2)/(object$ir[covar.code[2],event.code[h]]^4)
    
    if ( ci.fun == "lin" ){
      conf.lower[i] <- irr[i]-qnorm(quant)*sqrt(var[i])
      conf.upper[i] <- irr[i]+qnorm(quant)*sqrt(var[i])
    }
    if ( ci.fun == "log" ){
      conf.lower[i] <- irr[i]*exp(-qnorm(quant)*sqrt(var[i])/irr[i])
      conf.upper[i] <- irr[i]*exp(+qnorm(quant)*sqrt(var[i])/irr[i])
    }

    
    i <- i+1
  }
  res <- list (
               covar.code=covar.code,
               covar.lab=covar.lab,
               irr=irr,
               var=var,
               conf.lower=conf.lower,
               conf.upper=conf.upper,
               event.code = event.code,
               event.lab = event.lab,
               ci.fun=ci.fun,
               ci.level=ci.level)

  class(res) <- "irates.ratio"
  return ( res )
}

