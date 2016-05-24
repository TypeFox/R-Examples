tegarch <-
function(y, asym=TRUE, skew=TRUE, components=1,
  initial.values=NULL, lower=NULL, upper=NULL, hessian=TRUE,
  lambda.initial=NULL, c.code=TRUE, logl.penalty=NULL,
  aux=NULL, ...)
{
y <- as.zoo(y)
y <- na.trim(y)
y.index <- index(y)
y <- coredata(y)

#xts related:
if(is.matrix(y)){
  if (NCOL(y) > 1)
    stop("Dependent variable not a 1-dimensional matrix")
  y <- y[, 1]
}
y <- as.numeric(y)

aux$asym <- asym
aux$skew <- skew
aux$iN <- length(y)
aux$signnegy <- sign(-y)
aux$u <- rep(0,aux$iN)

if(components==1){
  if(is.null(initial.values)){
    initial.values <- c(0.02,0.95,0.05,0.01,10,0.98)
  }
  if(is.null(lower)){
    lower <- c(-Inf,-1+.Machine$double.eps,-Inf,-Inf,
      2+.Machine$double.eps,.Machine$double.eps)
  }
  if(is.null(upper)){
    upper <- c(Inf,1-.Machine$double.eps,Inf,Inf,Inf,Inf)
  }
  if(!aux$skew){
    initial.values <- initial.values[-6]
    lower <- lower[-6]
    upper <- upper[-6]
  }
  if(!aux$asym){
    initial.values <- initial.values[-4]
    lower <- lower[-4]
    upper <- upper[-4]
  }
  if(is.null(logl.penalty)){
    logl.penalty <- tegarchLogl(y, initial.values,
      lower=lower, upper=upper, lambda.initial=lambda.initial,
      logl.penalty=-1e+100, c.code=c.code, aux=aux)
  }
  objective.f <- function(pars, x=y){f <- -tegarchLogl(x,
    pars, lower=lower, upper=upper,
    lambda.initial=lambda.initial, logl.penalty=logl.penalty,
    c.code=c.code, aux=aux); f}
}else{
  if(is.null(initial.values)){
    initial.values <- c(0.02,0.95,0.9,0.001,0.01,0.005,10,0.98)
  }
  if(is.null(lower)){
    lower <- c(-Inf,-1+.Machine$double.eps,
      -1+.Machine$double.eps,-Inf,-Inf,-Inf,
      2+.Machine$double.eps,.Machine$double.eps)
  }
  if(is.null(upper)){
  upper <- c(Inf,1-.Machine$double.eps,
    1-.Machine$double.eps,Inf,Inf,Inf,Inf,Inf)
  }
  if(!aux$skew){
    initial.values <- initial.values[-8]
    lower <- lower[-8]
    upper <- upper[-8]
  }
  if(!aux$asym)
    stop("For identification, asym=TRUE is required when components=2")
#  if(!aux$asym){
#    asym=TRUE
#  }
  if(is.null(logl.penalty)){
    logl.penalty <- tegarchLogl2(y, initial.values,
      lower=lower, upper=upper, lambda.initial=lambda.initial,
      logl.penalty=-1e+100, c.code=c.code, aux=aux)
  }
  objective.f <- function(pars, x=y){f <- -tegarchLogl2(x,
    pars, lower=lower, upper=upper, lambda.initial=lambda.initial,
    logl.penalty=logl.penalty, c.code=c.code, aux=aux); f}
}

#estimate:
est <- nlminb(initial.values, objective.f, lower=lower,
  upper=upper, x=y, ...)
est$objective <- -est$objective
#sic <- -2*est$objective/aux$iN + length(est$par)*log(aux$iN)/aux$iN
#est <- c(list(sic=sic), est)

#compute Hessian:
if(hessian){
  hessian.numeric <- -optimHess(est$par, objective.f)
  est <- c(list(hessian=hessian.numeric), est)
}

#type of model:
model <- c(components, asym, skew)
names(model) <- c("components", "asym", "skew")
est <- c(list(y=zoo(y, order.by=y.index), date=date(),
  initial.values=initial.values, lower=lower, upper=upper,
  lambda.initial=lambda.initial, model=model), est)

#add names:
if(components==1){
  parnames <- c("omega", "phi1", "kappa1", "kappastar", "df", "skew")
  if(!aux$skew){ parnames <- parnames[-6] }
  if(!aux$asym){ parnames <- parnames[-4] }
}else{
  parnames <- c("omega", "phi1", "phi2", "kappa1", "kappa2",
    "kappastar", "df", "skew")
  if(!aux$skew){ parnames <- parnames[-8] }
#  if(!aux$asym){
#    est$NOTE <- "2 comp spec without leverage/asymmetry not available"
#  }
}
names(est$par) <- parnames
if(hessian){
  colnames(est$hessian) <- parnames
  rownames(est$hessian) <- parnames
}
names(est$initial.values) <- parnames
names(est$lower) <- parnames
names(est$upper) <- parnames

#out:
class(est) <- "tegarch"
return(est)
}
