vario.fit <- function (vario, bins, weights=rep(1, length(vario)),
                       type=c("spherical", "gaussian", "nugget", "linear",
                              "exponential", "sill", "periodic", "hole"),
                       start.vals=list(c0=0, 
                                       c1=max(vario), 
                                       a=max(bins)/4,
                                       b=0.1,
                                       c=0.1), control=list(maxit=10000)) {
  
  locs=which(!is.na(vario))
  data=data.frame(vario=vario[locs], bins=bins[locs])
  weights=weights[locs]
  types=c("spherical", "gaussian", "nugget", "linear", "exponential", "sill", 
          "periodic", "hole")
  type=match.arg(tolower(type), types)
  
  ## Linear model used to determine initial parameters for finicky NLS
  vario.lin=lm(vario ~ bins, weights=weights, data=data)
  
  if (type=="nugget") {
    vario.mod=lm(vario ~ 1, weights=weights, data=data)
    names=c("nugget")    
    success=TRUE
    vario.mod$convergence=0
  }
  else if (type=="linear") {
    names=c("c0", "b")
    vario.mod=try(nls(vario ~ c0+b*bins, 
                      weights=weights, 
                      lower=c(c0=0, b=-Inf),
                      algorithm="port",
                      start=list(c0=start.vals$c0,
                                 b=coef(vario.lin)[2]), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        coef(vario.lin)[2]),
                      rmse.lin, weights=weights, control=control,
                      data=data)
    }
    else{
      success=TRUE
    }
  }
  else if (type=="sill") {
    names=c("c0", "c1", "a", "b")
    vario.mod=try(nls(vario ~ (bins <= a)*(c0+b*bins) + (bins > a)*(c0+c1), 
                      weights=weights, 
                      lower=c(c0=0, c1=0, a=0, b=-Inf),
                      algorithm="port",                      
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a,
                                 b=coef(vario.lin)[2]), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a,
                        coef(vario.lin)[2]),
                        rmse.sill, weights=weights, control=control,
                      data=data)
    }
    else {
      success=TRUE
    }
  }
  else if (type=="exponential") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ c0+c1*(1-exp(-bins/a)),
                      weights=weights, 
                      lower=c(c0=0, c1=0, a=0),
                      algorithm="port",
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)
  
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a),
                        rmse.expo, weights=weights, control=control,
                      data=data)
    }
    else {
      success=TRUE  
    }
  }
  else if (type=="spherical") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ (bins <= a)*(c0+c1*((3*bins)/(2*a)-(1/2)*(bins/a)^3))+(bins > a)*(c0+c1),
                      weights=weights, 
                      lower=c(c0=0, c1=0, a=0),
                      algorithm="port",
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a),
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a),
                        rmse.sphere, weights=weights, control=control,
                      data=data)
    }
    else {
      success=TRUE
    }
  }
  else if (type=="gaussian") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ c0+c1*(1-exp(-3*(bins^2)/(a^2))), 
                      weights=weights, 
                      lower=c(c0=0, c1=0, a=0),
                      algorithm="port",                      
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)

    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a),
                        rmse.gauss, weights=weights, control=control,
                      data=data)
    }
    else {
      success=TRUE
    }
  }
  else if (type=="periodic") {
    names=c("a", "b", "c")
    vario.mod=try(nls(vario ~ a*cos(b*pi/max(bins)+c)*bins, 
                      weights=weights,
                      start=list(a=start.vals$a,
                                 b=start.vals$b,
                                 c=start.vals$c), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$a, 
                        start.vals$b,
                        start.vals$c),
                        rmse.period, weights=weights, control=control,
                      data=data)
    }
    else {
      success=TRUE
    }  
  }
  else if (type=="hole") {
    names=c("c0", "c1", "a")
    vario.mod=try(nls(vario ~ c0+c1*(1-(a*sin(bins/a))/bins), 
                      weights=weights, 
                      lower=c(c0=0, c1=0, a=0),
                      algorithm="port",                      
                      start=list(c0=start.vals$c0,
                                 c1=start.vals$c1,
                                 a=start.vals$a), 
                      data=data), silent=TRUE)
    if (class(vario.mod)=="try-error") {
      success=FALSE
      vario.mod=optim(c(start.vals$c0, 
                        start.vals$c1,
                        start.vals$a),
                      rmse.hole, weights=weights, control=control,
                      data=data)
    }
    else {
      success=TRUE
    }  
  }
  opt=vario.stats(data, vario.mod, type, names, success)
  converge=ifelse(opt$convergence==0, TRUE, FALSE)
  results=list(vario=vario, bins=bins, AIC=opt$AIC, rmse=opt$rmse, 
               params=opt$params, fit=opt$fit, model=type, nls.success=opt$nls.success,
               convergence=converge)
  class(results)="variofit"
  return (results)
}

vario.stats <- function (data, opt, type, names, success) {
  vario=data$vario
  bins=data$bins
  N=length(vario)
  if (success) {
    fit=predict(opt)
    params=coef(opt)
    rmse=sqrt(mean((predict(opt)-vario)^2))
  }
  else {
    names(opt$par)=names
    params=opt$par
    rmse=opt$value
    
    if (type=="sill") {
      fit=ifelse (bins <= opt$par["a"], 
                  opt$par["c0"]+opt$par["b"]*bins,
                  opt$par["c0"]+opt$par["c1"])
    }
    else if (type=="linear") {
      fit=opt$par["c0"]+opt$par["b"]*bins
    }
    else if (type=="exponential") {
      fit=opt$par["c0"]+opt$par["c1"]*(1-exp(-bins/opt$par["a"]))
    }
    else if (type=="spherical") {
      fit=ifelse (bins <= opt$par["a"],
                  opt$par["c0"]+
                    opt$par["c1"]*((3*bins)/(2*opt$par["a"])-
                                        0.5*(bins/opt$par["a"])^3),
                  opt$par["c0"]+opt$par["c1"])
    }
    else if (type=="gaussian") {
      fit=opt$par["c0"]+opt$par["c1"]*(1-exp(-3*(bins)^2/(opt$par["a"]^2)))
    }
    else if (type=="periodic") {
      fit=opt$par["a"]*cos(opt$par["b"]*pi*(bins/max(bins))+opt$par["c"])
    }
    else if (type=="hole") {
      fit=opt$par["c0"]+opt$par["c1"]*(1-(opt$par["a"]*sin(bins/opt$par["a"]))/bins)      
    }    
  }
  
  mod.aic=N*log(rmse^2)+2*(length(params))
  names(params)=names
  return (list(AIC=mod.aic, rmse=rmse, params=params, fit=fit, nls.success=success, 
               convergence=opt$convergence))
} 

rmse.lin <- function (x, weights, data) {
  c0=x[1]; b=x[2]
  vario=data$vario; bins=data$bins
  
  if (c0 >= 0) {
    variohat=c0+b*bins 
    rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
  }
  else
    rmse=Inf
}

rmse.period <- function (x, weights, data) {
  a=x[1]; b=x[2]; c=x[3]
  vario=data$vario; bins=data$bins
  
  variohat=a*cos(b*pi*(bins/max(bins))+c)
  rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
}

rmse.hole <- function (x, weights, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a >= 0 & a <= max(bins) & c0 >= 0 & c1 >= 0) {
    variohat=c0+c1*(1-(a*sin(bins/a))/bins)
    rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
  }
  else
    rmse=Inf
}

rmse.sill <- function (x, weights, data) {
  c0=x[1]; c1=x[2]; a=x[3]; b=x[4]
  vario=data$vario; bins=data$bins
  
  if (a <= max(bins) & a >= 0 & b >= 0 & c1 >= 0 & c0 >= 0) {
    variohat=ifelse (bins <= a, 
                     c0+b*bins, 
                     c0+c1)
    rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
  }
  else
    rmse=Inf
}

rmse.expo <- function (x, weights, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a <= max(bins) & a >= 0 & c1 >= 0 & c0 >= 0) {
    variohat=c0+c1*(1-exp(-bins/a))
    rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
  }
  else
    rmse=Inf
}

rmse.sphere <- function (x, weights, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a <= max(bins) & a >= 0 & c1 >= 0 & c0 >= 0) {
    variohat=ifelse (bins < a, 
                     c0+c1*(3*bins/(2*a)-0.5*(bins/a)^3), 
                     c0+c1)
    rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
  }
  else
    rmse=Inf
}

rmse.gauss <- function (x, weights, data) {
  c0=x[1]; c1=x[2]; a=x[3]
  vario=data$vario; bins=data$bins
  
  if (a <= max(bins) & a >= 0 & c1 >= 0 & c0 >= 0) {
    variohat=c0+c1*(1-exp(-3*bins^2/(a^2)))
    rmse=sqrt(weighted.mean((vario-variohat)^2, weights))
  }
  else
    rmse=Inf
}
