fitPGH <- function(x,                          #E 
                   y,                          #Quantum Efficiency, rETR or P
                   normalize=FALSE,            #Should curve be normalized to E (Default=TRUE for modeling Quantum Efficiency)
                   lowerlim=c(0,0,0),          #Lower bounds of parameter estimates (alpha,Beta,Ps)
                   upperlim=c(100,1000,1000),  #Upper bounds of parameter estimates (alpha,Beta,Ps)
                   fitmethod=c("Nelder-Mead")) #Fitting method passed to modFit  
{
  
  #If normalize =T, assign E = 0 to very small number
  if (normalize==T)  x[x==0] <- 1e-9       
  
  #Remove NA values
  ind   <- is.finite(x) & is.finite(y)
  res   <- rep(NA,length(x))
  x     <- x[ind==T]
  y     <- y[ind==T]
  
  #Intitial Parameter Estimates
  if (normalize==T){ 
    alpha <- max(y)
    beta  <- 0
    ps    <- max(x*y)
  }
  if (normalize==F){ 
    PE    <- y/x
    alpha <-  max(PE[is.finite(PE)])
    beta  <- 0
    ps    <- max(y)
  }
  
  #Load the model
  PGH     <- function(p,x) return(data.frame(x = x, y = p[3]*(1-exp(-1*p[1]*x/p[3]))*exp(-1*p[2]*x/p[3])))
  PGH.E   <- function(p,x) return(data.frame(x = x, y = p[3]*(1-exp(-1*p[1]*x/p[3]))*exp(-1*p[2]*x/p[3])/x))
  if (normalize==F) model.1 <- function(p) (y - PGH(p, x)$y)
  if (normalize==T) model.1 <- function(p) (y - PGH.E(p, x)$y)
  
  
  #In case of non-convergence, NAs are returned
  if (class(try(modFit(f = model.1,p = c(alpha,beta,ps),method = fitmethod, 
                       lower=lowerlim,upper=upperlim, 
                       hessian = TRUE),silent=T))=="try-error"){
    fit <- list(alpha=NA,beta=NA,ps=NA,ssr=NA,residuals=rep(NA,c(length(x))))
  }else{
    fit <- modFit(f = model.1,p = c(alpha,beta,ps),method = fitmethod, 
                  lower=lowerlim,upper=upperlim, hessian = TRUE)
    fit <- list(alpha=summary(fit)$par[1,],beta=summary(fit)$par[2,],ps=summary(fit)$par[3,],
                ssr=fit$ssr,residuals=fit$residuals,model="PGH",normalize=normalize)
  }
  
  return(fit)
  
}
