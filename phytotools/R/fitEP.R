fitEP  <- function(x,                          #E 
                   y,                          #Quantum Efficiency, rETR or P
                   normalize=FALSE,            #Should curve be normalized to E (Default=TRUE for modeling Quantum Efficiency)
                   lowerlim=c(0,0,0),          #Lower bounds of parameter estimates (alpha,Eopt,Ps)
                   upperlim=c(100,2000,2000),  #Upper bounds of parameter estimates (alpha,Eopt,Ps)
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
    eopt  <- mean(x)
    ps    <- max(x*y)
  }
  if (normalize==F){ 
    PE    <- y/x
    alpha <- max(PE[is.finite(PE)])
    eopt  <- mean(x)
    ps    <- max(y)
  }
  
  #Load the model
  EP     <- function(p,x) return(data.frame(x = x, y = x/((1/(p[1]*p[2]^2))*x^2+(1/p[3]-2/(p[1]*p[2]))*x+(1/p[1]))))
  EP.E   <- function(p,x) return(data.frame(x = x, y = 1/((1/(p[1]*p[2]^2))*x^2+(1/p[3]-2/(p[1]*p[2]))*x+(1/p[1]))))
  if (normalize==F) model.1 <- function(p) (y - EP(p, x)$y)
  if (normalize==T) model.1 <- function(p) (y - EP.E(p, x)$y)
  
  
  #In case of non-convergence, NAs are returned
  if (class(try(modFit(f = model.1,p = c(alpha,eopt,ps), method = fitmethod, 
                       lower=lowerlim, upper=upperlim, 
                       hessian = TRUE),silent=T))=="try-error"){
    fit <- list(alpha=NA,eopt=NA,ps=NA,ssr=NA,residuals=rep(NA,c(length(x))))
  }else{
    fit <- modFit(f = model.1,p = c(alpha,eopt,ps),method = fitmethod, 
                  lower=lowerlim, upper=upperlim, hessian = TRUE)
    fit <- list(alpha=summary(fit)$par[1,],eopt=summary(fit)$par[2,],ps=summary(fit)$par[3,],
                ssr=fit$ssr,residuals=fit$residuals,model="EP",normalize=normalize)
  }
  
  return(fit)
  
}
