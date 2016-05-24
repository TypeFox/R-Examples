fitWebb <- function(x,                          #E 
                    y,                          #Quantum Efficiency, rETR or P
                    normalize=FALSE,            #Should curve be normalized to E (Default=FALSE)
                    lowerlim=c(0,1),            #Lower bounds of parameter estimates (alpha,Ek)
                    upperlim=c(100,1000),       #Upper bounds of parameter estimates (alpha,Ek)
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
    ek    <- mean(range(x))
  }
  if (normalize==F){ 
    PE    <- y/x
    alpha <-  max(PE[is.finite(PE)])
    ek    <- mean(range(x))
  }
  
  #Load the model
  Webb   <- function(p,x) return(data.frame(x = x, y = p[1]*p[2]*(1-exp(-1*x/p[2]))))  
  Webb.E <- function(p,x) return(data.frame(x = x, y = p[1]*p[2]*(1-exp(-1*x/p[2]))/x))
  if (normalize==F) model.1 <- function(p) (y - Webb(p, x)$y)
  if (normalize==T) model.1 <- function(p) (y - Webb.E(p, x)$y)
  
  
  #In case of non-convergence, NAs are returned
  if (class(try(modFit(f = model.1,p = c(alpha,ek),method = fitmethod, 
                       lower=lowerlim,upper=upperlim, 
                       hessian = TRUE),silent=T))=="try-error"){
    fit <- list(alpha=NA,ek=NA,ssr=NA,residuals=rep(NA,c(length(x))))
  }else{
    fit <- modFit(f = model.1,p = c(alpha,ek),method = fitmethod, 
                  lower=lowerlim,upper=upperlim, hessian = TRUE)
    fit <- list(alpha=summary(fit)$par[1,],ek=summary(fit)$par[2,],
                ssr=fit$ssr,residuals=fit$residuals,model="Webb",normalize=normalize)
  }
  
  return(fit)
  
}
