
# Function of estimation
#----------------------------------------------------------------------------------------

flexrsurv.glm.fit <- function(formula, data, 
                           Spline, degree.Bh=3, knots.Bh, log.Bh=FALSE, model, 
                           control=list(epsilon = 1e-8, maxit = 100, trace = TRUE),
                           Min_T, Max_T, name.runningtime=".t", 
                           start=NULL) {

  control <- do.call("glm.control", control)
  

  
  tik <- data$tik
  data$exp_nbevent <-  tik*data$rate
    
  
  # Definition of the link function
  #----------------------------------------------------------------------------------------
  
  linkfun  <- function(mu){

    #log(mu - data$exp_nbevent)
    ifelse(mu < data$exp_nbevent,  - log(.Machine$double.xmax)/2, log(mu - data$exp_nbevent))
  }
  # inverse of the link function
  linkinv  <- function(eta){
    pmax(exp(eta), .Machine$double.eps/2) + data$exp_nbevent 
  }
  # dmu/deta
  mu.eta   <- function(eta){
    pmax(exp(eta), .Machine$double.eps/2)
  }
  valideta <- function(eta){
    TRUE
  }
  flexrsurv.link <- structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                              valideta = valideta, name = "flexrsurv"), class = "link-glm")
  
 
  
  
  if (length(knots.Bh)>=2) {
    k<-(knots.Bh)
    for (l in 2:length(k)) {
      k[l] <- paste(k[l-1],k[l],sep=",") 
    } 
    knots.Bh <- paste("c(",k[length(k)],")",sep="")  # consider multiple knots
  }
  
  gamma0_t <- as.character(paste("NLL(", name.runningtime , ', Spline="', Spline,'", Knots=', knots.Bh,
                                 ", Degree=", degree.Bh, ", Log=", log.Bh,
                                 ", Intercept=TRUE, Boundary.knots=c(", Min_T, ", ", Max_T, "))", sep=""))
 
  glmformula <- as.formula(paste(".fail ~", gamma0_t, "+", 
                                     deparse(terms(formula)[[3]], width.cutoff = 500L),
                                     "+ offset(log(tik)) -1", sep =" ", collapse = NULL))  

  estim <- glm(formula=glmformula, family=poisson(link=flexrsurv.link), data=data, control=control, start=start)

  expetahat <- pmax(exp(predict(estim)), .Machine$double.eps)

  ll <- sum(-expetahat + data$.fail*log(expetahat/data$tik + data$rate))


   cat("ll in GLM init: ",format(ll, digit=15),"\n")
  estim$loglik <- ll
      
  estim$method <- "flexrsurv.glm.fit"


  estim$df.null <- NULL
  estim$deviance <- NULL
  estim$aic <- NULL
  estim$weights <- NULL     
  estim$prior.weights <- NULL 

  
  class(estim) <- c("flexrsurv.glm", class(estim))

  return(estim)  
} 

