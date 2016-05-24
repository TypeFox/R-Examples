

 ########################
 #### function dbglm ####
 ########################
 
 ## Description:  Weighted generalized distance-based regression. There are 
 ##        4 different calls to the function dbglm:
 ##           dbglm.formula, dbglm.default, dbglm.D2 and dbglm.dist 
 ##
 ##    Inputs: family: binomial()
 ##                    gaussian()
 ##                    Gamma()
 ##                    inverse.gaussian()
 ##                    poisson()
 ##                    quasi()
 ##                    quasibinomial()
 ##                    quasipoisson()
 ##            maxiter: 100 iterations
 ##            eps1: stopping criterion 1, "Devstat"
 ##            eps2: stopping criterion 2, "mustat"
 ##            metric: "euclidean", "gower", "manhattan"
 ##            weights: Default (all individuals have the same weights)
 ##            mustart: started values: default NULL
 ##            rel.gvar: relative variability used in the linear predictions.
 ##    Outputs: an object of class dbglm.
 ##  
  


    ################################
    #### dbglm of class formula ####
    ################################
    
 #generic function with a commun parameter (y).
 dbglm<-function(...)  UseMethod("dbglm")
 
 dbglm.formula <-function(formula,data,family=gaussian, method ="GCV", full.search=TRUE,...,metric="euclidean",
            weights,maxiter=100,eps1=1e-10,eps2=1e-10,rel.gvar=0.95,eff.rank=NULL,
            offset,mustart=NULL, range.eff.rank ) 
{
  # call dbglm
  mf <- match.call(expand.dots = FALSE)
  
  if (missing(range.eff.rank))
   range.eff.rank <- NULL
   
  
  # control metric. See the auxiliar function
  metric<-control_metric(metric)    
  
  if (missing(data))
    data <- environment(formula)
    
  # recover z and y of the formula
  zy <- formula_to_zy(formula,data,mf,"dbglm",metric)
  
  # y and z are defined--> pass to default method (try for avoid the program crash).
  try(ans <- dbglm.yz(y=zy$y,z=zy$z,metric=metric,family=family, method = method, full.search =full.search,
            weights=weights,maxiter=maxiter,eps1=eps1,eps2=eps2,rel.gvar=rel.gvar,eff.rank=eff.rank,
            offset=offset,mustart=mustart,range.eff.rank=range.eff.rank))        
 
  if (class(ans)[1]=="try-error") 
   return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
  
  # hidden attribute 'call' (used in print and summary)
  ans$call <- mf
  attr(ans,"zs") <- zy$zini
  return(ans)   
}



    ################################
    #### default dbglm (y,z) #######
    ################################

dbglm.yz <- function(y,z,family=gaussian,metric="euclidean", method ="GCV",
        full.search=TRUE, weights,maxiter=100,eps1=1e-10,eps2=1e-10,
        rel.gvar=0.95,eff.rank=NULL,offset = rep(0, length(y)),
        mustart=NULL, range.eff.rank,...)
{
  # See if z or distance matrix is defined by the user.
  # require(cluster)

   # control metric. See the auxiliar function
   metric <- control_metric(metric)
   
   if (missing(range.eff.rank))
    range.eff.rank <- NULL
    
   # call z_to_dist to pass the explanatory variables to an object of class dist 
   dist_and_way <- z_to_dist(z,metric)
   D <- dist_and_way$D
   way <- dist_and_way$way
  
   # if metric=gower. the distance matrix D is already the squared.
   if (metric=="gower"){
     D2 <-as.matrix(D)
     class(D2) <- "D2"
   }else
     D2 <-disttoD2(D)
     
   # y and Distance are defined--> pass to dist method (try for avoid the program crash). 
   try(ans<-dbglm.D2(y=y,D2=D2,family=family, method =method, full.search=full.search, weights=weights,maxiter=maxiter,
          eps1=eps1,eps2=eps2,eff.rank=eff.rank,rel.gvar=rel.gvar,
          offset = offset,mustart=mustart,range.eff.rank=range.eff.rank))     
    
   if (class(ans)[1]=="try-error") 
    return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
 
   # hidden attributes 'metric', 'call', 'zs' and 'way'  
   attr(ans,"metric")<-metric
   attr(ans,"zs")<-z
   ans$call  <- match.call(expand.dots = FALSE)
   attr(ans,"way")<-way
  
   return(ans)
}


    #################################
    ####    dbglm with Dist or   ####
    ####  dissimilarity distance ####
    #################################

dbglm.dist <- function(distance,y,family=gaussian, method = "GCV", full.search=TRUE,weights,maxiter=100,
                eps1=1e-10,eps2=1e-10,rel.gvar=0.95,eff.rank=NULL, offset,mustart=NULL, range.eff.rank,...)
{ 
  
   # program controls: distance must be of class D2 dist or dissimilarity.
   if (missing(distance)||is.null(distance))
    stop("distance must be defined")
    
   if (missing(range.eff.rank))
    range.eff.rank <- NULL
   
   # dist to D2
   # dist to D2
   Delta <- disttoD2(distance)     
  
   try(ans <- dbglm.D2(D2=Delta,y=y,family=family, method =method, full.search=full.search, weights=weights,maxiter=maxiter,
         eps1=eps1,eps2=eps2,rel.gvar=rel.gvar,eff.rank=eff.rank,
         offset = offset,mustart=mustart,range.eff.rank=range.eff.rank))

  if (class(ans)[1]=="try-error")
   return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
   
  # hidden attributes 'call' and 'way'   
  ans$call <-match.call(expand.dots = FALSE)
  attr(ans,"way")<-"D2"
  attr(ans,"zs")<-Delta

   return(ans)                             
}


    #################################
    ####  dbglm with D2 distance ####
    #################################

dbglm.D2<-function(D2,y,...,family=gaussian, method ="GCV", full.search=TRUE, weights,maxiter=100,eps1=1e-10,
              eps2=1e-10,rel.gvar=0.95,eff.rank=NULL,offset,mustart=NULL, range.eff.rank)
{ 
  # D2 squared distances must be of class "D2"
  if (!any (class(D2)=="D2")) 
    stop("for a dbglm.D2 method the class of the distance matrix D2 must be 'D2'")

   Delta <- D2
   # number of observations
   nobs <- nrow(as.matrix(y)) 
   n <- nobs 

   # Program controls: see de auxiliar function 
   controls <- controls_dbglm(Delta,weights,offset,rel.gvar,maxiter,
                        eps1,eps2,y, method)
   
   method <-control_method(method,"dbglm")  
   
   if (missing(range.eff.rank))
    range.eff.rank <- NULL
   
   weights  <- controls$weights 
   #weights <- weights/sum(weights)
   offset   <- controls$offset 
   rel.gvar <- controls$rel.gvar 
   maxiter  <- controls$maxiter 
   eps1 <- controls$eps1
   eps2 <- controls$eps2
       
   # start and etastart <- NULL (internal parameters of glm that need initialize
   start<-NULL
   etastart<-NULL
  
   # get the function of the family (can be a character or function)  
   family <-control_family(family)

   # functions variance and linkinv according to the family of y
   variance <- family$variance
   linkinv <- family$linkinv
   
   # control: variance and linkinv must be a functions. 
   if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
        call. = FALSE)
   
   # another important functions: dev.resids, aic, mu.eta, valideta and validmu     
   dev.resids <- family$dev.resids # deviance
   aic <- family$aic               # aic
   mu.eta <- family$mu.eta         # wdmu: derivative function(eta)

   valideta <-family$valideta      # check if starting values mu and eta are valid 
   validmu <- family$validmu  
   
   # possibles links
   problem.links <- c("inverse","identity","logit","1/mu^2") 
      
   # Iterative dblm:  Distance-Based Linear Model (with weights)
   # Step 0: Initial values, assigned by the user, or set the intrinsic values 
   # in the initialize function
   class(Delta)<-"D2"
   if (is.null(mustart)) {
        eval(family$initialize)
   }
   else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
   }
   
   # eta = the link function.
   eta<-family$linkfun(mustart)
   
   # mu= the inverse of the link function.
   mu <- linkinv(eta)
   
   # check if starting values mu and eta are valid with validmu and valideta 
   # functions. 
   if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
   
   ### select eff.rank new!!!
   GCV           <- NULL
   BIC           <- NULL
   AIC           <- NULL        
   old.criterion <- NULL

 
   dev.resids <- family$dev.resids # deviance
   aic <- family$aic               # aic
   mu.eta <- family$mu.eta         # wdmu: derivative function(eta)

   valideta <-family$valideta      # check if starting values mu and eta are valid 
   validmu <- family$validmu  
   call    <- match.call(expand.dots = FALSE)
     
   if (all(method != c("eff.rank","rel.gvar")))
   {
     if (is.null(range.eff.rank))
              stop("'range.eff.rank' must be specified")
     else 
     {
       if(length(range.eff.rank)!=2 || range.eff.rank[2]- range.eff.rank[1]<0)
            stop("'range.eff.rank' is not correctly specified")
     }   
              
    GCV           <- rep(0,range.eff.rank[2]- range.eff.rank[1] + 1)
    BIC           <- rep(0,range.eff.rank[2]- range.eff.rank[1] + 1)
    AIC           <- rep(0,range.eff.rank[2]- range.eff.rank[1] + 1)        
    old.criterion <- Inf
      
     if(full.search)
     { 
        for (eff.rank in range.eff.rank[1]:range.eff.rank[2])
        {
           dbglm.aux    <- dbglm.iteration(y = y, mu = mu, weights = weights, nobs = nobs, eta = eta, Delta=Delta, method=method, 
                                           offset = offset, n = n, eff.rank = eff.rank, rel.gvar = rel.gvar, dev.resids=dev.resids,
                                           aic=aic, mu.eta = mu.eta, valideta=valideta, validmu=validmu, family=family,
                                           variance=variance,linkinv=linkinv, problem.links=problem.links, eps1=eps1,eps2=eps2,maxiter=maxiter)
           
          if (method == "GCV" && dbglm.aux$gcv.model < old.criterion)
          {
             dbglm.opt     <- dbglm.aux
             old.criterion <- dbglm.aux$gcv.model 
          }
          if (method == "AIC" && dbglm.aux$aic < old.criterion)
          {
             dbglm.opt     <- dbglm.aux
             old.criterion <- dbglm.aux$aic.model 
          }
          if (method == "BIC" && dbglm.aux$bic.model < old.criterion)
          {
             dbglm.opt     <- dbglm.aux
             old.criterion <- dbglm.aux$bic.model 
          }
           
           AIC[eff.rank] <- dbglm.aux$aic.model
           BIC[eff.rank] <- dbglm.aux$bic.model 
           GCV[eff.rank] <- dbglm.aux$gcv.model
        }
     }
     else
     {
        f <- function(rk, y, mu, weights,  nobs, eta, Delta, method, offset, n, rel.gvar, dev.resids,
                                         aic, mu.eta, valideta, validmu, family,variance, linkinv, problem.links,
                                         eps1,eps2,maxiter)
        {
       
           dbglm.aux    <- dbglm.iteration(y = y, mu = mu, weights = weights, nobs = nobs, eta = eta, Delta=Delta, method=method, 
                                           offset = offset, n = n, eff.rank = rk, rel.gvar = rel.gvar, dev.resids=dev.resids,
                                           aic=aic, mu.eta = mu.eta, valideta=valideta, validmu=validmu, family=family,
                                           variance=variance,linkinv=linkinv, problem.links=problem.links, eps1=eps1,eps2=eps2,maxiter=maxiter)
           
          if (method == "GCV")
          {
             criterion <- dbglm.aux$gcv.model 
          }
          if (method == "AIC")
          {
             criterion <- dbglm.aux$aic.model 
          }
          if (method == "BIC")
          {
             criterion <- dbglm.aux$bic.model 
          }
          return(criterion)
      }         
      criterion.opt <- optimise(f = f, range.eff.rank ,y = y, mu = mu, weights = weights, nobs = nobs, eta = eta, Delta=Delta, method=method, 
                         offset = offset, n = n,  rel.gvar = rel.gvar, dev.resids=dev.resids,
                         aic=aic, mu.eta = mu.eta, valideta=valideta, validmu=validmu, family=family,
                         variance=variance,linkinv=linkinv, problem.links=problem.links, eps1=eps1,eps2=eps2,maxiter=maxiter)
    
      eff.rank  <- round(criterion.opt$minimum)
      criterion <- criterion.opt$objective

      dbglm.opt    <- dbglm.iteration(y = y, mu = mu, weights = weights,  nobs = nobs, eta = eta, Delta=Delta, method=method, 
                                         offset = offset, n = n, eff.rank = eff.rank, rel.gvar = rel.gvar, dev.resids=dev.resids,
                                         aic=aic, mu.eta = mu.eta, valideta=valideta, validmu=validmu, family=family,
                                         variance=variance,linkinv=linkinv, problem.links=problem.links, eps1=eps1,eps2=eps2,maxiter=maxiter)

     }   
   } 
   else
   {
         dbglm.opt    <- dbglm.iteration(y = y, mu = mu, weights = weights,  nobs = nobs, eta = eta, Delta=Delta, method=method, 
                                         offset = offset, n = n, eff.rank = eff.rank, rel.gvar = rel.gvar, dev.resids=dev.resids,
                                         aic=aic, mu.eta = mu.eta, valideta=valideta, validmu=validmu, family=family,
                                         variance=variance,linkinv=linkinv, problem.links=problem.links, eps1=eps1,eps2=eps2,maxiter=maxiter)
   }      
   
  
   # return a list with the following attributes 
   ans <- dbglm.opt
   ans$call <- call                      
   class(ans) <- c(ans$class, c("dbglm", "dblm"))
   attr(ans,"last_dblm")<- dbglm.opt$dblm_aux
   attr(ans,"ori_weights")<-weights
   attr(ans,"way")<-"D2"
   attr(ans,"G")<-D2toG(Delta,weights)
   attr(ans,"zs")<-Delta
   attr(ans,"eta")<-eta  
   attr(ans,"aics")<-AIC  
   attr(ans,"bics")<-BIC  
   attr(ans,"gcvs")<-GCV  
   attr(ans,"full.search")<-full.search
   attr(ans,"method")<-method
   attr(ans,"range.eff.rank")<-range.eff.rank
   
    
  return(ans)
}



    #################################
    ####  dbglm with G products  ####
    #################################

dbglm.Gram<-function(G,y,...,family=gaussian, method ="GCV", full.search=TRUE, weights,maxiter=100,eps1=1e-10,
              eps2=1e-10,rel.gvar=0.95,eff.rank=NULL,offset,mustart=NULL, range.eff.rank = c(1,ncol(G)-1))
{ 
  # D2 squared distances must be of class "D2"
  if (class(G)[1]!="Gram") 
    stop("for a dbglm.Gram method the class of the inner products matrix G must be 'Gram'")
  else{
  
   D2<-GtoD2(G)
   
   try(ans <- dbglm.D2(D2=D2,y=y,family=family, method = method, full.search=full.search, weights=weights,maxiter=maxiter,
         eps1=eps1,eps2=eps2,rel.gvar=rel.gvar,eff.rank=eff.rank,
         offset = offset,mustart=mustart,range.eff.rank=range.eff.rank))
  }
  if (class(ans)[1]=="try-error")
   return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
  
  # hidden attributes 'call' and 'way'   
  ans$call <-match.call(expand.dots = FALSE)
  attr(ans,"way")<-"G"
  ans$G <- G    
  return(ans)
}







