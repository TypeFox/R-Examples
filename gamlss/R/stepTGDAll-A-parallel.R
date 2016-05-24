# trying to make the function parallel
#-------------------------------------------------------------------------------
stepTGDAll.A <- function(object, scope = NULL,
                                 newdata=NULL,
                                 steps=1000,
                            sigma.scope = NULL, 
                               nu.scope = NULL, 
                              tau.scope =  NULL, 
                                 mu.try = TRUE, 
                              sigma.try = TRUE, 
                                 nu.try = TRUE, 
                                tau.try = TRUE, 
                               parallel = c("no", "multicore", "snow"),
                                  ncpus = 1L, 
                                     cl = NULL, 
                          ...)
{
#-------------------------------------------------------------------------------
#--------------- PARALLEL-------------------------------------------------------
#----------------SET UP PART----------------------------------------------------
if (missing(parallel)) 
    parallel <- "no"
    parallel <- match.arg(parallel)
     have_mc <- have_snow <- FALSE
if (parallel != "no" && ncpus > 1L) 
{
  if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
  else if (parallel == "snow") 
    have_snow <- TRUE
  if (!have_mc && !have_snow) 
         ncpus <- 1L
  loadNamespace("parallel")
}
# -------------- finish parallel------------------------------------------------
#------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------
    ## make sure that the object is visible 
     objectAll<- object
# on.exit(rm(objectAll, envir=.GlobalEnv)) #' delete on exit
#    env <- attach(NULL, name="Object_All_Env")
#    assign("objectAll", object, envir=env)
#    on.exit(detach(Object_All_Env))
    All.anova <-list()
    ## if defferent scope is required
    sigma.scope <- if (is.null(sigma.scope)) scope else sigma.scope
       nu.scope <- if (is.null(   nu.scope)) scope else    nu.scope
      tau.scope <- if (is.null(  tau.scope)) scope else   tau.scope  
    # get the mu model
    #-----------------
    cat("---------------------------------------------------", "\n")
     if ("mu" %in% object$par && mu.try==TRUE)
         {
     current.par <- "mu"
         iferror <- try( assign("objectAll", stepTGD(objectAll, scope=scope, newdata=newdata, direction="forward", parameter = "mu",  parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 { 
                 cat("---------------------------------------------------", "\n")
                 cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu forward.", "\n",
                     "Maybe the mu model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))
                 cat("---------------------------------------------------", "\n")    
                  return(objectAll)
                 }
          # saving the anova
          All.anova <-list(mu.anova.for=objectAll$anova)
         }
     #eval(expression(objectAll$PPP<-objectAll$anova), envir=.GlobalEnv)
    # get the sigma model
    #-------------------- 
     if ("sigma" %in% object$par && sigma.try==TRUE)
        {
         cat("---------------------------------------------------", "\n")
          current.par <- "sigma"
        iferror <- try( assign("objectAll", stepTGD(objectAll, scope=sigma.scope, newdata=newdata, direction="forward", parameter = "sigma", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 {
                  cat("---------------------------------------------------", "\n") 
                  cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma forward.", "\n",
                     "Maybe the sigma model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))    
                  cat("---------------------------------------------------", "\n")
                   return(objectAll) 
                 } 
            All.anova$sigma.anova.for<-objectAll$anova      
        }
     # get the nu model
     #------------------ 
     if (  "nu" %in% object$par && nu.try==TRUE)
        {
         cat("---------------------------------------------------", "\n") 
          current.par <- "nu"
        iferror <- try( assign("objectAll", stepTGD(objectAll, scope=nu.scope, newdata=newdata,direction="forward", parameter = "nu",parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 { 
                    cat("---------------------------------------------------", "\n")
                    cat(paste("ERROR: stepGAICAll has failed trying to fit the model for nu forward.", "\n",
                     "Maybe the nu model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))    
                    cat("---------------------------------------------------", "\n")
                   return(objectAll) 
                 }  
            All.anova$nu.anova.for <- objectAll$anova      
        }
     # get the tau model
     #------------------           
     if ( "tau" %in% object$par && tau.try==TRUE)
        {
        cat("---------------------------------------------------", "\n") 
           current.par <- "tau"
        iferror <- try( assign("objectAll", stepTGD(objectAll, scope=tau.scope,newdata=newdata, direction="both", parameter = "tau", parallel = parallel,  ncpus = ncpus, cl = cl, ...)), silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 { 
                   cat("---------------------------------------------------", "\n")
                   cat(paste("ERROR: stepGAICAll has failed trying to fit the model for tau forward.", "\n",
                     "Maybe the tau model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))    
                   cat("---------------------------------------------------", "\n")
                   return(objectAll) 
                 }
            All.anova$tau.anova.for<-objectAll$anova       
        }   
     # get the nu model
     #------------------ 
     if (  "nu" %in% object$par && nu.try==TRUE && current.par!="nu")
        { 
         cat("---------------------------------------------------", "\n") 
          current.par <- "nu"
        iferror <- try(assign("objectAll", stepTGD(objectAll, scope=nu.scope, newdata=newdata, direction="backward", parameter = "nu", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 { 
                 cat("---------------------------------------------------", "\n")
                  cat(paste("ERROR: stepGAICAll has failed trying to fit the model for nu backward.", "\n",
                     "Maybe the nu model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))    
                   cat("---------------------------------------------------", "\n")
                   return(objectAll) 
                 }  
            All.anova$nu.anova.back<-objectAll$anova
        }
     # get the sigma model
     #------------------ 
      if ("sigma" %in% object$par && sigma.try==TRUE && current.par!="sigma")
        {
         cat("---------------------------------------------------", "\n")
          current.par <- "sigma"
        iferror <- try( assign("objectAll", stepTGD(objectAll, scope=sigma.scope, newdata=newdata, direction="backward", parameter = "sigma", parallel = parallel,  ncpus = ncpus, cl = cl, ...)) , silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 { cat("---------------------------------------------------", "\n")
                   cat(paste("ERROR: stepGAICAll has failed trying to fit the model for sigma backward.", "\n",
                     "Maybe the sigma model is too complicated for the data.", "\n",
                     "The model given is the final before the crush. \n"))    
                   cat("---------------------------------------------------", "\n")
                   return(objectAll) 
                 }  
            All.anova$sigma.anova.back <- objectAll$anova
        }
     # get the final mu model
     #------------------ 
      if ("mu" %in% object$par && mu.try==TRUE && current.par!="mu")
         {
          cat("---------------------------------------------------", "\n")  
           current.par <- "mu"
      iferror <- try( assign("objectAll", stepTGD(objectAll, scope=scope, newdata=newdata,  direction="backward",  parameter = "mu", parallel = parallel,  ncpus = ncpus, cl = cl, ...))  , silent = TRUE)
             if (any(class(iferror)%in%"try-error"))
                 { 
                  cat("---------------------------------------------------", "\n") 
                  cat(paste("ERROR: stepGAICAll has failed trying to fit the model for mu backward.", "\n",
                     "Maybe the mu model is too complicated for the data,", "\n",
                     "The model given is the final before the crush. \n"))
                  cat("---------------------------------------------------", "\n")    
                   return(objectAll) 
                 } 
            All.anova$mu.anova.back<-objectAll$anova        
         }
       cat("---------------------------------------------------", "\n")    

objectAll$anovaAll<-All.anova
objectAll
}
