 
  #########################
  #### dbplsr functions ###
  #########################





    ################################
    #### dbplsr of class formula ###
    ################################

     #generic function with a commun parameter (y).
 dbplsr<-function(...)  UseMethod("dbplsr")
 
 dbplsr.formula <- function(formula,data,...,metric="euclidean",method="ncomp",
                      weights,ncomp) 
{
  # call dbglm
  mf <- match.call(expand.dots = FALSE)

  # control metric. See the auxiliar function
  metric<-control_metric(metric)  
  
  if (missing(data))
    data <- environment(formula)
    
  # recover z and y of the formula
  zy <- formula_to_zy(formula,data,mf,"dblm",metric)

  # y and z are defined--> pass to default method (try for avoid the program crash). 
  try(ans <- dbplsr.yz(z=zy$z,y=zy$y,metric=metric,weights=weights,method=method,
        ncomp=ncomp))  
  
  if (class(ans)=="try-error") 
    return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
  
  # call dbglm
  ans$call <- mf
  
  return(ans)   
}




    ################################
    #### default dbplsr (y,z) ######
    ################################

dbplsr.yz <- function(y,z,metric="euclidean",weights,ncomp,method="ncomp",...)
{                    
 
  call <- match.call(expand.dots = FALSE)
  # See if z or distance matrix is defined by the user.
  #  require(cluster)
  # control metric. See the auxiliar function
  metric<-control_metric(metric)
  
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
    
  # y and Distance are defined--> pass to dist method (try for avoid crash). 
  try(ans <- dbplsr.D2(D2=D2,y=y,weights=weights,ncomp=ncomp,method=method)) 
 
  # y and Distance are defined--> pass to dist method (try for avoid crash). 
  if (class(ans)=="try-error") 
    return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
 
  ans$call <- call

  attr(ans,"metric") <- metric
  attr(ans,"zs") <- z
  attr(ans,"way") <- way
  
  return(ans)
}




    #################################
    ####    dbplsr with Dist or   ###
    ####  dissimilarity distance ####
    #################################


dbplsr.dist <- function(distance,y,...,weights,ncomp=ncomp,method="ncomp")
{    

   call <- match.call(expand.dots = FALSE)                   
   # program controls: distance must be of class D2 dist or dissimilarity.
   if (missing(distance)||is.null(distance))
    stop("distance must be defined")
   
   # dist to D2
   Delta <- disttoD2(distance)     
   
   try(ans<-dbplsr.D2(D2=Delta,y=y,weights=weights,method=method,ncomp=ncomp))
   if (class(ans)=="try-error")
     return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
   
   ans$call <- call
   ans$distance <- distance
   attr(ans,"way")<-"D2"
   return(ans)

 }


    #################################
    ####  dbplsr with D2 distance ###
    #################################
    
dbplsr.D2 <- function(D2,y,...,weights,ncomp=ncomp,method="ncomp")
{  
                     
   if (!any (class(D2)=="D2")) 
     stop("for a dblm.D2 method the class of the distance matrix D2 must be 'D2'")
    
       # control method. See the auxiliar function
   method<-control_method(method,"dblm")  
   
   # another program controls: See the auxiliar function
   y <- as.matrix(y)
   n <- nrow(y)
  
   # program controls: weights
   if (missing(weights)||is.null(weights))
    weights<-rep(1,nrow(as.matrix(y)))
   if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
   if (sum(weights<0)>0)
    stop("Weights array weights is not valid: some weights are negative")
   if (sum(weights)==0)
    stop("Weights array weights is not valid: sum(weights)=0")
   
   # calculing weights, Dw, y0 and rel.gvar
   ori_weights <-weights           # originals weights !!
   weights <- weights/sum(weights) # percent weights !! 
   
   # G: inner products matrix (symmetrical G)
   G <- Gcalc(n,weights,Delta=D2) 
   class(G)<-"Gram"
         
   try(ans<-dbplsr.Gram(G=G,y=y,weights=ori_weights,method=method,ncomp=ncomp))
    ans$call<-match.call(expand.dots = FALSE)
    attr(ans,"way")<-"D2"
    
    return(ans)
}


    ################################
    ####    dbplsr with Gram    ####
    ####       object           ####
    ################################

dbplsr.Gram <- function(G,y,...,weights,ncomp=ncomp,method="ncomp")
{                                       
   
   method <- control_method(method,"dbplsr")  
   
   # another program controls: See the auxiliar function
   weights <- controls_dbplsr(G,weights,ncomp,y)
   ori_weights <-weights           # originals weights !!
   weights <- weights/sum(weights) # percent weights !!   
   Dw   <- diag(weights)          # diagonal matrix with the weights
   Dsqw <- diag(sqrt(weights))
         
   # y 
   y <- as.matrix(y)
   n <- nrow(y)
   y0 <- y - sum(weights*y)       # centered response varaible(y) 
   
   
   # initial G
   G0 <- G
   Gini <- G0
   gvec <- diag(G0)
   gvar <- t(weights)%*%as.matrix(gvec)
   
   # define the fitted and the residuals
   yhat <- vector("list", ncomp+1)
   ytit <- vector("list", ncomp+1)
   names(yhat) <- c("mean(y)",paste(c(1:ncomp),rep("comp",ncomp),sep = " ")) 
   names(ytit) <- c("1-mean(y)",paste(c(1:ncomp),rep("comp",ncomp),sep = " ")) 
   
   # F
   fk <- vector("list", ncomp)
   names(fk)<-paste(rep("comp",ncomp),c(1:ncomp),sep = "") 
   
   # initial fitted and residuals 
   ytit[[1]] <- y - sum(weights*y)
   yhat[[1]] <- sum(weights*y)
    
   Hk <-0
   bk <-array(0,ncomp) 
   ocv <-array(0,ncomp)
   gcv <-array(0,ncomp)
   aic <-array(0,ncomp)
   bic <-array(0,ncomp)
   gvar.iter <-array(0,ncomp)

   
   # iterations
   for (j in 1:ncomp){
   
     fk[[j]] <- G0%*%Dw%*%ytit[[j]]
     fk2 <- as.numeric(t(fk[[j]])%*%Dw%*%fk[[j]])
     
     Pk <- (fk[[j]]%*%t(fk[[j]])%*%Dw)/fk2
     Qk <- diag(rep(1,n))-Pk
     
     bk[j] <- (t(ytit[[j]])%*%Dw%*%fk[[j]])/fk2 
     
     Gk <- Gproduct(fk[[j]],Dw,G0)  #call internal function Gproduct
     G0 <- Gk
        
     ytit[[j+1]] <- ytit[[j]]-fk[[j]]*bk[j]
     yhat[[j+1]] <- yhat[[j]]+fk[[j]]*bk[j] 
      
     Hk<-Hk+Pk 
                   
     # calculing the ordinary cross-validation estimator
     ocv[j] <- sum(weights*((y-yhat[[j+1]])/(1-diag(Hk)))^2)   
     # calculing the generalized cross-validation estimator
     gcv[j] <- sum(weights*(yhat[[j+1]]-y)^2)/(n*(1-mean(diag(Hk)))^2)
     # aic                                         
     rss <- sum(ori_weights*(ytit[[j+1]])^2)/n 
     aic[j] <- 2*j+n*log(rss) 
     # bic
     bic[j] <- n*log(rss)+j*log(n)
     gvar.iter[j] <- weights%*%diag(G0)      
     # gvar.iter[j] <- weights%*%svd(G0)$d     
     # gvar.iter[j] <- weights%*%(eigen(G0)$values)
                              
   }       
   
  	ncomp.opt <- switch(method,
		  "OCV"=which.min(ocv),
		  "GCV"=which.min(gcv),
		  "AIC"=which.min(aic),
		  "BIC"=which.min(bic),
      "ncomp"= ncomp)
   
                                      
   ans<-list(residuals=ytit,fitted.values=yhat,fk=fk,bk=bk,Pk=Pk,ncomp=ncomp,
             ncomp.opt=ncomp.opt,weights=ori_weights,method=method,y=y,H=Hk,
             G0=Gini,Gk=Gk,gvar=gvar,gvec=gvec,gvar.iter=gvar.iter,ocv=ocv,
             gcv=gcv,aic=aic,bic=bic)
          
   class(ans)<-"dbplsr"
   attr(ans,"way")<-"G"                                     
   ans$call <- match.call(expand.dots = FALSE)
   return(ans)
}


 #generic function with a commun paramatre (y).
 #dbplsr<-function(y,...)  
 # UseMethod("dbplsr") 
