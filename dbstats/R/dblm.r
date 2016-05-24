

 ##########################
 ####  dblm functions  ####
 ##########################

 ## Description: 
 ##    weighted distance-based Regression. Generic function with 4 different 
 ##    ways to compute the dblm: 
 ##        dblm.formula: the most similar style to the function lm.
 ##        dblm.default: needs a dependent variable y and z's.
 ##        dblm.D2: response variable y and a 'D2' object.
 ##        dblm.dist: response varaible y and a dist or dissimilarity object. 
 ##    
 ##    Inputs:  metric: "euclidean", "gower" or "manhattan"; 
 ##             method: "OCV", "GCV", "AIC","BIC","eff.rank" and "rel.gvar";
 ##             weights: individuals's weights;
 ##             eff.rank: effective rank. Only if method=eff.rank;
 ##             rel.gvar: relative geometric variability (default 0.95)
 ##    Outputs: an object of class dblm
 ##   



    ###############################
    #### dblm of class formula ####
    ###############################

     #generic function with a commun parameter (y).
 dblm <- function(...)  UseMethod("dblm")
 
dblm.formula <- function(formula,data,...,metric="euclidean",method="OCV",
                    full.search=TRUE,weights,rel.gvar=0.95,eff.rank) 
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
  try(ans <- dblm.yz(y=zy$y,z=zy$z,metric=metric,weights=weights,
        eff.rank=eff.rank,method=method,rel.gvar=rel.gvar,full.search=full.search))  
  
  if (class(ans)=="try-error") 
    return(paste("the program failed. Try to read the help. If the error persists attempts to communicate with us "))
  
  # call dbglm
  ans$call <- mf
  attr(ans,"zs") <- zy$zini
  
  return(ans)   
}



    #############################
    #### aux dblm (y,z) #####
    #############################

dblm.yz <- function(y,z,metric="euclidean",method="OCV",full.search=TRUE,
                  weights,rel.gvar=0.95,eff.rank,...)
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
  try(ans<-dblm.D2(y=y,D2=D2,weights=weights,eff.rank=eff.rank,method=method,
               rel.gvar=rel.gvar,full.search=full.search)) 
   
   
  # y and Distance are defined--> pass to dist method (try for avoid crash). 
  #try(ans <- dblm.dist(y,D,weights=weights,eff.rank=eff.rank,method=method,
  #           rel.gvar=rel.gvar,full_search=full_search)) 
  if (class(ans)=="try-error") 
   return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
  
  ans$call <- call

  attr(ans,"metric") <- metric
  attr(ans,"zs") <- z
  attr(ans,"way") <- way        
  
  return(ans)
}

    #################################
    ####    dblm with Dist or   ####
    ####  dissimilarity distance ####
    #################################

dblm.dist <- function(distance,y,...,method="OCV",full.search=TRUE,weights,
                rel.gvar=0.95,eff.rank)
{   
   call <- match.call(expand.dots = FALSE)                   
   # program controls: distance must be of class D2 dist or dissimilarity.
   if ((missing(distance)||is.null(distance)))
    stop("distance must be defined")
   
   # dist to D2
   Delta <- disttoD2(distance)     
   
   try(ans<-dblm.D2(D2=Delta,y=y,weights=weights,eff.rank=eff.rank,
                  method=method,rel.gvar=rel.gvar,full.search=full.search))
 
   if (class(ans)=="try-error")
     return(paste("the program failed. Try to read the help. If the error persists attempts to communicate with us "))
   
   ans$call <- call
   attr(ans,"way")<-"D2"   

   return(ans)
}


    #################################
    ####  dblm with D2 distance ####
    #################################
    
dblm.D2 <- function(D2,y,...,method="OCV",full.search=TRUE,weights,
            rel.gvar=0.95,eff.rank)
{  
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
         
   try(ans<-dblm.Gram(G=G,y=y,weights=ori_weights,eff.rank=eff.rank,
                  method=method,rel.gvar=rel.gvar,full.search=full.search))
                  
 
   if (class(ans)=="try-error")
     return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
    
    ans$call<-match.call(expand.dots = FALSE)
    attr(ans,"way")<-"D2"
    
    return(ans)
}


    #############################
    ####    dblm with Gram   ####
    #############################
 
 
 dblm.Gram <- function(G,y,...,method="OCV",full.search=TRUE,weights,
              rel.gvar=0.95,eff.rank)
{    
  
   if (!is.Gram(G))
     stop("for a dblm.Gram method the class of the G matrix must be 'Gram'")
    
   # another program controls: See the auxiliar function
   controls <- controls_dblm(G,weights,eff.rank,rel.gvar,method,y)
   eff.rank <- controls$eff.rank
   ini_eff.rank <- eff.rank 
   weights <- controls$weights


   # not expained variability 1-rel.gvar (for method rel.gvar)
   epsilon<-1-rel.gvar
    
   # calculing weights, Dw, y0 and rel.gvar
   ori_weights <-weights           # originals weights !!
   weights <- weights/sum(weights) # percent weights !!  
   Dw   <- diag(weights)          # diagonal matrix with the weights
   Dsqw <- diag(sqrt(weights))
   y <- as.matrix(y)
   n <- nrow(y)
   y0 <- y - sum(weights*y)       # centered response varaible(y) 
   
   g <- diag(G)

   # important! definim el paramatre intern rk, té us sentit similar al eff.rank, 
   # l'unic és que facilita el càlcul del nou rank si s'usa un dels 4 metodes. 
   # mirar en les funcions internes el funcionament.
    
   threshold <- 0
   # Ordinary Cross- validation to choose the effective rank
   if (method=="OCV"){
    # the limit. Number of components such that at least, takes the 99% variability.   
    threshold <- HwProject(G,Dsqw,rk=eff.rank,epsilon=epsilon,cvyes=TRUE)$eff.rank 
  
    # find the optimal ocv.
    f <- function (rk,G, n, Dsqw, weights, epsilon,
        y, y0, cvyes, ori_weights,method)
    {
       rk<-round(rk)
       auxHwyhat <- Hwyhat(G, n, Dsqw, weights, rk, epsilon,
        y, y0, cvyes = TRUE, ori_weights = ori_weights)
       return(auxHwyhat$ocv)
    }
    
    if (!full.search){
     ocv_opt <- optimise(f = f,c(1,threshold),G=G, n=n, Dsqw=Dsqw, weights=weights,
      epsilon=epsilon, y=y, y0=y0, cvyes = TRUE, ori_weights=ori_weights,method=method,
      tol =1)
    
     eff.rank <- round(ocv_opt$minimum)
     ocv <- ocv_opt$objective
    }
    if (full.search){
     ocvs<-array(0,threshold)
     ocvs<- apply(as.matrix(1:threshold),1,function(i){f(rk=i,Dsqw=Dsqw,G=G,
              weights=weights,y=y,y0=y0,n=n,ori_weights=ori_weights,method=method)})
     ocv <-min(ocvs)
     eff.rank <-which.min(ocvs)
    }           
    
   }
   # Generalized Cross- validation to choose the effective rank
   if (method=="GCV"){

     threshold<- HwProject(G,Dsqw,rk=eff.rank,epsilon=epsilon,cvyes=TRUE)$eff.rank
      
     # find the optimal ocv.
     f <- function (rk,G, n, Dsqw, weights, epsilon,
        y, y0, cvyes, ori_weights,method)
     {
       rk<-round(rk)
       auxHwyhat <- Hwyhat(G, n, Dsqw, weights, rk, epsilon,
        y, y0, cvyes = TRUE, ori_weights = ori_weights)
       return(auxHwyhat$gcv)
     }
   
     if (!full.search){
      gcv_opt<-optimise(f = f,c(1,threshold),G=G, n=n, Dsqw=Dsqw, weights=weights,
       epsilon=epsilon, y=y, y0=y0, cvyes = TRUE, ori_weights=ori_weights,method=method,
       tol =1)
    
      eff.rank <- round(gcv_opt$minimum)
      gcv <- gcv_opt$objective
     }
     if (full.search){
      gcvs<-array(0,threshold)
      gcvs<- apply(as.matrix(1:threshold),1,function(i){f(rk=i,Dsqw=Dsqw,G=G,
              weights=weights,y=y,y0=y0,n=n,ori_weights=ori_weights,method=method)})
      gcv <-min(gcvs)
      eff.rank <-which.min(gcvs)
    }         
   } 
  
  
   # Aikaike and Bayesian criterium to choose the effective rank
   if (method=="AIC" || method=="BIC"){
     threshold <- HwProject(G,Dsqw,rk=eff.rank,epsilon=epsilon,cvyes=TRUE)$eff.rank
     
     # find the optimal ocv.
     f <- function (rk,G, n, Dsqw, weights, epsilon,
        y, y0, cvyes, ori_weights,method)
     {
       rk<-round(rk)
       auxHwyhat <- Hwyhat(G, n, Dsqw, weights, rk, epsilon,
        y, y0, cvyes = TRUE, ori_weights = ori_weights)
       rss<-(auxHwyhat$resStand.err)^2*auxHwyhat$rdf/n    
       if (method=="AIC") b_aic<-2*(auxHwyhat$eff.rank+1)+n*log(rss)  
       if (method=="BIC") b_aic<-n*log(rss)+(auxHwyhat$eff.rank+1)*log(n) 
       return(b_aic)
     }
     
     if (!full.search){
      aic_opt<-optimise(f = f,c(1,threshold),G=G, n=n, Dsqw=Dsqw, weights=weights,
       epsilon=epsilon, y=y, y0=y0, cvyes = TRUE, ori_weights=ori_weights,method=method,
       tol =1)
    
      eff.rank <- round(aic_opt$minimum)
      if (method=="AIC") aic <- aic_opt$objective
      if (method=="BIC") bic <- aic_opt$objective
     }
     if (full.search){
      aics <-array(0,threshold)
      aics<- apply(as.matrix(1:threshold),1,function(i){f(rk=i,Dsqw=Dsqw,G=G,
              weights=weights,y=y,y0=y0,n=n,ori_weights=ori_weights,method=method)})
      if (method=="AIC") aic <- min(aics)
      if (method=="BIC") bic <- min(aics)
      eff.rank <-which.min(aics)
     }         
   }      
          
   # call to hwyhat function to calculate the hat matrix Hw, with the 
   # apropiate eff.rank. Also gives the ocv and  gcv estimators.
    if (method!="rel.gvar"){
      hwyhat<-Hwyhat(G,n,Dsqw,weights,eff.rank,epsilon,y,y0,
                ori_weights=ori_weights)  
    }else      
       hwyhat<-Hwyhat(G,n,Dsqw,weights,rk=0,epsilon,y,y0,
                ori_weights=ori_weights) 
               
   Hw<- hwyhat$Hw
   Fwplus<- hwyhat$Fwplus
   eff.rank<- hwyhat$eff.rank
   ocv<- hwyhat$ocv
   gcv<- hwyhat$gcv
   used_rel.gvar<-hwyhat$used_rel.gvar
    
   # fitted.values (fitted values), R2, adjusted R2 and the residual standard error
   fitted.values <- hwyhat$yhat 
   resStand.err<-hwyhat$resStand.err
   
   rss<-resStand.err^2*hwyhat$rdf/n  
   aic <- 2*eff.rank+n*log(rss) 
   bic <- n*log(rss)+eff.rank*log(n)  


   call<- match.call(expand.dots = FALSE)

     
   # return a list with the ocv and gcv estimator, the used eff.rank, the fitted
   # values yhat, R2 and adjusted R2, the residual standard error,
   # residual degrees of freedom (rdf), residuals(y-fitted values), 
   # the hat matrix Hhat, ocvs, gcvs, aics and bics with the esitmators of 
   # ocv,... for all the eff.ranks down of the 99% of variability. aic and bic 
   # (if one of this two method are used). And the last one, the original 
   # weights (not percents)

   ans<-list(residuals=(y-fitted.values),fitted.values=fitted.values,
             df.residuals=(n-1-eff.rank),weights=ori_weights,y=y,H=Hw,
             call=call,rel.gvar=used_rel.gvar,eff.rank=eff.rank,ocv=ocv,
             gcv=gcv,aic=aic,bic=bic)
   
   
   # attributes generaly not call by the user(for plot,print and predict functions)
   attr(ans,"sigma")<-resStand.err
   attr(ans,"full.search") <- full.search
   attr(ans,"G")<-G
   attr(ans,"method")<-method
   attr(ans,"Fwplus")<-Fwplus
   attr(ans,"way")<-"G"   
   attr(ans,"threshold")<- threshold   
   
   if (method=="rel.gvar") attr(ans,"ini_rel.gvar")<- rel.gvar
   if (method=="eff.rank") attr(ans,"ini_eff.rank")<- ini_eff.rank
   if (method=="OCV"&&full.search) attr(ans,"ocvs")<- ocvs
   if (method=="GCV"&&full.search) attr(ans,"gcvs")<- gcvs
   if (method=="AIC"&&full.search) attr(ans,"aics")<- aics
   if (method=="BIC"&&full.search) attr(ans,"bics")<- aics
       
   class(ans)<-"dblm"   
   return(ans) 
}

