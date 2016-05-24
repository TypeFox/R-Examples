

 ########################
 #### ldblm function ####
 ########################
 
 ## Description:
 ##     local linear distance-based Generalized Regression. Fit a polynomial 
 ##     surface determined by one or more numerical predictors, using local 
 ##     distance fitting. Cut short the range of the entry variables n times to 
 ##     approximate at each individual a linear model, with which estimate each 
 ##     Yi prediction. Modern regression methods are designed to address 
 ##     situations in which the classical procedures do not perform well or 
 ##     cannot be effectively applied without undue labor.
 ##        
 


    ################################
    #### dbglm of class formula ####
    ################################
 ldblm<-function(...)  UseMethod("ldblm")
     
ldblm.formula<-function(formula,data,...,kind.of.kernel=1,
              metric1="euclidean",metric2=metric1,method.h="GCV",weights,
              user.h=NULL,h.range=NULL,noh=10,k.knn=3,rel.gvar=0.95,
              eff.rank=NULL)
{ 
  # call dbglm
  mf <- match.call(expand.dots = FALSE)
  
  # control metric. See the auxiliar function
  metric1 <- control_metric(metric1)
  metric2 <- control_metric(metric2)  
  
  if (missing(data))
    data <- environment(formula)
  
  # recover z and y of the formula                            
  if (metric1=="gower"||metric2=="gower")
    zy <- formula_to_zy(formula,data,mf,"dblm","gower")
  else
    zy <- formula_to_zy(formula,data,mf,"dblm",metric1)
    
  # y and z are defined--> pass to default method (try for avoid the program crash). 
  try(ans<-ldblm.yz(y=zy$y,z=zy$z,kind.of.kernel=kind.of.kernel,
        method.h=method.h,weights=weights,metric1=metric1,metric2=metric2,
        user.h=user.h,h.range=h.range,noh=noh,k.knn=k.knn,rel.gvar=rel.gvar,
        eff.rank=eff.rank)) 
   
  if (class(ans)=="try-error") 
     return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
 
  # call dbglm
  ans$call<-mf
  attr(ans,"zs") <- zy$zini
 
  return(ans)  
}
      


    ################################
    #### default dbglm (y,z) #######
    ################################

ldblm.yz <- function(y,z,kind.of.kernel=1,metric1="euclidean",
        metric2=metric1,method.h="GCV",weights,user.h=NULL,h.range=NULL,
        noh=10,k.knn=3,rel.gvar=0.95,eff.rank=NULL,...)
{
  
   # See if z or distance matrix is defined by the user.
   # require(cluster)
  
   # control metric. See the auxiliar function
   metric1 <- control_metric(metric1)
   metric2 <- control_metric(metric2)
  
   # call z_to_dist to pass the explanatory variables to an object of class dist 
   dist_and_way <- z_to_dist(z,metric1)
   dist1 <- dist_and_way$D
   dist2 <-  z_to_dist(z,metric2)$D
   way <- dist_and_way$way
  
   # if metric=gower. the distance matrix D is already the squared.
   if (metric1=="gower")
    D2.1 <-as.matrix(dist1)
   else
    D2.1 <-as.matrix(dist1)^2 
   if (metric2=="gower")
    D2.2 <-as.matrix(dist2)
   else
    D2.2 <-as.matrix(dist2)^2
     
   class(D2.1) <- "D2"
   class(D2.2) <- "D2"
   
   try(ans <- ldblm.D2(D2.1=D2.1,D2.2=D2.2,y=y,kind.of.kernel=kind.of.kernel,
     method.h=method.h,weights=weights,user.h=user.h,h.range=h.range,noh=noh,
     k.knn=k.knn,rel.gvar=rel.gvar,eff.rank=eff.rank)) 
  
   if (class(ans)=="try-error") 
    return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
  
   ans$call<-match.call(expand.dots = FALSE)
  
   # hidden attributes
   if (!missing(metric1))
    attr(ans,"metric1") <- metric1
   if (!missing(metric2))
    attr(ans,"metric2") <- metric2
   
   attr(ans,"zs")<-z
   attr(ans,"way")<- way
  
   # return ans
   return(ans)
}
                      

    #################################
    ####    dbglm with Dist or   ####
    ####  dissimilarity distance ####
    #################################

ldblm.dist <- function(dist1,dist2=dist1,y,kind.of.kernel=1,method.h="GCV",
           weights,user.h=quantile(dist1,.25),
           h.range=quantile(as.matrix(dist1),c(.05,0.5)),noh=10,k.knn=3,
           rel.gvar=0.95,eff.rank=NULL,...){
          
   # stop if class of distance matrix is not dist
   if (!any (class(dist1)=="dist")) 
    stop("for a ldblm.dist method the class of the distance matrix dist1 must be 'dist'")
   if (!any (class(dist2)=="dist")) 
    stop("for a ldblm.dist method the class of the distance matrix dist2 must be 'dist'")
   
   # dist to D2
   Delta1 <- disttoD2(dist1)     
   Delta2 <- disttoD2(dist2)     

      
   # y and Distance are defined--> pass to dist method (try for avoid the program crash). 
   try(ans <- ldblm.D2(D2.1=Delta1,D2.2=Delta2,y=y,kind.of.kernel=kind.of.kernel,
              method.h=method.h,weights=weights,user.h=user.h,h.range=h.range,
              noh=noh,k.knn=k.knn,rel.gvar=rel.gvar,eff.rank=eff.rank))    
   if (class(ans)=="try-error")
    return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
   
   ans$call <- match.call(expand.dots = FALSE) 
   # hidden attributes  
   attr(ans,"way") <- "D2"  
   return(ans)
   

 }
  
  

    #################################
    ####  ldblm with D2 distance ####
    #################################

ldblm.D2<-function(D2.1,D2.2=D2.1,y,kind.of.kernel=1,method.h="GCV",weights,
         user.h=quantile(D2.1,.25)^.5,h.range=quantile(as.matrix(D2.1),c(.05,0.5))^.5,
         noh=10,k.knn=3,rel.gvar=0.95,eff.rank=NULL,...){
   
    # control method. See the auxiliar function
    method.h <- control_method(method.h,"ldblm")  
   
    # another controls: see the auxiliar function
    controls <- controls_ldblm(D2.1,D2.2,user.h,method.h,h.range,noh,k.knn,
            kind.of.kernel,y,weights)
    user.h <- controls$user.h
    h.range <- controls$h.range
    weights <- controls$weights
    ori_weights <- weights
    weights <- weights/sum(weights)
     
   # sequence of bandwidth to be evaluate
   if (method.h!="user.h"){
     h_low <- h.range[1]
     h_up <- h.range[2]
     h_vec <- exp(seq(log(h_low),log(h_up),length=noh))
   }
   n <- length(y) 

   # k.knn: three nearest neigbourh
   h.knn<-h.knn.funct(D2.1^.5,k=k.knn)
   
   # compute the model for each method.h
   if (method.h!="user.h"){

    if (method.h=="OCV") OCV <- rep(0,noh) # OCV's for each h(only if method==OCV)
    if (method.h=="GCV") GCV <- rep(0,noh) # GCV's for each h(only if method==GCV)
    if (method.h=="AIC") AIC <- rep(0,noh) # AIC's for each h(only if method==AIC)
    if (method.h=="BIC") BIC <- rep(0,noh) # BIC's for each h(only if method==BIC)
    
    i <- 0
    for (h in h_vec){
      # fitted values and Shat for each bandwidth (h) 
      aux <-pred.train.sample(y,D2.1,D2.2,n,h,h.knn,kind.of.kernel,ori_weights,
                    rel.gvar,eff.rank) 
      i <- i+1
      S <- aux$Shat
      fitted.values <- aux$yhat
      diagS <- diag(S)
      nu <- sum(diagS)
      
      # Ordinary cross validation criterium to choose the best bandwidth.
      if (method.h=="OCV"){
        OCV[i] <-  sum(weights*((y-fitted.values)/(1-diagS))^2 ) # ocv formula
         if (is.nan(OCV[i]))
          stop(paste("OCV for the bandwidth ",round(h,5), " is a NaN. Try to use another method.h or h.range"))
         
         if (i==1){
         OCV_opt <- OCV[i]
         h.opt <- h
         S_opt <- S
         yhat_opt<-fitted.values
        }else{
         if (OCV_opt > OCV[i]){ # improve the best ocv
            OCV_opt <- OCV[i]   # edit the ocv optim
            h.opt <- h          
            S_opt <- S
            yhat_opt<-fitted.values
         }
        }
     }else OCV_opt<-NULL
        
     # Generalized cross validation criterium to choose the best bandwidth.
      if (method.h=="GCV"){
      
        GCV[i] <-sum(weights*(fitted.values-y)^2)/(n*(1-nu/n)^2)   # gcv formula
        if (is.nan(GCV[i]))
          stop(paste("GCV for the bandwidth ",round(h,5), " is a NaN. Try to use another method.h or h.range"))
        if (i==1){
         GCV_opt <- GCV[i]
         h.opt <- h
         S_opt <- S
         yhat_opt<-fitted.values
        }else{
         if (GCV_opt > GCV[i]){   # improve the best gcv
            GCV_opt <- GCV[i]     # edit the gcv optim
            h.opt <- h
            S_opt <- S
            yhat_opt<-fitted.values
         }
        }
     }else GCV_opt<-NULL

     # Aikaike criterium to choose the best bandwidth.
      if (method.h=="AIC"){
        rss<-sum(ori_weights*(fitted.values-y)^2)/n   # residual standard desviation
        AIC[i]<-2*sum(diag(S))+n*log(rss)   # aic formula

        if (i==1){
         AIC_opt <- AIC[i]
         h.opt <- h
         S_opt <- S
         yhat_opt<-fitted.values
        }else{
         if (AIC_opt > AIC[i]){  # improve the best aic
            AIC_opt <- AIC[i]    # edit the aic optim
            h.opt <- h
            S_opt <- S
            yhat_opt<-fitted.values
         }
        }
     }else AIC_opt<-NULL

     # Bayesian criterium to choose the best bandwidth.
      if (method.h=="BIC"){
        rss<-sum(ori_weights*(fitted.values-y)^2)/n       # residual standard desviation
        BIC[i]<-n*log(rss)+log(n)*sum(diag(S))  # bic formula

        if (i==1){
         BIC_opt <- BIC[i]
         h.opt <- h
         S_opt <- S
         yhat_opt<-fitted.values
        }else{
         if (BIC_opt > BIC[i]){      # improve the best bic
            BIC_opt <-BIC[i]         # edit the bic optim
            h.opt <- h
            S_opt <- S
            yhat_opt<-fitted.values
         }
        }
       }else BIC_opt<-NULL
     }
    }
    
    # if method.h=user.h --> the model is estimated with the h bandwidth defined by the user
    if (method.h=="user.h"){
      # fitted values and Shat for user.h bandwidth
      aux <-pred.train.sample(y,D2.1,D2.2,n, h=user.h, h.knn,kind.of.kernel,ori_weights,
              rel.gvar,eff.rank)    
      S <- aux$Shat
      yhat_opt <- aux$yhat
      h.opt<-user.h
      
      OCV_opt<-NULL
      GCV_opt<-NULL
      AIC_opt<-NULL
      BIC_opt<-NULL
    }


    call<- match.call(expand.dots = FALSE)

    # return the next attributes 
    ans<-list(residuals=y-yhat_opt,fitted.values=yhat_opt,h.opt=h.opt,S=S,
              y=y,weights=weights,call=call,dist1=D2.1,dist2=D2.2) 
    
    attr(ans,"kind.of.kernel")<-kind.of.kernel             
    attr(ans,"method.h")<-method.h
    attr(ans,"dist1")<-D2.1
    attr(ans,"dist2")<-D2.2
    attr(ans,"OCV_opt")<-OCV_opt
    attr(ans,"GCV_opt")<-GCV_opt
    attr(ans,"AIC_opt")<-AIC_opt
    attr(ans,"BIC_opt")<-BIC_opt   
    attr(ans,"noh")<-noh  
    attr(ans,"way")<-"D2"   
    attr(ans,"rel.gvar")<-rel.gvar
    attr(ans,"eff.rank")<-eff.rank
    
    if (method.h!="user.h") attr(ans,"h_vec")<-h_vec  
    if (method.h=="OCV") attr(ans,"OCV")<-OCV
    if (method.h=="GCV") attr(ans,"GCV")<-GCV
    if (method.h=="AIC") attr(ans,"AIC")<-AIC
    if (method.h=="BIC") attr(ans,"BIC")<-BIC
   
    class(ans)<-"ldblm"
    return(ans)
    
}
     
    ##########################
    ####  ldblm with Gram ####
    ##########################

ldblm.Gram <- function(G1,G2=G1,y,kind.of.kernel=1,method.h="GCV",weights,
         user.h=NULL,h.range=NULL,noh=10,k.knn=3,rel.gvar=0.95,
         eff.rank=NULL,...){
    
   # stop if class of distance matrix is not D2      
   if (class(G1)[1]!="Gram")
    stop("for a ldblm.Gram method the class of the distance matrix G1 must be 'Gram'")
   if (class(G2)[1]!="Gram")
    stop("for a ldblm.Gram method the class of the distance matrix G2 must be 'Gram'")   
   
   # converts G to D2
   D2.1 <- GtoD2(G1)
   D2.2 <- GtoD2(G2) 
    
   # y and Distance are defined--> pass to dist method (try for avoid the program crash). 
   try(ans <- ldblm.D2(D2.1=D2.1,D2.2=D2.2,y=y,kind.of.kernel=kind.of.kernel,
              method.h=method.h,user.h=user.h,h.range=h.range,noh=noh,k.knn=k.knn,
              rel.gvar=rel.gvar,eff.rank=eff.rank))    
   if (class(ans)=="try-error")
    return(paste("the program failed.Tries to read the help. If the error persists attempts to communicate with us "))
   
   ans$call <- match.call(expand.dots = FALSE) 
   # hidden attributes  
   attr(ans,"way") <- "G" 
   attr(ans,"G1") <- G1 
   attr(ans,"G2") <- G2 
    
   return(ans)
}

