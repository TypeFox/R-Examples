#' Prediction Based on the Additive Hazards Model Fitted from Two-phase Sampling 
#' 
#' This function predicts a subject's overall hazard rates at given time points based on this subject's covariate
#' values. The prediction function is an object from \code{\link{ah.2ph}}. The  estimating procedures follow Hu (2014). 
#'  
#' @param object an object of class inhering from "ah.2ph".
#' @param newdata a dataframe of an individual's predictors. 
#' @param newtime  a given sequence of time points at which the prediction is performed. 
#' @param ...  further arguments passed to or from other methods.
#' 
#' @return A dataframe including the given time points, predicted hazards, their standard errors,
#'      their variances, the phase I component of the variance for predicted hazards
#'      and the phase II component of the variance.
#'
#' @seealso \code{\link{ah.2ph}} for fitting the additive hazards model with two-phase sampling and 
#' \code{\link{nwtsco}} for the description of nwtsco dataset
#'
#' @importFrom survival Surv
#' @importFrom stats delete.response model.matrix terms
#' @export
#'
#' @references
#' Jie Hu (2014) A Z-estimation System for Two-phase Sampling with Applications to Additive Hazards Models and 
#' Epidemiologic Studies. Dissertation, University of Washington.
#
#' @examples
#' library(survival)
#' ### load data
#' nwts <- nwtsco[1:100,]
#'
#' ### create strata based on  institutional histology and disease status
#' nwts$strt <- 1+nwts$instit
#' ### add a stratum containing all (relapsed) cases
#' nwts$strt[nwts$relaps==1] <- 3
#'
#' ### assign phase II subsampling probabilities
#' ### oversample unfavorable histology (instit =1) and cases
#' ### Pi = 0.5 for instit =0, Pi =1 for instit =1 and relaps =1 
#' nwts$Pi<-  0.5 * (nwts$strt == 1) + 1 * (nwts$strt == 2) + 1 * (nwts$strt == 3)
#'
#' ### generate phase II sampling indicators
#' N <- dim(nwts)[1]
#' nwts$in.ph2 <-  rbinom(N, 1, nwts$Pi)
#'
#' ### fit an additive hazards model to  two-phase sampling data without calibration
#' fit1 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts, R = in.ph2, Pi = Pi, robust = FALSE)
#'                                                            
#' 
#' 
#' ###  input the new data for prediction
#' newdata <- nwtsco[101,] 
#' ###  based on the fitted model fit1, perform prediction at time points t =3 and t= 5
#' predict(fit1, newdata, newtime = c(3,5))
#'
#' ### fit an additve hazards model to  two-phase sampling data with calibration
#' ### The calibration variable is stage
#' fit2 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts, R = in.ph2, Pi = Pi, 
#'                                    robust = FALSE, calibration.variables = stage)
#' 
#' ### based on the fitted model fit2, perform prediction at time points t =3 and t= 5
#' predict(fit2, newdata, newtime = c(3,5))





predict.ah.2ph<-function(object, newdata, newtime,...){

  # phase I complete data object
  object1 <-object$fit.pha1

  ######################### creat model matrix for the new data ##########
  tt <- terms(object1)
  if (!inherits(object, "ah.2ph")){
    warning("calling predict.2ph(<fake-ah-object>) ...")
  }
  
     Terms <- delete.response(tt)
    # identify factor variable
    factor.var.list<- names(object1$model)[sapply(object1$model,is.factor) == T]
    #assign levels to newdata variable according to the original data 
    if (length(factor.var.list) != 0){
      for (i in 1:length(factor.var.list)){
        factor.var <- factor.var.list[i]
        newdata[[factor.var]] <- factor(newdata[[factor.var]], 
                                           levels = levels(object1$model[[factor.var]]))
      }
    }
    # retrieve the new model matrix based on the new data
      pred <- model.matrix(Terms, newdata)
    # delete the intercept
      pred <- as.numeric(pred[,-1])
   ###################################################################################
     
    
  
    
   


    
    
   

  Pi.pha2<-object$Pi.pha2
  ## No calibration 
  if (object1$robust){
      print("If the additive hazards model does not hold, then stop prediction based on this model. No robust errors with be provided ")
   }else{
      if(!length(object$calibration.variables)){


         ingr<-ingredients(object1) 

             ############################ find the positions of a list of given newtime on the original timeline######################
    t.pos <- NULL
    t.unique<-ingr$t.unique
    
    for (i in 1:length(newtime)){
    # newtime s is less than t_1, then we report  
      if (newtime[i] < t.unique[1]){
           print(paste("The data used for building the ah model does not have enough information for predicting such a small t, t=", newtime[i]) )
      }else{
      
        # find the position of newtime[i] on the original timeline
        c <- rank(c(newtime[i],t.unique),ties.method="max")[1] 
        #when newtime  s equals one of the t.unique t_k, the above rank function returns k+1 on the original timelines
        #when newtime  s in between t_k and t_{k+1}, it returns k+1 on the original timelines
        #when newtime  s = t_1, it returns 1
        # thus 
        t.pos[i] <- ifelse(c==1,1,c-1)   
       
        
      }
    }

         ### predicted outcomes are calculated based on fitting an ah model to only phase II data using assigned weights
         pred.result <- predict.ah(object1, newdata, newtime, level = 0.95)
         L <- pred.result$L
         L.var.pha1 <- pred.result$L.se^2
         

         # Calculate the  phase II variance 
          L.resid <- do.call("cook.L.resid", c(ingr, list(time.pos = ingr$match.event)))



          
          ingr$wts <- sqrt(1-Pi.pha2)
          L.var.pha2<-do.call("L.rvar", c(ingr, object1, list(t.pos =t.pos, pred =pred, newtime=newtime, L.resid = L.resid)))
         
      }else{
      #Calibration 
        calibration.variables<-object$calibration.variables
        
        wts.pha2<-as.numeric(1/Pi.pha2)
        aux<-as.matrix(calibration.variables)
        P<-t(aux)%*%(aux)
        aux.pha2<-aux[object$R==1,]
        wts.cal<-cook.wts.cal(aux=aux,aux.pha2=aux.pha2,P=P,wts.pha2=wts.pha2)
        wts.pha2<-wts.cal
        
        
         #using the calibrated refit the additive hazards model 
         object1<-ah(object1$formula, object1$data, weights=wts.pha2, robust=object1$robust)
             

      
        ingr<-ingredients(object1) 

        t.pos <- NULL
        t.unique<-ingr$t.unique
    
    for (i in 1:length(newtime)){
    # newtime s is less than t_1, then we report  
      if (newtime[i] < t.unique[1]){
           print(paste("The data used for building the ah model does not have enough information for predicting such a small t, t=", newtime[i]) )
      }else{
      
        # find the position of newtime[i] on the original timeline
        c <- rank(c(newtime[i],t.unique),ties.method="max")[1] 
        #when newtime  s equals one of the t.unique t_k, the above rank function returns k+1 on the original timelines
        #when newtime  s in between t_k and t_{k+1}, it returns k+1 on the original timelines
        #when newtime  s = t_1, it returns 1
        # thus 
        t.pos[i] <- ifelse(c==1,1,c-1)   
       
        
      }
    }

         ### predicted outcomes are calculated based on fitting an ah model to only phase II data using  calibrated weights
        pred.result <- predict.ah(object1, newdata, newtime, level = 0.95)
        L <- pred.result$L

        #####################Calculate the variance of the predicted value 
        #See page 109 of Jie Hu' thesis for the detailed formula 

        ### Calculate the 1st component of the variance 
        L.var.pha1 <- pred.result$L.se^2
        
        ## retrieve the martingale integral, i.e.  part of \psi in the formula on page 109 of Jie Hu' thesis

        L.resid <- do.call("cook.L.resid", c(ingr, list(time.pos = ingr$match.event)))

         L.var.pha2<-do.call("L.rvar.calibration", c(ingr, object1, list(t.pos =t.pos, pred =pred, newtime=newtime, L.resid = L.resid,
                                                                        new.wts=sqrt((1-Pi.pha2)/(wts.cal*Pi.pha2)), aux.pha2=aux.pha2, P=P, wts.cal=wts.cal)))
      }
    }




     #L.var.pha2<-L.rvar(match.event=match.event,new.wts=sqrt(1-Pi.pha2),L.ncol=L.ncol, eta.cum=eta.cum, resid=resid, 
      #L.resid=L.resid,iA=iA, B=B)
     
     L.var=L.var.pha1+L.var.pha2
  
     L.se = sqrt(L.var)

    


      return(data.frame(L=L, L.se=L.se, L.var=L.var,L.var.pha1=L.var.pha1, L.var.pha2=L.var.pha2))
}
