#' Prediction Based on the Fitted Additive Hazards Model 
#'
#' This function predicts a subject's overall hazard rates at given time points based on this subject's covariate 
#' values.  The prediction function is an additive hazards model fitted using \code{\link{ah}}. 
#'  
#'
#' @param object an object of class inhering from \code{\link{ah}}.
#' @param newdata  a dataframe of an individual's predictors. 
#' @param newtime  a given sequence of time points at which the prediction is performed. 
#'  The time should be on the same scale as the survival time in \code{\link[survival]{Surv}}.
#' @param ...  further arguments passed to or from other methods.
#'
#' @return A dataframe including the time points for prediction, predicted values and their standard errors.
#'
#' @seealso \code{\link{ah}} for fitting the additive hazards model, \code{\link{nwtsco}} for 
#' the description of nwtsco dataset
#'
#' @importFrom survival Surv
#' @importFrom stats delete.response model.matrix terms
#' @export
#'
#' 
#'
#' 
#' @examples
#' library(survival)
#' ###  fit the additive hazards model to the data 
#' nwts<- nwtsco[1:100,]
#' fit <- ah(Surv(trel,relaps) ~ age + instit, data = nwts, robust = FALSE)
#' 
#' ###  see the covariate names in the prediction function 
#' fit$call
#' ###  the newdata should be a dataframe
#' ###  the variable names are the same as the covariate names in 
#' ###  the prediction function
#' newdata <- data.frame(age=60, instit =1) 
#'
#' ###  an alternative way to give the newdata 
#' newdata <- nwtsco[101,] 
#' 
#' ###  based on this subject's covariate values, the function predicts  individual specific
#' ###  hazard rates at time points 3 and 5 
#' predict(fit, newdata, newtime = c(3,5))
#' 








predict.ah <- function(object, newdata, newtime,...){

  	# get all the middle step results obtained in the model fitting step 
    ingr <-ingredients(object)
    

    #names(ingredients.list)
    #[1] "Z"           "t.fail"      "censor"      "n.obs"       "wts"        
    #[6] "t.unique"    "t.diff"      "t1.unique"   "n.par"       "theta"      
    #[11] "effects"     "eta.ncol"    "match.eta"   "eta"         "eta.den.vec"
    #[16] "eta.cum"     "idx.c"       "match.event" "n.death"     "z.death"    
    #[21] "dLambda1"    "lambda0"     "Lambda0" 
   
   tt <- terms(object)
  if (!inherits(object, "ah")) 
    warning("calling predict.ah(<fake-ah-object>) ...")
  #if (missing(newdata) || is.null(newdata)){
  # pred <- object$data$X  #predictor matrix
  #  pred <- as.matrix(pred) 
  # }else{
    ## retrieve the columns in newdata that match the names of the predictors in the model 
  

  
     Terms <- delete.response(tt)
    # identify factor variable
    factor.var.list<- names(object$model)[sapply(object$model,is.factor) == T]
    #assign levels to newdata variable according to the original data 
    if (length(factor.var.list) != 0){
      for (i in 1:length(factor.var.list)){
        factor.var <- factor.var.list[i]
         newdata[[factor.var]] <- factor(newdata[[factor.var]], 
                                         levels = levels(object$model[[factor.var]]))
         }
    }
     
      pred <- model.matrix(Terms, newdata)
     #contrasts.arg = lapply(newdata[sapply(newdata, is.factor) ], constrats, , contrasts=FALSE
      pred <- as.numeric(pred[,-1])
     
    # predicted additive effect 
    effect <- sum(pred * object$coef)
     

    #  middle step results for obtaining predicted outcomes and variances
 	  #b <- do.call("cook.basics",ingredients.list)
    

    #cumulated  baseline hazards
    #Lambda0 <- b$Lambda0
   
    # search the location of the newtime in the original timelines in the fitted model
    t.pos <- NULL
    t.unique<-ingr$t.unique
  
    # Store the predicted value
    L <- NULL
    
  	for (i in 1:length(newtime)){
    # newtime s is less than t_1, then we report  
      if (newtime[i] < t.unique[1]){
        print(paste("The data used for building the ah model does not have enough information for predicting such a small t, t=", newtime[i],".") )
        print(paste("remove, t=", newtime[i]," and run the prediction again") )
      }else{
      
        # find the position of newtime[i] on the original timeline
        c <- rank(c(newtime[i],t.unique),ties.method="max")[1] 
        #when newtime  s equals one of the t.unique t_k, the above rank function returns k+1 on the original timelines
        #when newtime  s in between t_k and t_{k+1}, it returns k+1 on the original timelines
        #when newtime  s = t_1, it returns 1
        # thus 
        t.pos[i] <- ifelse(c==1,1,c-1)   
        # The predicted valued
        # the additive part I 
        L[i] <- effect*newtime[i] + ingr$Lambda0[t.pos[i]]
        
    	}
    }
	
    L.var <-NULL
    if (!object$robust){
       # variance 
      # L.var<-do.call("L.nvar", c(ingredients.list, list(n.death=n.death, z.death=z.death, eta.den.vec=eta.den.vec, eta=eta,
       # dLambda1=dLambda1, eta.cum=eta.cum, iA=object$iA, B=object$var, time.pos=t.pos, bhaz=FALSE, newtime=newtime)))
      L.var<-do.call("L.nvar", c(ingr, object,  list(t.pos= t.pos, pred = pred, newtime = newtime)))

    
      }else{  
       L.resid <-do.call("cook.L.resid", c(ingr, list(time.pos = ingr$match.event)))

       # robust variance
       L.var<-do.call("L.rvar", c(ingr, object, list(t.pos =t.pos, pred =pred, newtime=newtime, L.resid = L.resid)))
    }
    L.se <- sqrt(L.var)
  	return(data.frame(time = newtime, L= L, L.se = L.se))

}







