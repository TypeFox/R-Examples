predict.sprmda <-
  function(object, newdata, ...)
    ## 14/10/09 IH
    ## 
    ## predicts class membership based on a Sparse Partial Robust M regression model of class sprmda
    ##
    ## uses ldafitfun
    ##
    ## Inputs: 
    ## object .... a "sprmda" class Partial Robust M regression object 
    ## newdata ... an optional data matrix or frame containing a new set of cases 
    ##
  {    
    
    if(!(class(object)=="sprmda")){stop("The SPRM predict function only applies to sprmda class objects")}
    b <- coef(object)
    if(missing(newdata)){
      Xn <- object$input$X0
    } else {
      if(ncol(newdata)==length(b)){
        Xn <- newdata
      } else {
        formula <- object$input$formula
        mt <- terms(formula, data=object$input$X0)
        ic <- attr(mt, "intercept")
        if (ic==0){
          Xn <- model.matrix(mt, newdata)
        } else{
          Xn <- model.matrix(mt, newdata)[,-1]
        } 
      }  
    }
    
    scores <- scale(Xn,center=object$XMeans,scale=object$Xscales) %*% object$R
    
    ldafit <- apply(scores, 1, ldafitfun, object$ldamod$cov, object$ldamod$m1, object$ldamod$m2, object$input$prior)
    pred <- (apply(ldafit, 2, which.max)-1.5)*(-2)  
    
    return(pred)
  }
