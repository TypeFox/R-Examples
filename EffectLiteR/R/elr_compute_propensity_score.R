
computePropensityScore <- function(input){
  
  propscore <- input@vnames$propscore
  
  ## propensity score
  if(!is.null(propscore)){
    
    x <- input@vnames$x
    d <- input@data
    ng <- input@ng
    
    if(is(propscore, "formula")){      
      form <- propscore
      environment(form) <- environment()
      
    }else{
      form <- as.formula(paste0(x, " ~ ", paste0(propscore, collapse=" + ")))
    }
    
    mprop <- nnet::multinom(form, data=d, na.action="na.omit", trace=FALSE)
    
    ## save output
    resprop <- summary(mprop)
    outprop <- list()
    outprop$formula <- paste0(deparse(form))
    outprop$coef <- resprop$coefficients
    outprop$se <- resprop$standard.errors
    outprop$tval <- resprop$coefficients/resprop$standard.errors
    input@outprop <- outprop
    
    ## fitted values
    dprop <- fitted(mprop)
    if(input@ng > 2){dprop <- dprop[,-1]}
    dprop <- apply(dprop,2,car::logit)       
    
    if(any(diag(var(dprop)) < 0.05)){
      warning(paste("very small variance of propensity scores \n ",
                    diag(var(dprop))))
    }
    
    dprop <- dprop[match(row.names(d), row.names(dprop)),] ## for missings    
    dprop <- as.data.frame(dprop) 
    names(dprop) <- paste0("logprop",1:(ng-1))
    input@data <- cbind(d,dprop)
    input@vnames$z <- c(input@vnames$z,paste0("logprop",1:(input@ng-1)))
    input@nz <- length(input@vnames$z)
    
  }
  
  return(input)
  
}


