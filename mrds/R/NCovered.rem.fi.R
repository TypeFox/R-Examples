# Compute abundance in covered region for rem.fi model object
#
# par    - parameter values (used when computing derivatives wrt parameter uncertainty)
# model  - ddf model object
# group  - if TRUE computes group abundance and if FALSE individual abundance
#
# result - abundance estimate
NCovered.rem.fi <- function(par=NULL,model,group=TRUE,...){
  if(!is.null(par)){
    model$mr$coefficients <- par
    fitted <- predict(model,compute=TRUE,integrate=TRUE)$fitted
  }else{
    fitted <- model$fitted
  }

  if(!group){
    size <- model$data$size[model$data$observer==1 &
                            model$data$object %in%
                              as.numeric(names(model$fitted))]
    Nhat <- sum(compute.Nht(fitted,FALSE,size))
  }else{
    Nhat <- sum(compute.Nht(fitted,TRUE,size=NULL))
  }

  return(Nhat)
}
